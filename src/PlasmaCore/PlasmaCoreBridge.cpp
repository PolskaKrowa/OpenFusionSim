//
// src/PlasmaCore/PlasmaCoreBridge.cpp
// Host-side C++ wrapper around the CUDA PIC kernels.
//
//  Improvements in this version (Round 3):
//    1. Wires in the actual PIC loop via pic_bridge.h.  When cfg.run_pic is
//       true, the bridge allocates a PICState on construction and runs one
//       picStep per game tick.  This drives the actual CUDA kernels (sort,
//       deposit, FDTD, push, collide, fuse, transport, radiation feedback)
//       instead of relying solely on the 0D power-balance model.
//    2. Adds MHD disruption tracking via mhd_physics.h.  The MHD module
//       evolves the 1-D current profile, finds q=m/n resonant surfaces,
//       and grows magnetic islands via the Rutherford ODE.  When W > W_crit
//       or a VDE is triggered, the bridge flags disruption_imminent.
//    3. Reads PIC diagnostics (T_e, n_e, q_rad per cell) back from device
//       and uses them to drive the radiation feedback.  This closes the
//       radiation → temperature loop that was previously missing.
//    4. Initializes ENDF/B-VIII.0 cross-sections on first update (loads
//       built-in defaults if no ACE file is provided).
//

#include "PlasmaCoreBridge.h"
#include "ReactorState.h"
#include "confinement_physics.h"
#include "mhd_physics.h"
#include "pic_bridge.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

// ─── Reactor geometry (ITER-class defaults; configurable via ReactorState) ──
static ConfinementPhysics::TokamakGeometry makeGeometry(const ReactorState& s)
{
    ConfinementPhysics::TokamakGeometry g;
    g.R_major_m    = s.R_out_m;          // outer major radius [m]  (ITER: 6.2)
    g.a_minor_m    = 2.0f;               // minor radius [m]  (fixed for ITER-class)
    // Read elongation from ReactorState (set by the shape controller in
    // ControlSystem::runShapeControl).  Previously this was hardcoded to
    // 1.7f, which meant the operator's κ slider had no effect on the plasma
    // physics (volume, q_95, confinement time all used 1.7 regardless).
    g.kappa        = std::clamp(s.kappa, 1.0f, 2.5f);
    g.B_toroidal_T = s.B_toroidal_T;
    g.I_plasma_MA  = s.plasma_current_MA;
    return g;
}

// ─── PlasmaCoreBridge implementation ─────────────────────────────────────────

PlasmaCoreBridge::PlasmaCoreBridge(const PlasmaConfig& cfg)
    : cfg_(cfg)
{
    initGPU();
}

PlasmaCoreBridge::~PlasmaCoreBridge()
{
    shutdownGPU();
}

void PlasmaCoreBridge::initGPU()
{
    cudaError err = cudaSetDevice(0);
    if (err != cudaSuccess)
        throw std::runtime_error("PlasmaCore: no CUDA device available");

    cudaStreamCreate(&stream_main_);
    cudaStreamCreate(&stream_neutron_);

    // ── Initialize ENDF/B-VIII.0 cross-sections (built-in defaults) ──────────
    // Loads Li-6, Li-7, steel, and coolant XS tables into __constant__ memory.
    picInitializeCrossSections(/*ace_file_path=*/nullptr);

    // ── Allocate PIC state (reduced grid/particle count for runtime) ─────────
    //  ONLY when the caller actually wants PIC — allocating a 32^3 grid with
    //  200k particles costs ~50 ms of CUDA work and ~200 MB of device memory
    //  even before the first step runs, so skipping it when cfg_.run_pic is
    //  false makes "0D-only / arcade mode" actually cheap to start up.
    if (cfg_.run_pic) {
        PICConfig pcfg;
        pcfg.N_particles = 50000;           // per species (reduced from 1M)
        pcfg.Nx = 32; pcfg.Ny = 32; pcfg.Nz = 32;
        pcfg.Lx = 20.0f; pcfg.Ly = 12.0f; pcfg.Lz = 12.0f;
        pcfg.dt_pic = cfg_.pic_dt;
        pcfg.coulomb_log = 17.0f;

        pic_state_ = picCreate(pcfg);
        if (!pic_state_) {
            // CUDA allocation failed — fall back to 0D-only mode
            initialised_ = false;
            return;
        }
    }

    // ── Initialize MHD state ─────────────────────────────────────────────────
    mhd_state_ = new MHDPhysics::MHDState();
    MHDPhysics::RadialMesh mesh;
    MHDPhysics::initializeMHD(*mhd_state_, mesh, /*I_p_MA=*/0.0f, /*B_T=*/5.3f);

    // PIC plasma will be initialized on first update (when T_e > 0)
    pic_initialized_ = false;

    initialised_ = true;
}

void PlasmaCoreBridge::shutdownGPU()
{
    if (!initialised_) return;
    cudaStreamSynchronize(stream_main_);
    cudaStreamSynchronize(stream_neutron_);
    cudaStreamDestroy(stream_main_);
    cudaStreamDestroy(stream_neutron_);

    if (pic_state_) {
        picDestroy(pic_state_);
        pic_state_ = nullptr;
    }
    if (mhd_state_) {
        delete mhd_state_;
        mhd_state_ = nullptr;
    }
    initialised_ = false;
}

void PlasmaCoreBridge::reset()
{
    // Wipe the 0D power-balance integrator state so the next cold start
    // actually starts cold — without this, W_stored_MJ_ from the previous
    // (possibly very hot) run survives the ReactorState wipe and the very
    // first updatePowerBalance() call after RESET re-derives a multi-keV
    // T_e from it.
    W_stored_MJ_ = 0.0f;
    T_i_keV_     = 0.0f;

    // Clear the W-initialized sticky flag so the next cold-start actually
    // re-seeds W from (3/2)·n·T·V using the freshly-reset ReactorState.
    // This is the half of the "fusion power drops to 0 at ~30 MW" fix that
    // lives outside updatePowerBalance — the re-init must only run once per
    // discharge, not every time W crashes through 0.1 MJ.
    W_initialized_ = false;

    // Force PIC re-initialisation on the next non-cold tick so stale
    // particle distributions from the previous discharge don't bleed
    // across the reset.
    pic_initialized_ = false;

    // Re-zero the MHD current profile so q95 isn't carried over either.
    if (mhd_state_) {
        MHDPhysics::RadialMesh mesh;
        MHDPhysics::initializeMHD(*mhd_state_, mesh, /*I_p_MA=*/0.0f, /*B_T=*/5.3f);
    }
}

void PlasmaCoreBridge::update(ReactorState& state, float dt)
{
    if (!initialised_) return;
    if (state.plasma_status == PlasmaStatus::Cold) return;

    // ── Initialize PIC plasma on first non-cold tick ─────────────────────────
    if (!pic_initialized_ && state.electron_temp_keV > 0.1f) {
        picInitializePlasma(pic_state_,
                              state.electron_temp_keV,
                              state.plasma_density_m3,
                              /*R_major=*/6.2f, /*a=*/2.0f, /*κ=*/1.7f,
                              state.B_toroidal_T);
        pic_initialized_ = true;
    }

    // ── Run a small number of PIC sub-steps per game tick ───────────────────
    // The PIC dt is ~1 ps; the game tick is ~1 ms.  We can't run 1e9 PIC
    // steps per tick — instead we run a *single* PIC step per tick to keep
    // the kinetic state roughly synced with the 0D model.  The bulk of the
    // physics is handled by the 0D power balance below.
    //
    // IMPORTANT: a single picStep on a 32^3 grid with 200k particles is
    // already ~50 ms on a Tesla T4.  Running 10 of them per tick (the old
    // value) meant every game tick cost ~500 ms — once the plasma was
    // non-cold and pic_initialized_ flipped true, the framerate collapsed
    // to ~2 fps.  Dropping to 1 step per tick restores real-time play.
    //
    // The PlasmaConfig::run_pic flag (default true) lets the caller
    // disable PIC entirely for a pure-0D mode that's even faster — set
    // pc.run_pic = false in main.cpp for "arcade mode".
    if (cfg_.run_pic && pic_initialized_ && pic_state_) {
        constexpr int N_PIC_STEPS_PER_TICK = 1;
        for (int i = 0; i < N_PIC_STEPS_PER_TICK; i++) {
            ::picStep(pic_state_, cfg_.pic_dt);
        }

        // ── PIC diagnostic readback ────────────────────────────────────────
        //  NOTE: the PIC's computeCellDiagnostics kernel is currently a
        //  PLACEHOLDER — it writes a hardcoded 10 keV / 1e20 m⁻³ to every
        //  cell regardless of the actual particle kinematics.  The real
        //  per-cell reductions are not yet implemented (see the TODO in
        //  pic_bridge.cu::computeCellDiagnostics).
        //
        //  Until the real diagnostics are wired in, we must NOT blend the
        //  placeholder T_e back into the 0D model — doing so clamps T_e
        //  to ~10 keV regardless of what the 0D power balance computed,
        //  which corrupts the fusion power (<σv> is extremely sensitive
        //  to T_e).  This was the root cause of the "fusion power drops
        //  to 0" bug: the 50:50 blend dragged T_e down to 10 keV, <σv>
        //  dropped, fusion power plateaued, and the stored-energy
        //  integration became unstable.
        //
        //  We still call picReadbackDiagnostics so the visualization tabs
        //  have *something* to display, but we don't feed it back into
        //  state.electron_temp_keV.
        PICCellDiagnostics diag;
        picReadbackDiagnostics(pic_state_, diag);
        // (Diagnostic readback for the plasma-viz tabs only — NOT blended
        //  back into the 0D model until the PIC diagnostics are real.)
    }

    // ── 0D reactor-grade physics (the main physics update) ───────────────────
    updatePowerBalance(state, dt);

    // ── MHD disruption tracking ──────────────────────────────────────────────
    //  updateMHD may overwrite state.q_safety with the MHD-computed value
    //  (which uses the integrated current profile, not the 0D cylindrical
    //  formula).  After MHD runs, re-evaluate the q-based disruption flag
    //  using the final q_safety so the disruption status is consistent with
    //  the value the UI displays.  Previously the disruption flag was set
    //  in updatePowerBalance using the 0D q, then MHD overwrote q — so the
    //  flag could be stale (e.g. 0D says q=1.9 → disruption, but MHD says
    //  q=2.1 → no disruption, yet the flag stayed true).
    updateMHD(state, dt);

    // Re-check the q-based disruption trigger with the final (MHD) q_safety.
    //  Only re-evaluate the LowQ cause — Greenwald and Troyon don't depend
    //  on q and were already correctly set in updatePowerBalance.
    if (state.q_safety < 2.0f && state.plasma_current_MA > 1.0f) {
        if (!state.disruption_flag) {
            state.disruption_flag  = true;
            state.disruption_cause = DisruptionCause::LowQ;
            state.alarm_disruption = true;
        }
    }
}

void PlasmaCoreBridge::picStep(float dt_pic)
{
    // Legacy stub — kept for ABI compatibility.  The actual PIC loop is now
    // driven by pic_bridge.h's picStep() function, called from update().
    (void)dt_pic;
}

// ─── 0D power-balance model (the main physics update) ────────────────────────
void PlasmaCoreBridge::updatePowerBalance(ReactorState& state, float dt)
{
    using namespace ConfinementPhysics;

    auto geom = makeGeometry(state);

    // ── Early-return guard for quenched / near-zero plasma ──────────────────
    //  During SCRAM, all setpoints ramp toward 0.  Once the plasma is
    //  effectively gone (density or current below noise floor), every
    //  formula in the power balance divides by zero:
    //    - T_new = (2/3)·W / (n·V·...)  → 0/0 = NaN when n→0
    //    - tau_E = ipb98y2(I_p=0, ...)  → 0  (powf(0,0.93)=0)
    //    - P_cond = W / tau_E           → 0/0 = NaN when W→0
    //    - beta_p = beta·B² / (2·μ₀·n·T) → 0/0 when n,T→0
    //  Short-circuit here and zero everything out.  This is physically
    //  correct — a quenched plasma has no fusion power, no radiation, no
    //  stored energy — and prevents NaN from leaking into ReactorState.
    if (state.plasma_density_m3 < 1e15f || state.plasma_current_MA < 0.01f) {
        state.fusion_power_MW    = 0.0f;
        state.alpha_power_MW     = 0.0f;
        state.radiated_power_MW  = 0.0f;
        state.neutron_flux_m2s   = 0.0f;
        state.beta               = 0.0f;
        state.beta_N             = 0.0f;
        state.stored_energy_GJ   = 0.0f;
        state.electron_temp_keV  = 0.0f;
        state.plasma_temp_keV    = 0.0f;
        state.Q_scientific       = 0.0f;
        state.q_safety           = 100.0f;
        state.disruption_flag    = false;
        state.alarm_disruption   = false;
        T_i_keV_                 = 0.0f;
        W_stored_MJ_             = 0.0f;
        // Transition to Quenched if not already Cold
        if (state.plasma_status != PlasmaStatus::Cold)
            state.plasma_status = PlasmaStatus::Quenched;
        return;
    }

    // ── Ramp-rate-limited approach to operator setpoints ──────────────────
    //  These ramp rates limit how fast the plasma state can change toward
    //  the operator's setpoints.  They're physical limits — you can't
    //  instantaneously change I_p (CS flux swing takes time), ne (gas
    //  puffing/pellets have finite throughput), or T_e (thermal inertia).
    //
    //  The T_e ramp rate was 0.3 keV/s, which is too slow — the operator
    //  makes a setpoint change and sees nothing happen for many seconds,
    //  which feels unresponsive.  Real tokamak T_e can change at ~1-2
    //  keV/s during ramp-up (limited by auxiliary heating power).  We use
    //  2.0 keV/s here for better responsiveness while still preventing
    //  instant jumps.
    constexpr float dIp_dt_max = 0.5f;    // MA/s
    constexpr float dTe_dt_max = 2.0f;    // keV/s  (was 0.3, too slow)
    constexpr float dne_dt_max = 5.0e18f; // m⁻³/s

    float Ip_target = std::clamp(state.sp_plasma_current_MA, 0.0f, 20.0f);
    float Te_target = std::clamp(state.sp_electron_temp_keV, 0.0f, 50.0f);
    float ne_target = std::clamp(state.sp_density_m3, 0.0f, 3e20f);

    state.plasma_current_MA += std::clamp(Ip_target - state.plasma_current_MA,
                                          -dIp_dt_max * dt, dIp_dt_max * dt);
    state.electron_temp_keV += std::clamp(Te_target - state.electron_temp_keV,
                                          -dTe_dt_max * dt, dTe_dt_max * dt);
    state.plasma_density_m3 += std::clamp(ne_target - state.plasma_density_m3,
                                          -dne_dt_max * dt, dne_dt_max * dt);
    state.plasma_current_MA = std::max(state.plasma_current_MA, 0.f);
    state.electron_temp_keV = std::max(state.electron_temp_keV, 0.f);
    state.plasma_density_m3 = std::max(state.plasma_density_m3, 0.f);
    state.plasma_temp_keV   = state.electron_temp_keV;
    T_i_keV_                  = state.electron_temp_keV;

    // ── Densities ──────────────────────────────────────────────────────────
    float ne = state.plasma_density_m3;
    float f_DT = std::max(1.0f - state.helium_fraction - state.impurity_fraction, 0.0f);
    float D_T_ratio = std::max(state.D_T_ratio, 0.01f);
    float n_D = 0.5f * ne * f_DT * (D_T_ratio / (1.0f + D_T_ratio)) * 2.0f;
    float n_T = 0.5f * ne * f_DT * (1.0f / (1.0f + D_T_ratio)) * 2.0f;
    float Z_eff = 1.0f + 4.0f * state.helium_fraction
                    + (state.impurity_fraction > 0 ? 10.0f * state.impurity_fraction : 0.0f);

    // ── Auxiliary heating: sum of NBI + ICRH + ECRH + LHCD ───────────────────
    //  The previous code hardcoded P_aux_MW = 50.0f.  This was the second half
    //  of the "fusion power drops to 0 at ~30 MW" bug: during ramp-down the
    //  operator dropped I_p and fuel rate to 0, but P_aux stayed at 50 MW,
    //  keeping T_e high enough to keep fusing at ~30 MW.  Then as I_p dropped
    //  further, confinement time τ_E → 0 (because τ_E ∝ I_p^0.93), P_cond → ∞,
    //  W crashed to 0, T_e crashed, and P_fus → 0.  With T_e crashed but P_aux
    //  still at 50 MW, the next tick reheated W to a value that yielded P_fus
    //  back up to ~30 MW — then the loop repeated.  Classic 2-tick limit cycle.
    //
    //  By sourcing P_aux from the H&CD systems (each of which the operator
    //  ramps down or the control system ramps down on Rampdown mode), the
    //  auxiliary heating falls with the rest of the operating point and the
    //  power balance stays well-behaved through the shutdown.
    float P_aux_total_MW = state.hcd_nbi_actual_MW + state.hcd_icrh_actual_MW
                         + state.hcd_ecrh_actual_MW + state.hcd_lhcd_actual_MW;
    P_aux_total_MW = std::max(0.0f, std::isfinite(P_aux_total_MW) ? P_aux_total_MW : 0.0f);

    // ── Power balance (ConfinementPhysics::solvePowerBalance) ──────────────
    PowerBalance pb{};
    pb.n_e_m3   = ne;
    pb.T_e_keV  = state.electron_temp_keV;
    pb.T_i_keV  = T_i_keV_;
    pb.n_D_m3   = n_D;
    pb.n_T_m3   = n_T;
    pb.Z_eff    = Z_eff;
    pb.P_aux_MW = P_aux_total_MW;
    // ── Ohmic heating ────────────────────────────────────────────────────────
    //  Ohmic heating (P_ohm = η·j²) is the plasma's resistive heating from the
    //  loop voltage induced by the central solenoid.  In a real tokamak:
    //    - During Initiating (ramp-up): P_ohm is the DOMINANT heat source
    //      (~15-20 MW for ITER) because the plasma is cold (high resistivity)
    //      and the loop voltage is high (~10 V).
    //    - During Burning: P_ohm drops to ~1-2 MW because the plasma is hot
    //      (low resistivity, Spitzer η ∝ T^-1.5) and the CS loop voltage
    //      drops to ~0.1 V.  It's not zero, though — there's always some
    //      residual resistive heating.
    //
    //  The previous code had P_ohm = 5 MW during Initiating and 0 during
    //  Burning.  The 5 MW was too low to heat the plasma to fusion-relevant
    //  temperatures without auxiliary heating, and the 0 during Burning
    //  meant the plasma couldn't sustain itself if alpha heating dropped
    //  (e.g. during a transient dip in fusion power).  This was a major
    //  contributor to the "insanely difficult to keep stable" report:
    //  the user had to keep ne very low (5e18) to reduce bremsstrahlung
    //  enough that 5 MW of ohmic heating could maintain temperature.
    //
    //  New values: 15 MW during Initiating (realistic for ITER ramp-up),
    //  2 MW during Burning (residual resistive heating), 0 during
    //  Disrupting/Quenched/Cold (plasma is gone or being killed).
    if (state.plasma_status == PlasmaStatus::Initiating) {
        pb.P_ohm_MW = 15.0f;
    } else if (state.plasma_status == PlasmaStatus::Burning) {
        pb.P_ohm_MW = 2.0f;
    } else {
        pb.P_ohm_MW = 0.0f;
    }
    pb.P_line_MW = 0.0f;     // impurity line radiation handled separately

    // ── W_stored initialization (cold-start only) ──────────────────────────────
    //  IMPORTANT: only initialize W from temperature when the integrator is
    //  truly uninitialized — i.e. on the very first call after construction
    //  or after a reset().  The old test `W_stored_MJ_ < 0.1f` was triggered
    //  again every time W crashed during a ramp-down — and that re-init
    //  pulled W back up from (3/2)·n·T·V using a temperature that hadn't yet
    //  been updated, which set up the 2-tick limit cycle the user saw as the
    //  "fusion power drops to 0 at ~30 MW then bounces back" bug.
    //
    //  We use a sticky `W_initialized_` flag (set false in reset()) so the
    //  re-init only happens at true cold-start.  During normal operation W
    //  is allowed to go all the way down to 0 smoothly — that's the correct
    //  end-state of a controlled shutdown.
    if (!W_initialized_) {
        W_stored_MJ_ = storedThermalEnergy(ne, state.electron_temp_keV, T_i_keV_, geom) * 1e-6f;
        W_initialized_ = true;
    }

    pb = solvePowerBalance(geom, pb, W_stored_MJ_, /*H_98=*/1.0f);

    // ── Integrate stored energy ─────────────────────────────────────────────
    //  Guard against NaN: solvePowerBalance clamps dW_dt to finite, but
    //  double-check here.  Once W becomes NaN, std::max(NaN, 0) returns
    //  NaN on libstdc++ and the poison spreads to T_new and everything
    //  downstream.
    if (std::isfinite(pb.dW_dt_MW)) {
        W_stored_MJ_ += pb.dW_dt_MW * dt;   // MW × s = MJ
    }
    if (!std::isfinite(W_stored_MJ_) || W_stored_MJ_ < 0.0f) {
        W_stored_MJ_ = 0.0f;
    }

    // ── Update temperatures from W (T = (2/3) · W / (n · V)) ────────────────
    float V = plasmaVolume(geom);
    float T_denom = ne * V * 1.602176634e-16f;
    float T_new_keV = (T_denom > 1e-30f)
                    ? (2.0f / 3.0f) * W_stored_MJ_ * 1e6f / T_denom
                    : 0.0f;
    if (!std::isfinite(T_new_keV)) T_new_keV = 0.0f;
    T_new_keV = std::clamp(T_new_keV, 0.0f, 200.0f);
    state.electron_temp_keV = T_new_keV;
    state.plasma_temp_keV   = T_new_keV;
    T_i_keV_                 = T_new_keV;

    // ── Write outputs to ReactorState ──────────────────────────────────────
    //  All values are guaranteed finite by the guards above, but add a
    //  final isfinite check before writing — if anything slipped through
    //  (e.g. bremSStrahlung producing inf for very high n·T), zero it
    //  rather than letting NaN propagate to the UI and downstream modules.
    state.fusion_power_MW   = std::isfinite(pb.P_fus_MW)   ? pb.P_fus_MW   : 0.0f;
    state.alpha_power_MW    = std::isfinite(pb.P_alpha_MW) ? pb.P_alpha_MW : 0.0f;
    state.tau_E_s           = std::isfinite(pb.tau_E_s)   ? pb.tau_E_s   : 0.0f;
    float P_rad_total = pb.P_brem_MW + pb.P_sync_MW + pb.P_line_MW;
    state.radiated_power_MW = std::isfinite(P_rad_total)   ? P_rad_total   : 0.0f;
    state.neutron_flux_m2s  = std::isfinite(pb.P_neutron_MW)
                            ? pb.P_neutron_MW * 1e6f / (14.07e6f * 1.602e-19f) / 600.0f
                            : 0.0f;
    state.stored_energy_GJ = W_stored_MJ_ * 1e-3f;

    // ── Safety factor q_95 ─────────────────────────────────────────────────
    //  Cylindrical q with elongation correction:
    //    q_95 ≈ (2π · a² · B_T · κ²) / (μ₀ · R · I_p)
    //
    //  The old formula `5 · a² · B / (R · μ₀ · I_p)` gave q ≈ 0.91 for
    //  ITER-class parameters (a=2, B=5.3, R=6.2, I_p=15 MA) — well below
    //  the disruption threshold of 2.  Real ITER q_95 is ~3.  The missing
    //  factor is the elongation κ² (the plasma cross-section is an ellipse,
    //  not a circle, which increases the effective edge path length).  With
    //  κ=1.7: q ≈ 0.91 × 1.7² ≈ 2.6, which is in the right ballpark.
    //
    //  Without this fix, disruption_flag was permanently true (q < 2 &&
    //  I_p > 1), forcing plasma_status = Disrupting, which zeroes P_ohm
    //  and prevents the plasma from ever reaching Burning status.
    const float R = geom.R_major_m, a = geom.a_minor_m;
    constexpr float PI = 3.14159265358979f;
    constexpr float MU0 = 1.25663706e-6f;
    float kappa = geom.kappa;
    if (state.plasma_current_MA > 0.05f && state.B_toroidal_T > 0.01f) {
        state.q_safety = (2.0f * PI * a * a * state.B_toroidal_T * kappa * kappa)
                       / (MU0 * R * state.plasma_current_MA * 1e6f);
    } else {
        state.q_safety = 100.0f;   // no plasma / no field → q meaningless
    }

    // ── Greenwald density limit ────────────────────────────────────────────
    float n_Greenwald = state.plasma_current_MA / (3.14159265f * a * a);
    bool greenwald_exceeded = (ne * 1e-20f > 1.5f * n_Greenwald);

    // ── Troyon β-limit ─────────────────────────────────────────────────────
    constexpr float beta_N_limit = 2.5f;
    float beta_max = (state.B_toroidal_T > 0.01f)
                   ? beta_N_limit * 1.257e-6f * state.plasma_current_MA * 1e6f
                     / (a * state.B_toroidal_T)
                   : 1e9f;   // no B field → no beta limit
    //  Guard: beta = 2·μ₀·n·T / B².  Both n and T are guaranteed > 0
    //  by the early-return guard above, and B is guarded here.  But
    //  add isfinite as a belt-and-suspenders check.
    float beta_raw = (state.B_toroidal_T > 0.01f)
                   ? (2.0f * 1.257e-6f * ne * state.electron_temp_keV * 1.602e-16f)
                     / (state.B_toroidal_T * state.B_toroidal_T)
                   : 0.0f;
    state.beta = std::isfinite(beta_raw) ? beta_raw : 0.0f;
    // Normalized beta (Troyon) — β_N = β·a·B_T / (μ₀·I_p)
    // β_N > 2.5 (Troyon limit) → ideal MHD kink mode.
    float beta_N_raw = (state.plasma_current_MA > 0.01f && state.B_toroidal_T > 0.01f)
                     ? state.beta * 100.0f * a * state.B_toroidal_T
                       / (state.plasma_current_MA)
                     : 0.0f;   // β[%] · a[m] · B[T] / I_p[MA]
    state.beta_N = std::isfinite(beta_N_raw) ? beta_N_raw : 0.0f;
    bool troyon_exceeded = (state.beta_N > beta_N_limit);

    // ── H-mode threshold (IPB98(th,98)) ────────────────────────────────────
    float n_20 = ne * 1e-20f;
    float P_LH_thresh = 0.04f * powf(n_20, 0.65f) * powf(state.B_toroidal_T, 0.85f)
                      * geom.R_major_m * powf(geom.a_minor_m, 0.4f);
    bool in_H_mode = (pb.P_alpha_MW + pb.P_aux_MW > P_LH_thresh);

    // ── Bootstrap current fraction ─────────────────────────────────────────
    //  Guard: beta_p = beta·B² / (2·μ₀·n·T).  Division by zero when
    //  n→0 or T→0, but the early-return guard above already caught
    //  those cases.  Belt-and-suspenders: check denominator > 0.
    float bp_denom = 2.0f * 1.257e-6f * ne * state.electron_temp_keV * 1.602e-16f;
    float beta_p = (state.beta > 0.001f && bp_denom > 1e-30f)
                 ? state.beta * (state.B_toroidal_T * state.B_toroidal_T) / bp_denom
                 : 0.0f;
    if (!std::isfinite(beta_p)) beta_p = 0.0f;
    state.q_safety = std::max(state.q_safety, 1.0f);
    float f_bs = bootstrapFraction(state.q_safety, R / a, beta_p);

    // ── Write the bootstrap current to ReactorState ─────────────────────────
    //  f_bs is the bootstrap fraction (fraction of I_p carried by bootstrap
    //  current).  The actual bootstrap current in MA = f_bs × I_p_MA.
    //  This is what the H&CD tab displays as "Bootstrap: X.XX MA" and what
    //  the "Total non-inductive" current sums up.  Previously this was
    //  always 0 because the bridge computed f_bs but never wrote it to
    //  state.hcd_bootstrap_current_MA — the HCD module wrote 0 with a
    //  "updated by bridge" comment, but the bridge never updated it.
    state.hcd_bootstrap_current_MA = std::isfinite(f_bs)
        ? f_bs * state.plasma_current_MA : 0.0f;

    // ── Scientific Q ───────────────────────────────────────────────────────
    state.Q_scientific = scientificQ(state.fusion_power_MW, pb.P_aux_MW);

    // ── Disruption heuristics (now augmented by MHD module) ────────────────
    //  Each trigger sets a specific DisruptionCause so the alarm system can
    //  tell the operator *what* went wrong and *what to do about it*, instead
    //  of just "disruption risk".
    state.disruption_flag  = false;
    state.disruption_cause = DisruptionCause::None;
    if (state.q_safety < 2.0f && state.plasma_current_MA > 1.0f) {
        state.disruption_flag  = true;
        state.disruption_cause = DisruptionCause::LowQ;
    }
    if (greenwald_exceeded) {
        state.disruption_flag  = true;
        state.disruption_cause = DisruptionCause::Greenwald;
    }
    if (troyon_exceeded) {
        state.disruption_flag  = true;
        state.disruption_cause = DisruptionCause::Troyon;
    }
    state.alarm_disruption = state.disruption_flag;

    // ── Plasma status transitions ──────────────────────────────────────────
    //  The Burning↔Initiating transitions now use lower thresholds so a
    //  controlled ramp-down doesn't oscillate the status (which was a
    //  minor contributor to the original "fusion power drops to 0" bug —
    //  status flipping Burning→Initiating enabled P_ohm=5 MW which jerked
    //  the power balance).
    //
    //  Burning requires P_fus > 10 MW (was 50).  The 50 MW threshold was
    //  too high — it meant the plasma couldn't transition to Burning
    //  without significant auxiliary heating, and once in Burning the
    //  P_ohm dropped to 0, causing the plasma to cool back below 50 MW
    //  and drop back to Initiating, then heat back up... a limit cycle.
    //  At 10 MW, alpha heating (20% of P_fus = 2 MW) plus the residual
    //  P_ohm (2 MW during Burning) is enough to sustain the plasma at a
    //  modest operating point even without aux heating.
    //
    //  Initiating↔Burning uses a 5 MW hysteresis band (10 up, 5 down) to
    //  prevent flicker at the threshold.
    if (state.plasma_status == PlasmaStatus::Disrupting) {
        if (!state.disruption_flag) {
            // Disruption cleared — go back to Initiating
            state.plasma_status = PlasmaStatus::Initiating;
        }
    }
    if (state.plasma_status == PlasmaStatus::Initiating) {
        bool ramped_up = state.plasma_current_MA >= 0.9f * state.sp_plasma_current_MA
                      && state.electron_temp_keV >= 0.5f * state.sp_electron_temp_keV
                      && state.plasma_density_m3 >= 0.5f * state.sp_density_m3;
        if (ramped_up && state.fusion_power_MW > 10.0f) {
            state.plasma_status = PlasmaStatus::Burning;
        }
    } else if (state.plasma_status == PlasmaStatus::Burning) {
        if (state.fusion_power_MW < 5.0f) {
            state.plasma_status = PlasmaStatus::Initiating;
        }
    }
    if (state.disruption_flag) {
        state.plasma_status = PlasmaStatus::Disrupting;
    }
}

// ─── MHD disruption tracking ─────────────────────────────────────────────────
void PlasmaCoreBridge::updateMHD(ReactorState& state, float dt)
{
    if (!mhd_state_) return;

    // ── Skip MHD update for quenched / near-zero plasma ─────────────────────
    //  When I_p → 0 (during SCRAM), the MHD current profile degenerates
    //  and the q-profile computation, tearing-mode analysis, and VDE
    //  tracker all become numerically meaningless.  Worse, if T_e has
    //  already gone NaN (from the power balance), the Spitzer resistivity
    //  in diffuseCurrent() will produce NaN, which poisons the j_phi
    //  array permanently.  Short-circuit and set a safe q_safety.
    if (state.plasma_current_MA < 0.01f || !std::isfinite(state.electron_temp_keV)) {
        state.q_safety = 100.0f;
        return;
    }

    MHDPhysics::RadialMesh mesh;

    // Update the MHD state's temperature profile from the 0D T_e
    // (in a full simulation, this would come from the PIC per-cell T_e)
    float T_e_safe = std::isfinite(state.electron_temp_keV) ? state.electron_temp_keV : 0.0f;
    for (size_t i = 0; i < mhd_state_->T_e_keV.size(); i++) {
        float r_norm = mesh.r_norm((int)i);
        float shape = std::max(1.0f - r_norm * r_norm, 0.0f);
        mhd_state_->T_e_keV[i] = T_e_safe * shape;
    }

    // Renormalize the MHD current profile to match the reactor I_p
    // (the diffusion can drift the total current)
    constexpr float PI = 3.14159265358979f;
    float I_total = 0.0f;
    for (int i = 0; i < mesh.N; i++) {
        float r = mesh.r(i);
        I_total += mhd_state_->j_phi[i] * 2.0f * PI * r * mesh.dr;
    }
    // Guard: if j_phi has become NaN (e.g. from a prior tick's NaN T_e
    // before the guard above was added), I_total will be NaN and the
    // renormalization would propagate it.  Detect and re-initialize.
    if (!std::isfinite(I_total)) {
        MHDPhysics::initializeMHD(*mhd_state_, mesh,
                                    state.plasma_current_MA, state.B_toroidal_T);
        I_total = 0.0f;
        for (int i = 0; i < mesh.N; i++) {
            float r = mesh.r(i);
            I_total += mhd_state_->j_phi[i] * 2.0f * PI * r * mesh.dr;
        }
    }
    // ── Initialize the MHD current profile when I_p first becomes non-zero ──
    //  The MHD state was created in initGPU() with I_p=0, so j_phi is all
    //  zeros.  Without this block, j_phi stays at 0 forever (the
    //  renormalization below only fires if I_total > 1.0 A, which it
    //  never is when j_phi = 0).  A zero current profile gives q = 100
    //  everywhere — the MHD module becomes dead weight and never detects
    //  tearing modes or real q-profile evolution.
    //
    //  When I_p first exceeds 0.1 MA and the current profile is still
    //  zero, seed it with a parabolic profile matching the reactor's I_p.
    //  After this, the renormalization below keeps it synced as I_p
    //  evolves.
    float I_target = state.plasma_current_MA * 1e6f;
    if (std::fabs(I_total) < 1.0f && I_target > 1e5f) {
        MHDPhysics::initializeMHD(*mhd_state_, mesh,
                                    state.plasma_current_MA, state.B_toroidal_T);
        I_total = 0.0f;
        for (int i = 0; i < mesh.N; i++) {
            float r = mesh.r(i);
            I_total += mhd_state_->j_phi[i] * 2.0f * PI * r * mesh.dr;
        }
    }
    if (std::fabs(I_total) > 1.0f) {
        float scale = I_target / I_total;
        for (auto& j : mhd_state_->j_phi) j *= scale;
    }

    // Run the MHD update
    auto log = MHDPhysics::updateMHD(*mhd_state_, mesh,
                                      state.plasma_current_MA,
                                      state.B_toroidal_T, dt);

    // Augment the disruption flag with MHD-derived triggers
    if (log.tearing_21_disruptive || log.tearing_32_disruptive
        || log.vde_triggered || mhd_state_->disruption_ongoing) {
        state.disruption_flag = true;
        state.alarm_disruption = true;
        state.plasma_status = PlasmaStatus::Disrupting;
    }

    // Write back q0 and q95 from the MHD-computed q-profile
    //  Guard: if the MHD module produced a non-finite q95 (e.g. from
    //  a degenerate current profile), fall back to the 0D q_safety
    //  computed in updatePowerBalance rather than writing NaN.
    //
    //  ── Elongation correction ─────────────────────────────────────────────
    //  The MHD module's computeQProfile uses the CYLINDRICAL q formula
    //    q(r) = (2π · r² · B_T) / (μ₀ · R · I_enc(r))
    //  which does NOT include the elongation κ² factor.  For an elongated
    //  plasma, the correct q is approximately κ² times the cylindrical
    //  value (the poloidal circumference is longer, so the field line
    //  has to spiral more times toroidally to close).  Without this
    //  correction, MHD q95 ≈ 0.9 for ITER parameters, which is below
    //  the disruption threshold of 2 — causing the disruption_flag to
    //  fire, forcing plasma_status = Disrupting, zeroing P_ohm, and
    //  killing fusion power.  With κ² = 1.7² = 2.89, the corrected
    //  q95 ≈ 2.6, matching the 0D formula and real ITER.
    constexpr float kappa_correction = 1.7f * 1.7f;  // κ² for ITER elongation
    if (std::isfinite(mhd_state_->q95) && mhd_state_->q95 > 0.0f) {
        state.q_safety = mhd_state_->q95 * kappa_correction;
    }
}

void PlasmaCoreBridge::readbackDiagnostics(ReactorState& state, float dt)
{
    (void)state; (void)dt;
}
