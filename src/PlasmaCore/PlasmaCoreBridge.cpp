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

// ─── CS transformer flux budget ──────────────────────────────────────────────
//  Volt-seconds available for RESISTIVE flux consumption during the discharge
//  (the inductive part of the swing — L_p·I_p — is treated as supplied by the
//  PF system and recoverable).  ITER's CS provides roughly this much for the
//  burn; V_loop ≈ 0.075 V at burning temperatures gives the design ~1000 s
//  pulse.  Driving non-inductive current (NBI/LHCD/bootstrap) stretches it.
static constexpr float CS_FLUX_BUDGET_WB = 90.0f;

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
    W_e_MJ_      = 0.0f;
    W_i_MJ_      = 0.0f;
    T_i_keV_     = 0.0f;
    h_mode_      = false;
    sawtooth_clock_s_ = 0.0f;
    sawtooth_count_   = 0;
    elm_clock_s_      = 0.0f;
    elm_count_        = 0;
    flux_used_Wb_     = 0.0f;
    re_seeded_        = false;
    re_conv_frac_     = 0.0f;
    gw_violation_s_   = 0.0f;
    disrupt_flag_s_   = 0.0f;
    quench_committed_ = false;

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
        state.ion_temp_keV       = 0.0f;
        state.h_mode             = false;
        state.H98_factor         = 1.0f;
        state.P_LH_MW            = 0.0f;
        state.P_ohm_MW           = 0.0f;
        state.P_aux_actual_MW    = 0.0f;
        state.P_brem_MW          = 0.0f;
        state.P_sync_MW          = 0.0f;
        state.P_line_MW          = 0.0f;
        state.P_cond_MW          = 0.0f;
        state.dW_dt_MW           = 0.0f;
        state.loop_voltage_V     = 0.0f;
        state.greenwald_frac     = 0.0f;
        state.sawtooth_period_s  = 0.0f;
        state.elm_freq_Hz        = 0.0f;
        //  Clear the mitigation latch — the quench it governed is complete.
        //  (The DM tab's FIRED lamps and last-fire info keep the history.)
        state.dm_mitigated       = false;
        h_mode_                  = false;
        T_i_keV_                 = 0.0f;
        W_stored_MJ_             = 0.0f;
        W_e_MJ_                  = 0.0f;
        W_i_MJ_                  = 0.0f;
        elm_clock_s_             = 0.0f;

        // ── Release the quench commitment — the plasma is gone ───────────────
        //  BUG FIX: the committed-quench release used to live only at the
        //  bottom of this function, gated on I_p < 0.008 MA — but THIS guard
        //  returns early at I_p < 0.01 MA, so the release was unreachable.
        //  quench_committed_ then stayed latched forever: every re-INITIATEd
        //  discharge was forced straight back into Disrupting and died,
        //  soft-locking the game until a full cold RESET.
        quench_committed_ = false;
        disrupt_flag_s_   = 0.0f;

        // ── Runaway-electron wall strike ─────────────────────────────────────
        //  The current quench is complete.  If an RE beam formed during it,
        //  position control is lost the moment the bulk plasma vanishes and
        //  the beam strikes the first wall, dumping its magnetic energy
        //  E = ½·L_p·I_RE² (hundreds of MJ at multi-MA beam currents) into a
        //  localized wetted spot.  ThermalHydraulics consumes the pending
        //  energy and spikes the wall temperature.
        if (state.re_beam_MA > 0.5f) {
            constexpr float L_p_H = 1.25e-5f;   // ITER-class plasma inductance
            float I_re_A = state.re_beam_MA * 1e6f;
            state.re_wall_energy_MJ += 0.5f * L_p_H * I_re_A * I_re_A * 1e-6f;
            state.re_wall_strike = true;
            state.re_strike_count++;
            state.alarm_runaway  = true;
        }
        state.re_beam_MA = 0.0f;
        re_seeded_       = false;
        re_conv_frac_    = 0.0f;

        // ── Recharge the central solenoid between shots ──────────────────────
        //  With no plasma, the CS can be swung back to full pre-charge; the
        //  next discharge gets the full volt-second budget.
        flux_used_Wb_               = 0.0f;
        state.cs_flux_remaining_Wb  = CS_FLUX_BUDGET_WB;
        state.cs_flux_exhausted     = false;

        // Transition to Quenched if not already Cold
        if (state.plasma_status != PlasmaStatus::Cold)
            state.plasma_status = PlasmaStatus::Quenched;
        return;
    }

    // ── Disruption persistence (evaluated on LAST tick's final flag) ────────
    //  state.disruption_flag at this point is the flag as of the END of the
    //  previous tick — including sources set after the power balance ran
    //  (the MHD tearing/VDE tracker, mitigation systems).  Accumulating
    //  here, before the heuristics below reset the flag, means a
    //  disruption signalled by ANY module for >100 ms commits the quench.
    //  (Accumulating after the reset only ever saw this function's own
    //  triggers: an MHD-flagged disruption could force Disrupting status
    //  every tick while the setpoint ramps regrew the current — a zombie
    //  plasma, permanently "disrupting" yet burning.)
    //  Commitment threshold: 0.5 s of sustained flagging.  The MHD
    //  tearing tracker raises transient ~100 ms disruption_ongoing
    //  windows while the current profile settles during ramp-up — those
    //  are model noise and must not kill the discharge.  A genuine locked
    //  mode or density-limit collapse flags continuously and commits.
    //  (MGI/SPI firing commits instantly, independent of this timer.)
    if (state.disruption_flag) disrupt_flag_s_ += dt;
    else                       disrupt_flag_s_  = 0.0f;
    if (disrupt_flag_s_ > 0.5f) quench_committed_ = true;

    // ── Ramp-rate-limited approach to operator setpoints ──────────────────
    //  These ramp rates limit how fast the plasma state can change toward
    //  the operator's setpoints.  They're physical limits — you can't
    //  instantaneously change I_p (CS flux swing takes time) or ne (gas
    //  puffing/pellets have finite throughput).
    //
    //  REALISM UPGRADE: the T_e setpoint ramp has been REMOVED.  It used
    //  to nudge state.electron_temp_keV toward sp_electron_temp_keV every
    //  tick, which — once the two-temperature energy accounting was made
    //  consistent — turned out to be a hidden free-energy source of up to
    //  ~40 MW (1.5·n·ΔT·V per second, from nowhere).  A real tokamak has
    //  no temperature dial: T_e is EARNED through ohmic + H&CD heating
    //  against radiation and transport.  The sp_electron_temp_keV setpoint
    //  is retained as the operator's target indication for the UI and the
    //  control system, but the bridge no longer forces it.
    constexpr float dIp_dt_max = 0.5f;    // MA/s
    constexpr float dne_dt_max = 5.0e18f; // m⁻³/s

    float Ip_target = std::clamp(state.sp_plasma_current_MA, 0.0f, 20.0f);
    float ne_target = std::clamp(state.sp_density_m3, 0.0f, 3e20f);

    if (quench_committed_ && state.plasma_status == PlasmaStatus::Disrupting) {
        // ── Current quench ───────────────────────────────────────────────────
        //  A committed disruption decays I_p exponentially — the control
        //  system has lost the plasma and the setpoints are irrelevant.
        //  Mitigated quenches (MGI/SPI fired) are the slow, controlled kind
        //  (τ_CQ ≈ 100 ms → smaller eddy/halo forces); unmitigated quenches
        //  are violent (τ_CQ ≈ 30 ms).  This is the entire point of the
        //  mitigation systems.  Particles dump to the wall on a slower
        //  timescale.
        float tau_cq = state.dm_mitigated ? 0.10f : 0.03f;

        // ── Runaway-electron avalanche ───────────────────────────────────────
        //  The fast current decay induces a huge toroidal E-field (V_loop of
        //  hundreds of volts), which avalanche-multiplies seed electrons into
        //  a relativistic beam.  The conversion fraction is decided once, at
        //  quench commit: exponential avalanche gain ~ exp(2.5·I_p[MA]) makes
        //  this a high-current problem — negligible below ~4 MA, and up to
        //  ~60% of the pre-disruption current at ITER's 15 MA.  MGI/SPI
        //  mitigation raises the density enough that collisional drag kills
        //  the avalanche — the second reason those systems exist.
        if (!re_seeded_) {
            re_seeded_ = true;
            re_conv_frac_ = state.dm_mitigated
                ? 0.0f
                : std::clamp(0.05f * (state.plasma_current_MA - 4.0f), 0.0f, 0.6f);
        }
        float I_before = state.plasma_current_MA;
        state.plasma_current_MA *= expf(-dt / tau_cq);
        state.plasma_density_m3 *= expf(-dt / (4.0f * tau_cq));
        //  The beam picks up a fraction of every ampere the quench sheds.
        state.re_beam_MA += re_conv_frac_ * (I_before - state.plasma_current_MA);
        if (state.re_beam_MA > 0.5f) state.alarm_runaway = true;
    } else if (state.plasma_status == PlasmaStatus::Quenched
            || state.plasma_status == PlasmaStatus::Cold) {
        //  No plasma — the setpoints don't drive anything until the
        //  operator INITIATEs a new discharge.  (Previously the ramps ran
        //  unconditionally and re-grew the current of a quenched plasma
        //  toward the stale setpoint in the background.)
        state.plasma_current_MA = std::max(state.plasma_current_MA - 2.0f * dt, 0.0f);
        state.plasma_density_m3 = std::max(state.plasma_density_m3 * expf(-dt / 0.5f), 0.0f);
    } else {
        // ── CS flux exhaustion: no inductive drive left ─────────────────────
        //  Once the transformer swing is spent, the only current the machine
        //  can hold is what the non-inductive systems carry (NBI + LHCD
        //  current drive + bootstrap).  I_p relaxes toward that level — the
        //  classic end-of-pulse decay of a conventional tokamak.  Operators
        //  who invested in LHCD/NBI current drive and a high-bootstrap
        //  operating point keep burning; purely ohmic discharges wind down.
        if (state.cs_flux_exhausted) {
            float I_NI = state.hcd_nbi_current_drive_MA
                       + state.hcd_lhcd_current_drive_MA
                       + state.hcd_bootstrap_current_MA;
            Ip_target = std::min(Ip_target, std::max(I_NI, 0.0f));
        }
        state.plasma_current_MA += std::clamp(Ip_target - state.plasma_current_MA,
                                              -dIp_dt_max * dt, dIp_dt_max * dt);
        state.plasma_density_m3 += std::clamp(ne_target - state.plasma_density_m3,
                                              -dne_dt_max * dt, dne_dt_max * dt);
    }
    state.plasma_current_MA = std::max(state.plasma_current_MA, 0.f);
    state.electron_temp_keV = std::max(state.electron_temp_keV, 0.f);
    state.plasma_density_m3 = std::max(state.plasma_density_m3, 0.f);
    //  T_i is not slaved to T_e: the ion channel evolves purely from
    //  alpha/NBI/ICRH heating plus collisional equilibration.  On true
    //  cold start both channels are seeded together.
    if (!W_initialized_) T_i_keV_ = state.electron_temp_keV;

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
    // ── Auxiliary heating: NBI + ICRH + ECRH + LHCD ──────────────────────────
    //  Sourced from the H&CD systems (never a hardcoded 50 MW — see the
    //  Round-3 notes in git history for the limit-cycle bug that caused).
    //  The per-system powers are split by heated SPECIES when the
    //  PowerBalance struct is filled below, which is why there is no
    //  lumped total computed here any more.

    // ── Power balance (ConfinementPhysics::solvePowerBalance) ──────────────
    PowerBalance pb{};
    pb.n_e_m3   = ne;
    pb.T_e_keV  = state.electron_temp_keV;
    pb.T_i_keV  = T_i_keV_;
    pb.n_D_m3   = n_D;
    pb.n_T_m3   = n_T;
    pb.Z_eff    = Z_eff;
    // ── Species-resolved auxiliary heating ─────────────────────────────────
    //  Each H&CD system heats a different species, which is what makes the
    //  actuator choice matter physically:
    //    ECRH → electrons (electron cyclotron resonance)
    //    LHCD → electrons (Landau damping on electrons)
    //    ICRH → ions      (minority/second-harmonic ion heating)
    //    NBI  → ~50/50    (1 MeV beams slow on both; Stix split ≈ half each
    //                       at reactor T_e — good enough for a 0D model)
    pb.P_aux_e_MW = state.hcd_ecrh_actual_MW + state.hcd_lhcd_actual_MW
                  + 0.5f * state.hcd_nbi_actual_MW;
    pb.P_aux_i_MW = state.hcd_icrh_actual_MW + 0.5f * state.hcd_nbi_actual_MW;
    if (!std::isfinite(pb.P_aux_e_MW)) pb.P_aux_e_MW = 0.0f;
    if (!std::isfinite(pb.P_aux_i_MW)) pb.P_aux_i_MW = 0.0f;
    pb.P_aux_e_MW = std::max(pb.P_aux_e_MW, 0.0f);
    pb.P_aux_i_MW = std::max(pb.P_aux_i_MW, 0.0f);
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
    //  REALISM UPGRADE: the status-based step (15 MW Initiating / 2 MW
    //  Burning / 0 otherwise) is replaced by the actual resistive heating
    //  P_ohm = I_p² · R_plasma with neoclassical Spitzer resistivity.
    //  Because η ∝ T_e^−1.5, this naturally delivers tens of MW into a
    //  cold ~100 eV start-up plasma and decays smoothly to ≲1 MW at
    //  burning temperatures — no discontinuity when the status flips, so
    //  the old Burning↔Initiating limit-cycle mechanism is gone entirely.
    //  The loop voltage V = I_p·R_p is published for the UI.
    if ((state.plasma_status == PlasmaStatus::Initiating
      || state.plasma_status == PlasmaStatus::Burning)
     && !state.cs_flux_exhausted) {
        auto ohm = ohmicHeating(geom, state.plasma_current_MA,
                                state.electron_temp_keV, Z_eff);
        pb.P_ohm_MW          = ohm.P_ohm_MW;
        state.loop_voltage_V = ohm.V_loop_V;

        // ── CS volt-second consumption ──────────────────────────────────────
        //  The transformer only has to supply the flux for the INDUCTIVELY
        //  driven part of the current: dΨ/dt = V_loop · (1 − I_NI/I_p).
        //  Non-inductive current (NBI/LHCD CD + bootstrap, from the previous
        //  tick) directly reduces the burn's flux consumption — this is what
        //  makes the H&CD current-drive systems worth their electricity.
        float I_NI = state.hcd_nbi_current_drive_MA
                   + state.hcd_lhcd_current_drive_MA
                   + state.hcd_bootstrap_current_MA;
        float f_ind = 1.0f;
        if (state.plasma_current_MA > 0.1f)
            f_ind = std::clamp(1.0f - I_NI / state.plasma_current_MA, 0.0f, 1.0f);
        if (std::isfinite(ohm.V_loop_V) && ohm.V_loop_V > 0.0f)
            flux_used_Wb_ += ohm.V_loop_V * f_ind * dt;
        if (flux_used_Wb_ >= CS_FLUX_BUDGET_WB) {
            flux_used_Wb_ = CS_FLUX_BUDGET_WB;
            state.cs_flux_exhausted = true;
        }
    } else {
        pb.P_ohm_MW          = 0.0f;
        state.loop_voltage_V = 0.0f;
    }
    state.cs_flux_remaining_Wb = CS_FLUX_BUDGET_WB - flux_used_Wb_;

    // ── Impurity line radiation ────────────────────────────────────────────
    //  P_line = n_e²·(f_lowZ·L_lowZ(T) + f_W·L_W(T))·V.  The low-Z curve
    //  peaks below ~100 eV, which produces the classic start-up RADIATION
    //  BARRIER: a cold, dirty plasma radiates harder than it can be heated
    //  and collapses.  Boronization (which reduces impurity_fraction)
    //  now has a direct, physical payoff.  A small tungsten trace (0.5%
    //  of the impurity inventory, representing divertor sputtering) keeps
    //  a realistic residual core radiation at burning temperatures.
    {
        float f_lowZ = std::max(state.impurity_fraction, 0.0f);
        //  Tungsten trace: 0.1% of the low-Z impurity inventory, representing
        //  divertor sputtering.  At 2% total impurity this is a 2e-5 W
        //  fraction → ~15-20 MW of core radiation at burn — consistent with
        //  ITER W-accumulation projections, and a solid reason to keep the
        //  walls conditioned.
        float f_W    = 0.001f * f_lowZ;
        //  MGI/SPI radiator: injected Ne/Ar never fully strips at these
        //  temperatures, so it radiates like the flat high-Z curve.  A
        //  0.25-0.30 fractional dose drives hundreds of MW of line
        //  radiation — the radiative thermal quench that IS the mitigation.
        //  The inventory pumps out on a ~10 s timescale.
        float f_dm = std::max(state.dm_injected_impurity, 0.0f);
        if (f_dm > 0.0f) {
            state.dm_injected_impurity = f_dm * expf(-dt / 10.0f);
            if (state.dm_injected_impurity < 1e-4f) state.dm_injected_impurity = 0.0f;
        }
        float p_line_W = lineRadiationDensity(ne, state.electron_temp_keV,
                                              f_lowZ, f_W + f_dm)
                       * plasmaVolume(geom);
        pb.P_line_MW = std::isfinite(p_line_W) ? p_line_W * 1e-6f : 0.0f;
        //  Cap keeps the integrator stable for pathological inputs; the cap
        //  is raised while a mitigation radiator is present so the thermal
        //  quench can proceed at full violence.
        pb.P_line_MW = std::min(pb.P_line_MW, f_dm > 0.005f ? 5000.0f : 500.0f);
    }

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
        //  Seed both channels with equal temperatures at cold start.
        float W0 = storedThermalEnergy(ne, state.electron_temp_keV, T_i_keV_, geom) * 1e-6f;
        W_e_MJ_ = 0.5f * W0;
        W_i_MJ_ = 0.5f * W0;
        W_initialized_ = true;
    }
    //  (No re-derivation of W_e from a ramped T_e here — the stored-energy
    //   channels are the single source of truth; temperature is read back
    //   out of them after integration.)

    // ── L-H confinement regime (Martin-2008 threshold + hysteresis) ─────────
    //  The confinement quality genuinely switches: L-mode is roughly half
    //  the IPB98(y,2) H-mode value (H98 ≈ 0.55, i.e. ITER89-P territory);
    //  crossing the Martin scaling latches H-mode, and the plasma only
    //  drops back to L when the heating falls below ~half the threshold
    //  (experimentally observed hysteresis).  This is the make-or-break
    //  physics of every H-mode reactor: you must PAY the L-H power toll
    //  before the confinement needed for burn becomes available.
    float S_area = plasmaSurfaceArea(geom);
    state.P_LH_MW = martinLHThreshold(ne, state.B_toroidal_T, S_area);
    float P_heat_prev = pb.P_aux_e_MW + pb.P_aux_i_MW + pb.P_ohm_MW
                      + state.alpha_power_MW;   // alpha from previous tick
    if (!h_mode_ && P_heat_prev > state.P_LH_MW)        h_mode_ = true;
    else if (h_mode_ && P_heat_prev < 0.5f * state.P_LH_MW) h_mode_ = false;
    float H98 = h_mode_ ? 1.0f : 0.55f;
    state.h_mode     = h_mode_;
    state.H98_factor = H98;

    pb = solvePowerBalance(geom, pb, W_e_MJ_, W_i_MJ_, H98);

    // ── Integrate the two stored-energy channels ─────────────────────────────
    //  Guard against NaN: solvePowerBalance clamps the derivatives to
    //  finite, but double-check here.  Once W becomes NaN, std::max(NaN, 0)
    //  returns NaN on libstdc++ and the poison spreads everywhere.
    if (std::isfinite(pb.dWe_dt_MW)) W_e_MJ_ += pb.dWe_dt_MW * dt;  // MW·s = MJ
    if (std::isfinite(pb.dWi_dt_MW)) W_i_MJ_ += pb.dWi_dt_MW * dt;
    if (!std::isfinite(W_e_MJ_) || W_e_MJ_ < 0.0f) W_e_MJ_ = 0.0f;
    if (!std::isfinite(W_i_MJ_) || W_i_MJ_ < 0.0f) W_i_MJ_ = 0.0f;

    // ── Sawtooth oscillations (q₀ < 1 internal-kink relaxation) ─────────────
    //  When the on-axis safety factor drops below 1, the m=1/n=1 internal
    //  kink periodically flattens the core: temperature ramps up, crashes,
    //  ramps again — the iconic sawtooth trace.  0D proxy: estimate q₀
    //  from the MHD module (previous tick) or q95/3, and when q₀ < 1 run
    //  a crash clock.  Each crash expels a slice of core energy (counted
    //  as a conduction burst) and mixes the channels toward equal T.
    {
        //  q₀ estimate: the MHD module's on-axis cylindrical q — now that
        //  the axis discretisation is fixed, this is directly meaningful
        //  and self-consistent with the Kadomtsev crash proxy (which
        //  resets flat-core q₀ to 1.05 in the same convention).  Fall back
        //  to the q95/3 rule of thumb when the profile is degenerate.
        float q0 = state.q_safety / 3.0f;
        if (mhd_state_ && std::isfinite(mhd_state_->q0)
            && mhd_state_->q0 > 0.001f && mhd_state_->q0 < 50.0f) {
            q0 = mhd_state_->q0;
        }
        state.q0_estimate = q0;
        if (q0 < 1.0f && state.plasma_status == PlasmaStatus::Burning) {
            //  Period: resistive-timescale flavoured — longer when hotter.
            float period = std::clamp(0.5f + 0.075f * state.electron_temp_keV,
                                      0.5f, 4.0f);
            state.sawtooth_period_s = period;
            sawtooth_clock_s_ += dt;
            if (sawtooth_clock_s_ >= period) {
                sawtooth_clock_s_ = 0.0f;
                sawtooth_count_++;
                //  Kadomtsev-ish crash: ~7% of electron and ~4% of ion
                //  stored energy redistributed/lost through the mixing
                //  radius in a single tick.
                W_e_MJ_ *= 0.93f;
                W_i_MJ_ *= 0.96f;
                //  Kadomtsev reconnection: the crash expels core current
                //  OUTWARD until q₀ returns to ~1 — this is the feedback
                //  that regulates the core in a real tokamak.  (Merely
                //  averaging j inside the mixing radius conserves the core
                //  current and leaves q₀ pinned low.)  Set the mixing
                //  region to the flat current density that gives
                //  q₀ = q_target, and hand the expelled current to the
                //  annulus just outside, conserving total I_p exactly.
                if (mhd_state_ && !mhd_state_->j_phi.empty()) {
                    auto& j = mhd_state_->j_phi;
                    const size_t N = j.size();
                    const size_t n_mix = N / 3;
                    if (n_mix > 1 && 2 * n_mix < N) {
                        constexpr float MU0_ = 1.25663706e-6f;
                        constexpr float q_target = 1.05f;
                        //  q on axis for locally flat j (SHAPED units, to
                        //  match the q-profile the MHD module now reports):
                        //    q₀ = F · 2B/(μ₀·R·j₀)  →  j₀ = F·2B/(μ₀·R·q₀)
                        float F = std::clamp(mhd_state_->q_shape_factor, 1.0f, 4.0f);
                        float j_target = F * 2.0f * state.B_toroidal_T
                                       / (MU0_ * geom.R_major_m * q_target);
                        const float dr = geom.a_minor_m / (float)N;
                        float dI = 0.0f;   // expelled current [A]
                        for (size_t k = 0; k < n_mix; k++) {
                            float r = (k + 0.5f) * dr;
                            if (j[k] > j_target) {
                                dI += (j[k] - j_target) * 2.0f * 3.14159265f * r * dr;
                                j[k] = j_target;
                            }
                        }
                        if (dI > 0.0f) {
                            float A_ann = 0.0f;
                            for (size_t k = n_mix; k < 2 * n_mix; k++) {
                                float r = (k + 0.5f) * dr;
                                A_ann += 2.0f * 3.14159265f * r * dr;
                            }
                            if (A_ann > 1e-6f) {
                                float dj = dI / A_ann;
                                for (size_t k = n_mix; k < 2 * n_mix; k++) j[k] += dj;
                            }
                        }
                    }
                }
            }
        } else {
            state.sawtooth_period_s = 0.0f;
            sawtooth_clock_s_ = 0.0f;
        }
        state.sawtooth_count = sawtooth_count_;
    }

    // ── ELMs — Type-I edge-localized modes (H-mode only) ─────────────────────
    //  The H-mode transport barrier builds an edge pedestal holding roughly
    //  a third of the stored energy.  The pedestal pressure gradient
    //  periodically crosses the peeling-ballooning limit and crashes,
    //  expelling a few percent of W onto the divertor in a sub-ms burst.
    //  This is the price of H-mode: the same barrier that doubles your
    //  confinement fires MJ-scale heat pulses at your tungsten tiles.
    //
    //  Mitigation is the operator's job: firing the fuelling pellets FASTER
    //  than the natural ELM frequency paces the pedestal — each pellet
    //  triggers a premature, proportionally smaller crash (ΔW·f ≈ const,
    //  the experimentally observed constant-exhaust scaling).  ITER's
    //  baseline pacing target is ~45 Hz for exactly this reason.
    {
        float W_tot_MJ = W_e_MJ_ + W_i_MJ_;
        bool elm_active = h_mode_ && W_tot_MJ > 1.0f
                       && (state.plasma_status == PlasmaStatus::Burning
                        || state.plasma_status == PlasmaStatus::Initiating);
        if (elm_active) {
            //  Natural Type-I ELM size: ~8% of the pedestal energy
            //  (W_ped ≈ 0.3·W_tot) → ~2.4% of the total stored energy.
            //  At an ITER-grade burn (W ≈ 350 MJ) that's ~8 MJ per crash.
            float W_ped_MJ  = 0.3f * W_tot_MJ;
            float dW_nat_MJ = 0.08f * W_ped_MJ;
            //  Natural frequency from exhaust balance: ELMs carry ~30% of
            //  the power conducted through the separatrix.
            float P_sep = std::max(pb.P_cond_MW, 1.0f);
            float f_nat = std::clamp(0.3f * P_sep / std::max(dW_nat_MJ, 0.1f),
                                     0.2f, 20.0f);
            //  Pellet pacing: pellets arriving faster than f_nat set the ELM
            //  clock; the crash energy shrinks to keep ΔW·f constant.
            float f_pace = state.pellet_frequency_Hz;
            bool  paced  = (f_pace > f_nat) && state.fuel_T_inventory_g > 0.1f;
            float f_elm  = paced ? f_pace : f_nat;
            float dW_MJ  = dW_nat_MJ * (f_nat / f_elm);

            elm_clock_s_ += dt;
            if (elm_clock_s_ >= 1.0f / f_elm) {
                elm_clock_s_ = 0.0f;
                elm_count_++;
                //  Expel the crash energy from both channels in proportion
                //  (the pedestal holds electron and ion pressure alike).
                float fe = (W_tot_MJ > 1e-3f) ? W_e_MJ_ / W_tot_MJ : 0.5f;
                W_e_MJ_ = std::max(W_e_MJ_ - dW_MJ * fe, 0.0f);
                W_i_MJ_ = std::max(W_i_MJ_ - dW_MJ * (1.0f - fe), 0.0f);
                //  ~85% of the expelled energy lands on the divertor strike
                //  points (the rest is radiated in the scrape-off layer).
                state.elm_div_pulse_MJ += 0.85f * dW_MJ;
                state.elm_size_MJ = dW_MJ;
            }
            state.elm_freq_Hz = f_elm;
            state.elm_paced   = paced;
        } else {
            state.elm_freq_Hz = 0.0f;
            state.elm_paced   = false;
            elm_clock_s_      = 0.0f;
        }
        state.elm_count = elm_count_;
    }

    W_stored_MJ_ = W_e_MJ_ + W_i_MJ_;

    // ── Update temperatures from the channels (T_x = (2/3)·W_x/(n·V)) ───────
    float V = plasmaVolume(geom);
    float T_denom = ne * V * 1.602176634e-16f;
    float T_e_new = 0.0f, T_i_new = 0.0f;
    if (T_denom > 1e-30f) {
        T_e_new = (2.0f / 3.0f) * W_e_MJ_ * 1e6f / T_denom;
        T_i_new = (2.0f / 3.0f) * W_i_MJ_ * 1e6f / T_denom;
    }
    if (!std::isfinite(T_e_new)) T_e_new = 0.0f;
    if (!std::isfinite(T_i_new)) T_i_new = 0.0f;
    T_e_new = std::clamp(T_e_new, 0.0f, 200.0f);
    T_i_new = std::clamp(T_i_new, 0.0f, 200.0f);
    state.electron_temp_keV = T_e_new;
    state.plasma_temp_keV   = T_i_new;   // OVERVIEW labels this "Ion Temp"
    state.ion_temp_keV      = T_i_new;
    T_i_keV_                = T_i_new;

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

    // ── Publish the full power-balance ledger for the UI (all MW) ──────────
    auto fin = [](float v){ return std::isfinite(v) ? v : 0.0f; };
    state.P_ohm_MW        = fin(pb.P_ohm_MW);
    state.P_aux_actual_MW = fin(pb.P_aux_MW);
    state.P_brem_MW       = fin(pb.P_brem_MW);
    state.P_sync_MW       = fin(pb.P_sync_MW);
    state.P_line_MW       = fin(pb.P_line_MW);
    state.P_cond_MW       = fin(pb.P_cond_MW);
    state.dW_dt_MW        = fin(pb.dW_dt_MW);

    // ── Safety factor q_95 — Uckan shaping formula (ITER Physics Basis) ────
    //  Includes elongation κ, triangularity δ and the aspect-ratio
    //  correction, so BOTH shape knobs now change the edge safety factor.
    //  ITER-class inputs (a=2, B=5.3, R=6.2, I=15 MA, κ=1.7, δ=0.33)
    //  give q_95 ≈ 3.0 — the design value.
    const float R = geom.R_major_m, a = geom.a_minor_m;
    state.q_safety = q95Uckan(geom, state.delta);

    // ── Greenwald density limit ────────────────────────────────────────────
    //  n_G [10²⁰ m⁻³] = I_p[MA] / (π·a²).  Real tokamaks disrupt via a
    //  radiative edge collapse as the line-averaged density approaches
    //  ~1×n_G; the hard trigger fires slightly above it.
    float n_Greenwald = state.plasma_current_MA / (3.14159265f * a * a);
    state.greenwald_frac = (n_Greenwald > 1e-6f)
                         ? (ne * 1e-20f) / n_Greenwald : 0.0f;
    if (!std::isfinite(state.greenwald_frac)) state.greenwald_frac = 0.0f;
    //  Persistence-based trigger: the density-limit disruption is a
    //  radiative edge collapse, which takes seconds of sustained
    //  over-density to develop — not an instantaneous trip.  This also
    //  forgives the transient at plasma start-up while the current ramp
    //  is still catching up with the density ramp.  A gross excursion
    //  (>1.5×) still disrupts immediately.
    bool gw_established = (state.plasma_current_MA > 3.0f);
    if (gw_established && state.greenwald_frac > 1.05f) {
        gw_violation_s_ += dt;
    } else {
        gw_violation_s_ = std::max(gw_violation_s_ - 2.0f * dt, 0.0f);
    }
    bool greenwald_exceeded = gw_established
                            && ((state.greenwald_frac > 1.5f)
                                || (gw_violation_s_ > 2.0f));

    // ── Troyon β-limit ─────────────────────────────────────────────────────
    constexpr float beta_N_limit = 2.5f;
    //  β = 2·μ₀·n·k(T_e + T_i)/2 / B² — now uses BOTH species temperatures
    //  (the old form used T_e only, under-counting ion pressure whenever
    //  the two channels diverge).
    float T_sum_keV = 0.5f * (state.electron_temp_keV + state.ion_temp_keV) * 2.0f;
    float beta_raw = (state.B_toroidal_T > 0.01f)
                   ? (1.257e-6f * ne * T_sum_keV * 1.602e-16f)
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

    //  (The old ad-hoc H-mode threshold block that lived here has been
    //   superseded by the Martin-2008 scaling with hysteresis computed
    //   BEFORE solvePowerBalance, where it selects the H98 factor.)

    // ── Bootstrap current fraction ─────────────────────────────────────────
    //  Poloidal beta: β_p = β · (B_T/B_p)², with the edge poloidal field
    //  from Ampère's law around the (elongated) poloidal circumference:
    //  B_p = μ₀·I_p / (2π·a·√((1+κ²)/2)).  ITER (15 MA): B_p ≈ 1.1 T →
    //  β_p ≈ 0.65 at the design point.
    //
    //  BUG FIX: this used to divide the total pressure by the ELECTRON
    //  pressure — β·B²/(2μ₀·n·T_e) algebraically reduces to (T_e+T_i)/2T_e
    //  ≈ 1 regardless of the actual plasma, so β_p carried no information
    //  and the bootstrap fraction sat at its clamp.
    float beta_p = 0.0f;
    if (state.plasma_current_MA > 0.05f && state.beta > 1e-6f) {
        float kfac = sqrtf(0.5f * (1.0f + geom.kappa * geom.kappa));
        float B_p  = 1.25663706e-6f * state.plasma_current_MA * 1e6f
                   / (2.0f * 3.14159265f * a * kfac);
        if (B_p > 1e-3f) {
            float ratio = state.B_toroidal_T / B_p;
            beta_p = state.beta * ratio * ratio;
        }
    }
    if (!std::isfinite(beta_p)) beta_p = 0.0f;
    state.q_safety = std::max(state.q_safety, 1.0f);
    float f_bs = bootstrapFraction(a / R, beta_p);

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
        //  Burning gate: current ramped, density at least half of target,
        //  temperature in the fusion-relevant range, and real fusion power.
        //  (The old test compared T_e against the operator's T_e SETPOINT,
        //  which no longer drives the physics — use an absolute threshold.)
        bool ramped_up = state.plasma_current_MA >= 0.9f * state.sp_plasma_current_MA
                      && state.electron_temp_keV >= 5.0f
                      && state.plasma_density_m3 >= 0.5f * state.sp_density_m3;
        if (ramped_up && state.fusion_power_MW > 10.0f) {
            state.plasma_status = PlasmaStatus::Burning;
        }
    } else if (state.plasma_status == PlasmaStatus::Burning) {
        if (state.fusion_power_MW < 5.0f) {
            state.plasma_status = PlasmaStatus::Initiating;
        }
    }
    // ── Quench commitment ────────────────────────────────────────────────────
    //  Any MGI/SPI firing commits the quench immediately (the radiator is
    //  in the plasma — there's no going back).  Sustained disruption flags
    //  commit via the persistence accumulator at the top of this function.
    if (state.dm_injected_impurity > 0.02f && state.plasma_current_MA > 0.5f) {
        quench_committed_ = true;
    }
    if (quench_committed_) {
        // Hold the disruption until the plasma is gone, then land in
        // Quenched — the operator must re-INITIATE for a new discharge.
        // (Without the explicit transition, the setpoint ramps resurrected
        // the plasma the moment the commitment cleared, and a lingering
        // mitigation radiator re-committed it: a zombie limit cycle.)
        state.disruption_flag = true;
        state.plasma_status = PlasmaStatus::Disrupting;
        //  Release only once I_p is below EVERY module's activity floor
        //  (the MHD tracker runs down to 0.01 MA and could re-flag in the
        //  0.01-0.05 MA gap, resurrecting the discharge).
        if (state.plasma_current_MA < 0.008f) {
            quench_committed_ = false;   // plasma gone — commitment served
            disrupt_flag_s_   = 0.0f;
            state.disruption_flag = false;
            state.plasma_status = PlasmaStatus::Quenched;
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
        //  Separatrix floor: a real tokamak edge sits at ~50-100 eV, never
        //  at zero — a zero edge temperature drove the resistivity clamp
        //  (η = 1e9) into the current-diffusion scheme.
        mhd_state_->T_e_keV[i] = std::max(T_e_safe * shape, 0.05f);
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

    //  ── Geometry shaping factor (set BEFORE the MHD update) ────────────────
    //  computeQProfile multiplies the cylindrical q by this Uckan-form
    //  shaping factor (elongation, triangularity, aspect corrections), so
    //  the whole MHD subsystem — resonant surface radii, q0, q95, tearing
    //  analysis — works in REAL q units.  Previously the module ran in
    //  cylindrical q (≈2.6× low for ITER shapes: q95 ≈ 0.9 < the
    //  disruption threshold of 2), and q95 was patched with a κ² factor
    //  after the fact while the resonant surfaces stayed misplaced.
    {
        auto geom = makeGeometry(state);
        float q_cyl = 5.0f * geom.a_minor_m * geom.a_minor_m * state.B_toroidal_T
                    / (std::max(geom.R_major_m, 0.1f)
                       * std::max(state.plasma_current_MA, 0.05f));
        float uck = ConfinementPhysics::q95Uckan(geom, state.delta);
        mhd_state_->q_shape_factor = std::clamp(uck / std::max(q_cyl, 0.01f), 1.0f, 4.0f);
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

    // Write back q95 from the MHD-computed (already shaped) q-profile.
    //  Guard: if the module produced a non-finite or degenerate q95, keep
    //  the 0D Uckan value computed in updatePowerBalance.
    if (std::isfinite(mhd_state_->q95)
        && mhd_state_->q95 > 0.05f && mhd_state_->q95 < 500.0f) {
        state.q_safety = mhd_state_->q95;
    }
}

void PlasmaCoreBridge::readbackDiagnostics(ReactorState& state, float dt)
{
    (void)state; (void)dt;
}
