//
// include/mhd_physics.h
// Disruption MHD module: Rutherford island growth + 1D current diffusion.
//
//  Replaces the heuristic disruption flag in PlasmaCoreBridge with two
//  actual MHD models that predict disruption onset from first principles:
//
//    1. Rutherford ODE for magnetic island growth:
//       dW/dt = η' · Δ'(r_s) / μ_0
//       where W = island width, r_s = resonant surface radius, Δ' = tearing
//       mode stability index, η' = resistivity at r_s.
//       When W > W_crit (typically 0.1·a), the island locks to the wall and
//       a major disruption occurs.
//
//    2. 1-D current profile diffusion (cylindrical, axisymmetric):
//       ∂j/∂t = (1/μ_0 r) · ∂/∂r (r · ∂j/∂r) · η(T(r),r)
//       j(r,t) is the toroidal current density, evolved under resistive
//       diffusion.  When the q-profile crosses q=2 or q=3/2 at the edge, a
//       tearing mode is seeded and the Rutherford ODE above takes over.
//
//    3. VDE (Vertical Displacement Event) tracker: when the plasma vertical
//       position z_p drifts beyond the controllability margin (typically
//       κ·a/4), the halo current flows into the first wall and causes a
//       thermal quench.
//
//  References:
//    [1] Rutherford, Phys. Fluids 16, 1903 (1973) — island growth ODE.
//    [2] Wesson, "Tokamaks" 4th ed., Ch. 7 — tearing modes & disruptions.
//    [3] Zakharov, "VDE theory" — 1984 Nucl. Fusion 24, 1003.
//    [4] ITER Physics Basis 1999 §3.4 — disruption loads.
//

#pragma once
#include "ReactorState.h"
#include "confinement_physics.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace MHDPhysics {

// ─── Tokamak 1-D radial mesh ─────────────────────────────────────────────────
struct RadialMesh {
    int   N = 64;                  // number of radial cells
    float a = 2.0f;                // minor radius [m]
    float dr;                      // = a / N

    RadialMesh() : dr(a / N) {}

    float r(int i) const { return (i + 0.5f) * dr; }       // cell centre [m]
    float r_norm(int i) const { return r(i) / a; }         // 0..1
};

// ─── MHD state (carried between time steps) ──────────────────────────────────
struct MHDState {
    // 1-D current density profile j(r) [A/m²]
    std::vector<float> j_phi;
    // 1-D safety factor profile q(r)
    std::vector<float> q_profile;
    // 1-D temperature profile T(r) [keV] (from PlasmaViz or 0D)
    std::vector<float> T_e_keV;
    // Tearing mode state (one per resonant surface q=m/n)
    struct TearingMode {
        float m = 2.0f;            // poloidal mode number (2/1 most common)
        float n = 1.0f;            // toroidal mode number
        float r_s = 0.5f;          // resonant surface radius [m]
        float W = 0.0f;            // island width [m]
        float Delta_prime = 0.0f;  // tearing stability index [1/m]
        bool  locked = false;      // mode-locked to wall?
        bool  disruptive = false;  // W > W_crit?
    };
    TearingMode tm_21;             // q=2/1 mode (most dangerous)
    TearingMode tm_32;             // q=3/2 mode

    // VDE state
    float z_plasma = 0.0f;         // vertical position [m]
    float vz_plasma = 0.0f;        // vertical velocity [m/s]
    float halo_current_MA = 0.0f;  // halo current [MA]

    // Thermal quench + current quench timers
    float t_thermal_quench = -1.0f;  // -1 = not yet triggered
    float t_current_quench = -1.0f;
    float quench_duration_s = 1e-3f; // total disruption timescale

    // Diagnostics
    float q0 = 1.0f;               // on-axis q
    float q95 = 3.0f;              // q at 95% flux
    float li = 0.8f;               // internal inductance
    bool  vde_triggered = false;
    bool  disruption_imminent = false;
    bool  disruption_ongoing = false;

    void resize(int N) {
        j_phi.assign(N, 0.0f);
        q_profile.assign(N, 3.0f);
        T_e_keV.assign(N, 1.0f);
    }
};

// ─── Resistivity: Spitzer (parallel, electron) ───────────────────────────────
//
//  η_∥ = (m_e^{1/2} · e² · ln Λ) / (3 · (2π)^{3/2} · ε_0² · (k_B T_e)^{3/2})
//
//  In practical units:  η [Ω·m] ≈ 1.65e-9 · Z_eff · ln Λ / T_e[keV]^{3/2}
//
inline float spitzerResistivity(float T_e_keV, float Z_eff = 1.0f, float ln_Lambda = 17.0f)
{
    if (T_e_keV < 0.01f) return 1e9f;   // very high resistivity when cold
    return 1.65e-9f * Z_eff * ln_Lambda / std::pow(T_e_keV, 1.5f);
}

// ─── q-profile from current profile (cylindrical approximation) ──────────────
//
//  q(r) = (2π · r² · B_T) / (μ_0 · R · I_enclosed(r))
//  where I_enclosed(r) = ∫_0^r j(r') · 2π · r' dr'
//
inline void computeQProfile(MHDState& mhd, const RadialMesh& mesh,
                             float R_major, float B_T)
{
    float I_enc = 0.0f;
    constexpr float MU0 = 1.25663706e-6f;
    constexpr float PI = 3.14159265358979f;

    for (int i = 0; i < mesh.N; i++) {
        float r = mesh.r(i);
        float dr = mesh.dr;
        // Trapezoidal integration of j·2π·r·dr
        I_enc += mhd.j_phi[i] * 2.0f * PI * r * dr;
        if (I_enc > 1.0f) {
            mhd.q_profile[i] = (2.0f * PI * r * r * B_T) / (MU0 * R_major * I_enc);
        } else {
            mhd.q_profile[i] = 100.0f;
        }
    }

    // On-axis q (q0) and q at 95% flux (q95)
    mhd.q0  = mhd.q_profile[0];
    mhd.q95 = mhd.q_profile[(int)(0.95f * mesh.N)];

    // Internal inductance: li = 2 ∫_0^1 (B_p / B_p,a)^2 · ρ dρ
    //   ≈ 2 · <(B_p/B_p,a)^2>_avg
    float li_sum = 0.0f;
    for (int i = 0; i < mesh.N; i++) {
        float r = mesh.r(i);
        float r_norm = r / mesh.a;
        if (r_norm > 0.01f) {
            // B_p ∝ I_enc/r, B_p,a ∝ I_total/a
            float Bp_ratio = r_norm;  // for parabolic j, this is the rough scaling
            li_sum += Bp_ratio * Bp_ratio * 2.0f * r_norm * mesh.dr / mesh.a;
        }
    }
    mhd.li = std::clamp(li_sum, 0.3f, 2.0f);
}

// ─── Find resonant surface q = m/n ────────────────────────────────────────────
//
//  Returns the radius r_s where q(r) = m/n, or -1 if no such surface exists
//  in the plasma (i.e. q0 > m/n everywhere, or q_a < m/n).
//
inline float findResonantSurface(const std::vector<float>& q_profile,
                                  const RadialMesh& mesh,
                                  float m, float n)
{
    float q_target = m / n;
    for (int i = 1; i < mesh.N; i++) {
        if ((q_profile[i-1] - q_target) * (q_profile[i] - q_target) < 0) {
            // Linear interpolation
            float f = (q_target - q_profile[i-1]) / (q_profile[i] - q_profile[i-1]);
            return mesh.r(i-1) + f * mesh.dr;
        }
    }
    return -1.0f;
}

// ─── Δ' (tearing stability index) from current profile ────────────────────────
//
//  Δ' = [ψ'(r_s⁺) - ψ'(r_s⁻)] / ψ(r_s)
//  where ψ is the perturbed flux function satisfying the tearing-mode
//  eigenvalue equation.  For a step-function current profile j(r) = j_0 (1 - r²/a²),
//  Δ' can be computed analytically; for arbitrary j(r) we use the
//  Fitzpatrick approximation:
//
//    Δ' ≈ (m / r_s) · (1 - 2 · m · ∫_0^{r_s} (j(r')/j(r_s) - 1) (r'/r_s)^{m-1} dr'
//                          - (m+1) ∫_{r_s}^{a} (j(r')/j(r_s) - 1) (r'/r_s)^{m-1} dr')
//
//  A positive Δ' means the tearing mode is UNSTABLE (island grows).
//
inline float tearingDeltaPrime(const std::vector<float>& j_profile,
                                const RadialMesh& mesh,
                                float r_s, int m)
{
    if (r_s < mesh.dr) return 0.0f;
    constexpr float PI = 3.14159265358979f;

    // Find j at r_s (linear interp)
    int i_s = (int)(r_s / mesh.dr);
    if (i_s >= mesh.N - 1) i_s = mesh.N - 2;
    float j_s = j_profile[i_s];
    if (std::fabs(j_s) < 1e-3f) return 0.0f;

    // Inner integral: ∫_0^{r_s} (j/j_s - 1) (r/r_s)^{m-1} dr
    float I_inner = 0.0f;
    for (int i = 0; i < i_s; i++) {
        float r = mesh.r(i);
        float ratio = r / r_s;
        float w = std::pow(ratio, m - 1);
        I_inner += (j_profile[i] / j_s - 1.0f) * w * mesh.dr;
    }
    // Outer integral: ∫_{r_s}^{a} (j/j_s - 1) (r/r_s)^{m-1} dr
    float I_outer = 0.0f;
    for (int i = i_s + 1; i < mesh.N; i++) {
        float r = mesh.r(i);
        float ratio = r / r_s;
        if (ratio < 1.0f) continue;
        float w = std::pow(ratio, m - 1);
        I_outer += (j_profile[i] / j_s - 1.0f) * w * mesh.dr;
    }

    float Delta_prime = ((float)m / r_s) * (1.0f - 2.0f * (float)m * I_inner
                       - (float)(m + 1) * I_outer);
    return Delta_prime;
}

// ─── Rutherford ODE: dW/dt = (η/μ_0) · Δ' ────────────────────────────────────
//
//  The classical Rutherford equation for island width evolution.  When Δ' > 0
//  (unstable tearing), W grows linearly in time.  When W exceeds W_crit
//  (~10% of minor radius), the island is large enough to flatten the current
//  profile locally → seeds a chain reaction that ends in disruption.
//
inline void updateTearingMode(MHDState::TearingMode& tm,
                              const MHDState& mhd,
                              const RadialMesh& mesh,
                              float dt)
{
    if (tm.r_s < 0.0f) {
        tm.W = 0.0f;
        return;
    }

    // Local resistivity at r_s (Spitzer, using T_e at that radius)
    int i_s = (int)(tm.r_s / mesh.dr);
    if (i_s >= mesh.N) i_s = mesh.N - 1;
    float T_local = (i_s < (int)mhd.T_e_keV.size()) ? mhd.T_e_keV[i_s] : 1.0f;
    float eta = spitzerResistivity(T_local);

    constexpr float MU0 = 1.25663706e-6f;

    // Recompute Δ' (it changes as j(r) evolves)
    tm.Delta_prime = tearingDeltaPrime(mhd.j_phi, mesh, tm.r_s, (int)tm.m);

    // dW/dt = (η/μ_0) · Δ'    [Rutherford 1973]
    // Note: this is the LINEAR island width; nonlinear saturation is
    // approximated by capping W at 0.4·a.
    if (tm.Delta_prime > 0.0f) {
        tm.W += (eta / MU0) * tm.Delta_prime * dt;
    } else {
        // Stable: island decays
        tm.W *= std::exp(-dt / 0.1f);   // 100ms decay timescale
    }

    // Cap at 40% of minor radius (saturation)
    tm.W = std::min(tm.W, 0.4f * mesh.a);

    // Disruption criterion: W > 0.1·a AND not yet disruptive
    bool was_disruptive = tm.disruptive;
    tm.disruptive = (tm.W > 0.1f * mesh.a);
    if (tm.disruptive && !was_disruptive) {
        // Mode has just become disruptive — flag it
        tm.locked = true;
    }
}

// ─── 1-D current diffusion ───────────────────────────────────────────────────
//
//  ∂j/∂t = (1/μ_0) · (1/r) · ∂/∂r (r · η · ∂j/∂r)
//
//  Discretized as implicit Crank-Nicolson for stability:
//    j_i^{n+1} = j_i^n + (Δt/(μ_0 · Δr²)) · [η_{i+1/2}(j_{i+1}^{n+1} - j_i^{n+1})
//                                            - η_{i-1/2}(j_i^{n+1} - j_{i-1}^{n+1})]
//  We use a simple explicit forward-Euler step here for clarity, which is
//  stable for Δt < Δr²·μ_0/(2·η_max).
//
inline void diffuseCurrent(MHDState& mhd, const RadialMesh& mesh, float dt)
{
    constexpr float MU0 = 1.25663706e-6f;

    std::vector<float> j_new = mhd.j_phi;

    for (int i = 1; i < mesh.N - 1; i++) {
        float r = mesh.r(i);
        float eta_i   = spitzerResistivity(mhd.T_e_keV[i]);
        float eta_ip  = spitzerResistivity(mhd.T_e_keV[i+1]);
        float eta_im  = spitzerResistivity(mhd.T_e_keV[i-1]);
        float eta_half_p = 0.5f * (eta_i + eta_ip);
        float eta_half_m = 0.5f * (eta_i + eta_im);

        // (1/r) ∂/∂r (r · η · ∂j/∂r)  ≈  (1/r_i) · (r_{i+1/2}·η_{i+1/2}·(j_{i+1}-j_i)
        //                                  - r_{i-1/2}·η_{i-1/2}·(j_i-j_{i-1})) / Δr²
        float r_p = r + 0.5f * mesh.dr;
        float r_m = r - 0.5f * mesh.dr;
        float flux_p = r_p * eta_half_p * (mhd.j_phi[i+1] - mhd.j_phi[i]) / mesh.dr;
        float flux_m = r_m * eta_half_m * (mhd.j_phi[i] - mhd.j_phi[i-1]) / mesh.dr;
        float djdt = (1.0f / (MU0 * r)) * (flux_p - flux_m) / mesh.dr;

        j_new[i] = mhd.j_phi[i] + djdt * dt;

        // Stability check: if dT/dt > 10·j/dt, the timestep is too large
        // (shouldn't happen with proper sub-cycling in the bridge)
    }

    // Boundary conditions:
    //  - r=0:  ∂j/∂r = 0  (axis symmetry)  →  j[0] = j[1]
    //  - r=a:  j[a] = 0   (edge current carried by bootstrap/edge pedestal)
    j_new[0] = j_new[1];
    j_new[mesh.N - 1] = 0.0f;

    // Renormalize to preserve total current I_p (the diffusion can drift the
    // total due to discretization error)
    float I_total_new = 0.0f;
    float I_total_old = 0.0f;
    constexpr float PI = 3.14159265358979f;
    for (int i = 0; i < mesh.N; i++) {
        float r = mesh.r(i);
        I_total_new += j_new[i] * 2.0f * PI * r * mesh.dr;
        I_total_old += mhd.j_phi[i] * 2.0f * PI * r * mesh.dr;
    }
    if (std::fabs(I_total_new) > 1.0f) {
        float scale = I_total_old / I_total_new;
        for (int i = 0; i < mesh.N; i++) j_new[i] *= scale;
    }

    mhd.j_phi = std::move(j_new);
}

// ─── VDE (Vertical Displacement Event) tracker ───────────────────────────────
//
//  Simple 1-D vertical equation of motion:
//    m_p · d²z/dt² = F_vertical_stability + F_halo_drag
//  where:
//    F_z = -k_v · z  (vertical field gradient; k_v > 0 means stable, k_v < 0 unstable)
//    k_v = (μ_0 I_p² / R · a) · (3/2 · (1 + κ²) · n - 1)   (n = decay index)
//
//  When |z| > κ·a/4, the plasma hits the wall and halo currents flow.
//
struct VDEParams {
    float I_p_MA = 15.0f;
    float R_major = 6.2f;
    float a_minor = 2.0f;
    float kappa = 1.7f;
    //  Vertical stability.  n_decay_index > 0 means the external field
    //  configuration is vertically STABLE (the default for a well-controlled
    //  tokamak with active vertical stabilization).  n_decay_index < 0
    //  means unstable — the plasma will VDE on a millisecond timescale
    //  unless the VS system is actively counteracting it.
    //
    //  The old default of -0.5 caused the VDE to trigger within the first
    //  few simulation ticks: the vertical instability growth rate is
    //  ω ≈ sqrt(|k_v|/m_p) ≈ 3e5 rad/s, but the game tick is 1 ms, so z
    //  grows by exp(300) per tick → instant overflow to inf → NaN in the
    //  VDE integrator → latched vde_triggered = true forever.
    //
    //  Default is now +0.5 (stable).  The operator can trigger a VDE
    //  scenario by setting this negative externally.
    float n_decay_index = 0.5f;
    float halo_drag_coeff = 1.0f;  // s^-1
    //  Active vertical stabilization (VS) gain.  The VS coils provide a
    //  restoring force proportional to z.  With k_vs > |k_v|, the plasma
    //  is stable even if n_decay_index < 0.  This represents the fast
    //  feedback system that all modern tokamaks use.
    float k_vs = 0.0f;  // 0 = VS off (use n_decay_index to determine stability)
};

inline void updateVDE(MHDState& mhd, const VDEParams& p, float dt)
{
    constexpr float MU0 = 1.25663706e-6f;
    constexpr float PI  = 3.14159265358979f;

    // Vertical stability parameter (negative = unstable)
    // k_v ≈ (μ_0 · I_p² / (R · a)) · (3(1+κ²)n/2 - 1)
    float I_p = p.I_p_MA * 1e6f;
    float k_v = (MU0 * I_p * I_p / (p.R_major * p.a_minor))
              * (1.5f * (1.0f + p.kappa * p.kappa) * p.n_decay_index - 1.0f);

    // Active vertical stabilization: adds a restoring stiffness k_vs
    float k_total = k_v + p.k_vs;   // net stiffness (positive = stable)

    // Effective plasma mass (rough: m_p ≈ ρ · V_p; for DT plasma ρ ~ 1e-6 kg/m³)
    float V_p = 2.0f * PI * PI * p.R_major * p.a_minor * p.a_minor * p.kappa;
    float m_p = 1e-6f * V_p;   // ~1e-4 kg (very light — hence very fast VDE!)

    // Equation of motion: m_p · d²z/dt² = -k_total · z - drag · v
    //  Guard against NaN/inf in z_plasma from a prior tick's overflow.
    //  If z is already non-finite, the plasma has already "hit the wall"
    //  — just latch vde_triggered and skip the integration.
    float z_safe = std::isfinite(mhd.z_plasma) ? mhd.z_plasma : 0.0f;
    float vz_safe = std::isfinite(mhd.vz_plasma) ? mhd.vz_plasma : 0.0f;

    float F_restore = -k_total * z_safe;
    float F_drag = -p.halo_drag_coeff * m_p * vz_safe;
    float az = (F_restore + F_drag) / m_p;

    // Guard against NaN acceleration (from k_total being NaN etc.)
    if (!std::isfinite(az)) az = 0.0f;

    mhd.vz_plasma = vz_safe + az * dt;
    mhd.z_plasma  = z_safe + mhd.vz_plasma * dt;

    // Cap z to prevent numerical overflow.  The VDE growth rate can be
    // ~300,000 rad/s, but the timestep is 1 ms — so z can grow by
    // exp(300) per tick, which overflows to inf and then NaN.  Once
    // |z| exceeds z_max, the plasma has hit the wall: latch the VDE,
    // clamp z to z_max, and zero the velocity.  This prevents the
    // inf/NaN cascade while still correctly flagging the disruption.
    float z_max = p.kappa * p.a_minor * 0.25f;
    if (std::fabs(mhd.z_plasma) > z_max || !std::isfinite(mhd.z_plasma)) {
        mhd.z_plasma = std::copysign(z_max, mhd.z_plasma);
        mhd.vz_plasma = 0.0f;
        mhd.vde_triggered = true;
        mhd.halo_current_MA = p.I_p_MA * 0.5f;
    }
}

// ─── Top-level MHD update (called once per game tick) ────────────────────────
//
//  Combines:
//    1. Current profile diffusion (slow timescale ~ seconds)
//    2. Tearing mode Rutherford growth (fast, ~100ms)
//    3. VDE tracking (fast, ~10ms)
//    4. Disruption cascade logic
//
struct MHDEventLog {
    bool  tearing_21_unstable = false;
    bool  tearing_32_unstable = false;
    bool  tearing_21_disruptive = false;
    bool  tearing_32_disruptive = false;
    bool  vde_triggered = false;
    bool  thermal_quench = false;
    bool  current_quench = false;
    float W_21_m = 0.0f;
    float W_32_m = 0.0f;
    float Delta_prime_21 = 0.0f;
    float Delta_prime_32 = 0.0f;
};

inline MHDEventLog updateMHD(MHDState& mhd, const RadialMesh& mesh,
                              float I_p_MA, float B_T, float dt)
{
    MHDEventLog log{};

    // 1. Current diffusion
    diffuseCurrent(mhd, mesh, dt);

    // 2. Recompute q-profile (uses current profile + geometry)
    computeQProfile(mhd, mesh, /*R_major=*/6.2f, B_T);

    // 3. Find resonant surfaces and update tearing modes
    mhd.tm_21.r_s = findResonantSurface(mhd.q_profile, mesh, /*m=*/2.0f, /*n=*/1.0f);
    mhd.tm_32.r_s = findResonantSurface(mhd.q_profile, mesh, /*m=*/3.0f, /*n=*/2.0f);

    updateTearingMode(mhd.tm_21, mhd, mesh, dt);
    updateTearingMode(mhd.tm_32, mhd, mesh, dt);

    log.W_21_m = mhd.tm_21.W;
    log.W_32_m = mhd.tm_32.W;
    log.Delta_prime_21 = mhd.tm_21.Delta_prime;
    log.Delta_prime_32 = mhd.tm_32.Delta_prime;
    log.tearing_21_unstable = (mhd.tm_21.Delta_prime > 0.0f);
    log.tearing_32_unstable = (mhd.tm_32.Delta_prime > 0.0f);
    log.tearing_21_disruptive = mhd.tm_21.disruptive;
    log.tearing_32_disruptive = mhd.tm_32.disruptive;

    // 4. VDE tracking (only meaningful if I_p > 0.1 MA)
    if (I_p_MA > 0.1f) {
        VDEParams vp;
        vp.I_p_MA = I_p_MA;
        updateVDE(mhd, vp, dt);
        log.vde_triggered = mhd.vde_triggered;
    }

    // 5. Disruption cascade logic
    bool any_disruptive = mhd.tm_21.disruptive || mhd.tm_32.disruptive;
    bool any_vde = mhd.vde_triggered;

    mhd.disruption_imminent = any_disruptive || any_vde;

    if (mhd.disruption_imminent && mhd.t_thermal_quench < 0.0f) {
        // Thermal quench begins: T drops to ~10 eV in ~1 ms
        mhd.t_thermal_quench = 0.0f;
        log.thermal_quench = true;
    }

    if (mhd.t_thermal_quench >= 0.0f) {
        mhd.t_thermal_quench += dt;
        if (mhd.t_thermal_quench > 1e-3f && mhd.t_current_quench < 0.0f) {
            mhd.t_current_quench = 0.0f;
            log.current_quench = true;
        }
        if (mhd.t_current_quench >= 0.0f) {
            mhd.t_current_quench += dt;
            mhd.disruption_ongoing = (mhd.t_current_quench < 0.1f);  // 100ms total
        }
    }

    return log;
}

// ─── Initialize the MHD state from a 0D ReactorState ─────────────────────────
//
//  Sets the initial current profile to a parabolic shape
//  j(r) = j_0 · (1 - (r/a)²)^α_j  with α_j = 1.0 (parabolic, typical for inductive drive)
//  and q-profile accordingly.
//
inline void initializeMHD(MHDState& mhd, const RadialMesh& mesh,
                           float I_p_MA, float B_T)
{
    mhd.resize(mesh.N);

    // Parabolic current profile
    constexpr float PI = 3.14159265358979f;
    float I_p = I_p_MA * 1e6f;
    // Total current: I = ∫_0^a j_0 (1 - r²/a²)^α · 2π r dr
    // For α=1:  I = j_0 · π · a²/2  →  j_0 = 2I / (π a²)
    float j_0 = 2.0f * I_p / (PI * mesh.a * mesh.a);
    for (int i = 0; i < mesh.N; i++) {
        float rn = mesh.r_norm(i);
        mhd.j_phi[i] = j_0 * std::max(1.0f - rn * rn, 0.0f);
        // Initial temperature profile (will be updated by PlasmaViz data)
        mhd.T_e_keV[i] = 10.0f * std::max(1.0f - rn * rn, 0.0f);
    }

    computeQProfile(mhd, mesh, /*R_major=*/6.2f, B_T);

    // Find resonant surfaces
    mhd.tm_21.r_s = findResonantSurface(mhd.q_profile, mesh, 2.0f, 1.0f);
    mhd.tm_32.r_s = findResonantSurface(mhd.q_profile, mesh, 3.0f, 2.0f);

    mhd.tm_21.W = 0.0f;
    mhd.tm_32.W = 0.0f;
    mhd.disruption_imminent = false;
    mhd.disruption_ongoing = false;
    mhd.t_thermal_quench = -1.0f;
    mhd.t_current_quench = -1.0f;
}

} // namespace MHDPhysics
