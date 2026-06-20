//
// confinement_physics.h
// Lawson criterion, IPB98(y,2) H-mode confinement scaling, and 0D power-balance.
//
//  These functions provide the reactor-grade 0D model that complements the
//  kinetic PIC core.  The original PlasmaCoreBridge used a very rough
//  <σv> ≈ 3e-22 · T² approximation and ignored confinement losses entirely;
//  for reactor studies you need the actual ITER98(y,2) scaling plus the
//  full power-balance ODE:
//
//      dW/dt = P_α + P_aux + P_ohm − P_brem − P_sync − P_line − P_cond − P_dia
//
//  where:
//    W        = (3/2) ∫ n (T_e + T_i) dV           (stored thermal energy)
//    P_α      = 0.2 · P_fus                          (α heating)
//    P_aux    = NBI + ICRH + ECH input
//    P_ohm    = η · j²                               (only during ramp-up)
//    P_brem   = C_br · n_e² · Z_eff · √T_e · ḡ_B    (bremsstrahlung)
//    P_sync   = C_s · n_e · B² · T_e^{2.5}           (synchrotron, Albajar)
//    P_line   = n_e · Σ n_i · L_i(T_e)               (line radiation)
//    P_cond   = W / τ_E                              (conduction/convection)
//    P_dia    = particle losses, charge-exchange     (often small in DT)
//
//  τ_E is computed from IPB98(y,2):
//      τ_E [s] = 0.0562 · I_p^0.93 · B_T^0.15 · P_loss^-0.69
//                    · n_19^0.41 · M^0.19 · R^1.97 · ε^0.58 · κ^0.78
//  where:
//    I_p     [MA]    plasma current
//    B_T     [T]     toroidal field on axis
//    P_loss  [MW]    P_aux + P_α − dW/dt  (the loss power going through the edge)
//    n_19    [10¹⁹ m⁻³]  line-averaged density
//    M       [amu]   effective ion mass (2.5 for 50:50 D-T)
//    R       [m]     major radius
//    ε = a/R         inverse aspect ratio
//    κ                elongation
//
//  Ignition occurs when P_α ≥ P_brem + P_sync + P_line + P_cond  (with P_aux=0).
//
//  References:
//    [1] ITER Physics Basis 1999, Nucl. Fusion 39, 2137 (IPB98 scaling).
//    [2] Lawson, Proc. Phys. Soc. B 70, 6 (1957).
//    [3] NRL Plasma Formulary 2023, p. 36 (radiation coefficients).
//    [4] Post, J. Nucl. Mater. 220-222, 143 (1995) (impurity cooling factors).
//

#pragma once
#include <cmath>
#include <algorithm>

namespace ConfinementPhysics {

// ─── Geometry & operating-point bundle ────────────────────────────────────────
struct TokamakGeometry {
    float R_major_m;       // major radius [m]  (ITER: 6.2)
    float a_minor_m;       // minor radius [m]  (ITER: 2.0)
    float kappa;           // elongation (1.7 typical)
    float B_toroidal_T;    // toroidal field on axis [T]
    float I_plasma_MA;     // plasma current [MA]
};

// ─── IPB98(y,2) energy confinement time ──────────────────────────────────────
//
//  τ_E [s] = 0.0562 · I_p^0.93 · B_T^0.15 · P_loss^-0.69
//                 · n_19^0.41 · M^0.19 · R^1.97 · ε^0.58 · κ^0.78
//
//  Validity: H-mode, ITER-like.  H_98 factor multiplies the result; ITER
//  targets H_98 = 1.0-1.1, advanced scenarios achieve up to 1.4.
//
inline float ipb98y2(const TokamakGeometry& g,
                     float P_loss_MW,
                     float n_19,             // [10¹⁹ m⁻³]
                     float M_amu = 2.5f,     // effective ion mass (DT mix)
                     float H_98 = 1.0f)      // confinement enhancement factor
{
    if (P_loss_MW < 0.1f) return 0.0f;
    float eps = g.a_minor_m / g.R_major_m;
    float tau_E = 0.0562f
                * powf(g.I_plasma_MA, 0.93f)
                * powf(g.B_toroidal_T, 0.15f)
                * powf(P_loss_MW,    -0.69f)
                * powf(n_19,          0.41f)
                * powf(M_amu,         0.19f)
                * powf(g.R_major_m,   1.97f)
                * powf(eps,           0.58f)
                * powf(g.kappa,       0.78f);
    return H_98 * tau_E;
}

// ─── Lawson triple product (D-T ignition) ────────────────────────────────────
//
//  Ignition requires:    n · T · τ_E > 3 × 10²¹ keV·s/m³
//  (more precisely, minimum at T ≈ 13-25 keV; factor 1.3 for Z_eff > 1)
//
inline float lawsonTripleProduct(float n_m3, float T_keV, float tau_E_s)
{
    return n_m3 * T_keV * tau_E_s;   // [keV·s/m³]
}

inline bool isIgnited(float n_m3, float T_keV, float tau_E_s, float Z_eff = 1.0f)
{
    float tp = lawsonTripleProduct(n_m3, T_keV, tau_E_s);
    float threshold = 3.0e21f * (1.0f + 0.3f * (Z_eff - 1.0f));
    return tp > threshold;
}

// ─── Stored thermal energy ───────────────────────────────────────────────────
//
//  W = (3/2) · ∫ n · (T_e + T_i) dV    ≈    (3/2) · n · (T_e + T_i) · V_p
//
//  V_p = 2π² · R · a² · κ    (tokamak plasma volume, elongated ellipse)
//
inline float plasmaVolume(const TokamakGeometry& g)
{
    constexpr float PI = 3.14159265358979f;
    return 2.0f * PI * PI * g.R_major_m * g.a_minor_m * g.a_minor_m * g.kappa;
}

inline float storedThermalEnergy(float n_e_m3, float T_e_keV, float T_i_keV,
                                  const TokamakGeometry& g)
{
    float V = plasmaVolume(g);
    // (3/2) n (T_e + T_i) in Joules:  T[keV] · n[m⁻³] · 1.602e-16 J/keV
    return 1.5f * n_e_m3 * (T_e_keV + T_i_keV) * 1.602176634e-16f * V;
}

// ─── Radiation losses (host-side mirrors of radiation.cu) ────────────────────
//
//  These are per-volume power densities [W/m³].  Multiply by V_p for total.
//
inline float bremSStrahlung(float n_e_m3, float T_e_keV, float Z_eff)
{
    if (n_e_m3 < 1.0f || T_e_keV < 0.01f) return 0.0f;
    constexpr float C_br   = 5.34e-37f;
    constexpr float g_avg  = 1.30f;
    constexpr float m_e_c2 = 511.0f;
    float tau = T_e_keV / m_e_c2;
    float corr = 1.0f + 2.76f * tau + 2.04f * tau * tau;
    return C_br * n_e_m3 * n_e_m3 * Z_eff * sqrtf(T_e_keV) * g_avg * corr;
}

inline float synchrotron(float n_e_m3, float T_e_keV, float B_T, float R_wall = 0.6f)
{
    if (n_e_m3 < 1.0f || T_e_keV < 0.1f || B_T < 0.01f) return 0.0f;
    // ── Synchrotron (electron cyclotron) — Trubnikov/Albajar single-pass form.
    //
    //  P_sync [W/m³] = C_s · n_e · B² · T_e^{2.5} · F
    //
    //  IMPORTANT: the constant C_s here is dimensioned for n_e in m⁻³ and
    //  T_e in keV.  The Albajar 2001 fit gives a leading-order coefficient
    //  of ~1.6e-21 W·m³/(T²·keV^{2.5}) after the (T/m_e c²)^{1.5}
    //  prefactor is folded in — NOT 6.21e-17 (which is the *single-pass,
    //  non-reabsorbed, n_e-in-cm⁻³* constant from older NRL editions and
    //  overestimates reactor sync losses by ~4 orders of magnitude).
    //
    //  Sanity check at ITER conditions (n=1e20, B=5.3, T=20 keV, R_w=0.6):
    //      F = 0.4 / (1 + 2.4) = 0.118
    //      P_sync = 1.6e-21 · 1e20 · 28 · 1789 · 0.118
    //             ≈ 9.4e2 W/m³  →  total ≈ 0.8 MW over 840 m³
    //  which matches published ITER estimates (~1-5 MW).
    //
    //  The previous value (6.21e-17) produced ~30 GW of "radiation" at the
    //  same operating point — that was the root cause of the spurious
    //  overtemperature SCRAM seconds after plasma startup.
    constexpr float C_s = 1.6e-21f;
    float F = (1.0f - R_wall) / (1.0f + 0.12f * T_e_keV);
    return C_s * n_e_m3 * B_T * B_T * powf(T_e_keV, 2.5f) * F;
}

// ─── Fusion power density (D-T) ──────────────────────────────────────────────
//
//  P_fus [W/m³] = n_D · n_T · <σv>(T) · E_fus
//                E_fus = 17.59 MeV = 2.819e-12 J
//
//  We compute <σv>(T) via the Bosch-Hale fit on the host (mirroring the
//  CUDA kernel in fusion_data.cuh).  This is the proper Maxwellian-averaged
//  reactivity — NOT the beam-target σ(E) form the original code used.
//
inline float boschHaleSigmaV_DT(float T_keV)
{
    // Bosch-Hale 1992 coefficients for D-T (verified vs. UWFDM-1268 Table 5)
    constexpr float BG    = 34.3827f;
    constexpr float mrc2  = 1124656.0f;
    constexpr float C1    = 1.17302e-9f;
    constexpr float C2    = 1.51361e-2f, C3 = 7.51886e-2f;
    constexpr float C4    = 4.60643e-3f, C5 = 1.35000e-2f;
    constexpr float C6    = -1.06750e-4f, C7 = 1.36600e-5f;

    if (T_keV < 0.2f || T_keV > 100.0f) return 0.0f;

    float P_num = C2 + T_keV * (C4 + T_keV * C6);
    float P_den = 1.0f + T_keV * (C3 + T_keV * (C5 + T_keV * C7));
    float theta_inv = 1.0f - T_keV * P_num / P_den;
    if (std::fabs(theta_inv) < 1e-12f) return 0.0f;
    float theta = T_keV / theta_inv;

    float xi_arg = BG * BG / (4.0f * theta);
    if (xi_arg <= 0.0f) return 0.0f;
    float xi = std::cbrt(xi_arg);          // CUBE ROOT (do NOT omit!)

    float exponent = -3.0f * xi;
    if (exponent < -80.0f) return 0.0f;

    float sv_cm3 = C1 * theta
                 * std::sqrt(xi / (mrc2 * T_keV * T_keV * T_keV))
                 * std::exp(exponent);
    return sv_cm3 * 1e-6f;   // cm³/s → m³/s
}

inline float fusionPowerDensity(float n_D_m3, float n_T_m3, float T_keV)
{
    constexpr float E_fus_J = 17.59e6f * 1.602176634e-19f;
    float sigma_v = boschHaleSigmaV_DT(T_keV);
    return n_D_m3 * n_T_m3 * sigma_v * E_fus_J;
}

// ─── Power balance ODE (single-step RK4) ─────────────────────────────────────
//
//  dW/dt = P_α + P_aux + P_ohm − P_brem − P_sync − P_line − P_cond − P_neutron
//
//  where P_cond = W / τ_E(IPB98(y,2)).  The neutron power leaves the plasma
//  (carries 14.07 MeV out of the plasma volume; ~99% deposited in the blanket).
//
//  The output is the new stored energy W_new and the net heating rate dW/dt.
//  Temperature is then updated by T_new = (2/3) · W_new / (n · V).
//
struct PowerBalance {
    // Inputs (per-step)
    float n_e_m3;
    float T_e_keV;
    float T_i_keV;
    float n_D_m3;
    float n_T_m3;
    float Z_eff;
    float P_aux_MW;
    float P_ohm_MW;

    // Outputs
    float P_fus_MW;
    float P_alpha_MW;
    float P_brem_MW;
    float P_sync_MW;
    float P_line_MW;          // 0 here (impurity radiation handled separately)
    float P_cond_MW;
    float P_neutron_MW;
    float dW_dt_MW;
    float tau_E_s;
    float H_98;
};

inline PowerBalance solvePowerBalance(const TokamakGeometry& g,
                                      const PowerBalance& in,
                                      float W_stored_MJ,
                                      float H_98 = 1.0f)
{
    PowerBalance out = in;
    float V = plasmaVolume(g);
    float n_19 = in.n_e_m3 * 1e-19f;

    // Fusion power (volume-integrated)
    float p_fus_W_m3 = fusionPowerDensity(in.n_D_m3, in.n_T_m3, 0.5f * (in.T_e_keV + in.T_i_keV));
    out.P_fus_MW = p_fus_W_m3 * V * 1e-6f;

    // α heating = 3.5/17.6 · P_fus
    out.P_alpha_MW = out.P_fus_MW * (3.521f / 17.589f);

    // Neutron power leaves the plasma: 14.07/17.59 · P_fus
    out.P_neutron_MW = out.P_fus_MW * (14.070f / 17.589f);

    // Radiation
    float p_brem = bremSStrahlung(in.n_e_m3, in.T_e_keV, in.Z_eff) * V * 1e-6f;
    float p_sync = synchrotron(in.n_e_m3, in.T_e_keV, g.B_toroidal_T) * V * 1e-6f;
    out.P_brem_MW = p_brem;
    out.P_sync_MW = p_sync;
    out.P_line_MW = in.P_line_MW;

    // Loss power for τ_E: P_loss = P_aux + P_α − dW/dt
    // For steady-state dW/dt = 0, so P_loss = P_aux + P_α.
    float P_loss = in.P_aux_MW + out.P_alpha_MW;
    out.tau_E_s = ipb98y2(g, P_loss, n_19, /*M_amu=*/2.5f, H_98);
    out.H_98 = H_98;

    // Conduction loss: P_cond = W / τ_E
    // (W is in MJ; τ_E in s; P_cond in MW)
    //
    // Guard: when I_p → 0 (during SCRAM), ipb98y2 returns tau_E = 0
    // because powf(0, 0.93) = 0.  Without this guard, P_cond = W/0 = inf
    // (or 0/0 = NaN when W is also 0), which propagates through dW/dt
    // into the stored-energy integrator and turns every downstream
    // field (T_e, beta, radiation, fusion power, ...) into NaN.
    if (out.tau_E_s > 1e-6f && W_stored_MJ >= 0.0f && std::isfinite(W_stored_MJ)) {
        out.P_cond_MW = W_stored_MJ / out.tau_E_s;
    } else {
        out.P_cond_MW = 0.0f;
    }

    // dW/dt
    out.dW_dt_MW = in.P_aux_MW + out.P_alpha_MW + in.P_ohm_MW
                 - out.P_brem_MW - out.P_sync_MW - out.P_line_MW
                 - out.P_cond_MW - out.P_neutron_MW;

    // ── NaN / inf safety net ────────────────────────────────────────────────
    //  If any term in the power balance went non-finite (e.g. because a
    //  guarded division above produced inf that then got subtracted from
    //  another inf), clamp the result to 0 rather than letting NaN
    //  propagate into the stored-energy integrator.  Once W_stored_MJ
    //  becomes NaN, std::max(NaN, 0) returns NaN on libstdc++ (because
    //  `NaN < 0` is false, so max returns its first argument), and from
    //  there every downstream field goes NaN permanently.
    if (!std::isfinite(out.dW_dt_MW))  out.dW_dt_MW  = 0.0f;
    if (!std::isfinite(out.P_cond_MW)) out.P_cond_MW = 0.0f;
    if (!std::isfinite(out.P_fus_MW))  out.P_fus_MW  = 0.0f;

    return out;
}

// ─── Q (scientific gain) ──────────────────────────────────────────────────────
//
//  Q = P_fusion / P_aux
//
//  Q = 1: break-even.  Q = ∞: ignition (P_aux = 0).  ITER targets Q = 10.
//
inline float scientificQ(float P_fusion_MW, float P_aux_MW)
{
    if (P_aux_MW < 0.1f) return (P_fusion_MW > 1.0f) ? 1e9f : 0.0f;
    return P_fusion_MW / P_aux_MW;
}

// ─── Bootstrap current fraction (Nevins 14.1) ────────────────────────────────
//
//  f_bs ≈ 0.3 · q · (R/a)^0.5 · β_p
//
//  For advanced tokamak scenarios this can be > 0.5 (steady-state DEMO target).
//
inline float bootstrapFraction(float q_safety, float R_over_a, float beta_p)
{
    return 0.3f * q_safety * sqrtf(R_over_a) * beta_p;
}

} // namespace ConfinementPhysics
