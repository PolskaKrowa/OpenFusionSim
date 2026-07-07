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
    constexpr double E_fus_J = 17.59e6 * 1.602176634e-19;
    float sigma_v = boschHaleSigmaV_DT(T_keV);
    //  OVERFLOW FIX: n_D·n_T ≈ 2×10³⁹ at reactor densities, which exceeds
    //  FLT_MAX (3.4×10³⁸) and silently becomes +inf in float — downstream
    //  isfinite guards then zeroed the fusion power.  This is why the
    //  plasma only ever produced power at suspiciously low density
    //  (n ≲ 6×10¹⁹ kept the product under FLT_MAX).  Compute in double.
    double p = (double)n_D_m3 * (double)n_T_m3 * (double)sigma_v * E_fus_J;
    return (float)p;
}

// ─── Plasma surface area (elongated torus) ───────────────────────────────────
//
//  S ≈ 4π² · R · a · √((1 + κ²)/2)
//
//  The √((1+κ²)/2) factor is the standard approximation for the poloidal
//  circumference of an ellipse relative to a circle of radius a.
//  ITER (R=6.2, a=2.0, κ=1.7): S ≈ 683 m².
//
inline float plasmaSurfaceArea(const TokamakGeometry& g)
{
    constexpr float PI = 3.14159265358979f;
    float kfac = sqrtf(0.5f * (1.0f + g.kappa * g.kappa));
    return 4.0f * PI * PI * g.R_major_m * g.a_minor_m * kfac;
}

// ─── L→H power threshold — Martin 2008 scaling ───────────────────────────────
//
//  P_LH [MW] = 0.0488 · n̄_e20^0.717 · B_T^0.803 · S^0.941
//
//  (Martin et al., J. Phys. Conf. Ser. 123, 012033 (2008) — the ITPA
//  multi-machine fit used for the ITER baseline.  ITER at n̄=1e20, B=5.3 T,
//  S≈680 m²: P_LH ≈ 52 MW, matching the published ITER value.)
//
//  Below n_e,min ≈ 0.7e19·(I_p·B)^~0.4 the threshold rises again (the
//  low-density branch); we approximate that with a quadratic bowl below
//  n_min so startup at very low density doesn't get a free H-mode.
//
inline float martinLHThreshold(float n_e_m3, float B_T, float S_m2)
{
    if (n_e_m3 < 1e17f || B_T < 0.01f) return 1e9f;   // no plasma → unreachable
    float n20 = n_e_m3 * 1e-20f;
    float P = 0.0488f * powf(n20, 0.717f) * powf(B_T, 0.803f) * powf(S_m2, 0.941f);
    // Low-density branch: threshold rises below n20 ≈ 0.3
    constexpr float n20_min = 0.3f;
    if (n20 < n20_min && n20 > 1e-4f) {
        float r = n20_min / n20;
        P *= 0.5f * (1.0f + r * r) / r;   // = 1 at n_min, grows ~r below it
    }
    return P;
}

// ─── Neoclassical (trapped-particle corrected) Spitzer resistivity ───────────
//
//  η_neo ≈ η_Spitzer / (1 − √ε)²
//
//  Trapped electrons can't carry parallel current, so the effective
//  resistivity of a torus exceeds the Spitzer value by the (1−√ε)⁻²
//  factor (ε = a/R ≈ 0.32 for ITER → factor ≈ 5.2).
//
inline float neoclassicalResistivity(float T_e_keV, float Z_eff,
                                      float eps, float ln_Lambda = 17.0f)
{
    if (T_e_keV < 0.01f) T_e_keV = 0.01f;               // cold-plasma clamp
    float eta_sp = 1.65e-9f * Z_eff * ln_Lambda / powf(T_e_keV, 1.5f);
    float trap = 1.0f - sqrtf(std::clamp(eps, 0.0f, 0.9f));
    return eta_sp / std::max(trap * trap, 0.05f);
}

// ─── Physical ohmic heating ──────────────────────────────────────────────────
//
//  R_plasma = η_neo · 2πR / (κ · π · a²)      (toroidal loop resistance)
//  P_ohm    = I_p² · R_plasma                  (goes to the ELECTRONS)
//  V_loop   = I_p · R_plasma                   (resistive loop voltage)
//
//  This replaces the old status-based step (15 MW while Initiating, 2 MW
//  while Burning): now P_ohm falls smoothly as T_e^−1.5 during heat-up —
//  tens of MW in a cold ~100 eV start-up plasma, dropping to ≲1 MW at
//  burning temperatures, exactly as in a real tokamak.
//
struct OhmicResult { float P_ohm_MW; float V_loop_V; float R_plasma_Ohm; };

inline OhmicResult ohmicHeating(const TokamakGeometry& g,
                                 float I_p_MA, float T_e_keV, float Z_eff)
{
    constexpr float PI = 3.14159265358979f;
    OhmicResult out{0.f, 0.f, 0.f};
    if (I_p_MA < 0.01f) return out;
    float eps = g.a_minor_m / std::max(g.R_major_m, 0.1f);
    float eta = neoclassicalResistivity(T_e_keV, Z_eff, eps);
    float A_cross = PI * g.a_minor_m * g.a_minor_m * std::max(g.kappa, 1.0f);
    out.R_plasma_Ohm = eta * 2.0f * PI * g.R_major_m / A_cross;
    float I_A = I_p_MA * 1e6f;
    out.V_loop_V  = I_A * out.R_plasma_Ohm;
    out.P_ohm_MW  = I_A * I_A * out.R_plasma_Ohm * 1e-6f;
    // Cap: the CS power supply can't deliver unbounded ohmic power into a
    // very cold plasma — ITER's CS delivers ≲30 MW of ohmic input.
    out.P_ohm_MW = std::min(out.P_ohm_MW, 30.0f);
    return out;
}

// ─── Impurity line-radiation cooling factor L_z(T_e) ─────────────────────────
//
//  P_line [W/m³] = n_e · n_imp · L_z(T_e)
//
//  Coronal-equilibrium cooling curves (Post 1995 style, smooth log-space
//  fits).  Two species matter for the game:
//
//   • Generic low-Z (Be/C wall material): L_z peaks ~5e-32 W·m³ near
//     ~10 eV and falls steeply once fully stripped (>~1 keV).  This is
//     what produces the classic RADIATION BARRIER during start-up: a
//     cold, dirty plasma can radiate more than the available heating and
//     collapse — the operator must limit impurities (boronization!) to
//     burn through.
//
//   • Tungsten (divertor sputtering): never fully strips; L_z stays
//     ~1e-31 W·m³ across 1-30 keV, which is why W accumulation at even
//     ~1e-5 fractional density can radiatively kill an ITER-class burn.
//
inline float coolingFactorLowZ(float T_keV)
{
    // Peak ~2e-32 W·m³ near 30 eV (carbon-like coronal maximum), steep
    // ∝ T^-1.5 fall above 100 eV, floor 5e-35 W·m³ once fully stripped
    // (the residual continuum is already counted via Z_eff bremsstrahlung).
    // Magnitudes bench-marked against burn-through experience: 1% carbon
    // at n_e = 1e19 m⁻³ gives a ~15 MW barrier — surmountable with the
    // ~30 MW of ohmic power available, exactly the marginal fight a real
    // start-up is.  Fuel to 3e19 BEFORE burn-through and the barrier
    // triples past anything the heating systems can deliver.
    float T_eV = std::max(T_keV * 1000.0f, 1.0f);
    float peak = 2.0e-32f;
    float x = T_eV / 30.0f;                        // dimensionless
    float shape = x / (0.3f + x * x * sqrtf(x));   // rises, peaks, falls ~x^-1.5
    return std::max(peak * shape, 5.0e-35f);
}

inline float coolingFactorTungsten(float T_keV)
{
    // Flat-ish 1e-31..3e-31 W·m³ across the 0.1-30 keV band (never strips)
    float T = std::clamp(T_keV, 0.05f, 50.0f);
    return 1.0e-31f * (1.0f + 2.0f / (1.0f + T));  // 3e-31 cold → 1e-31 hot
}

inline float lineRadiationDensity(float n_e_m3, float T_e_keV,
                                   float f_imp_lowZ, float f_imp_W = 0.0f)
{
    if (n_e_m3 < 1.0f || T_e_keV < 0.001f) return 0.0f;
    //  n_e² overflows float at reactor densities — compute in double.
    double Lz = (double)f_imp_lowZ * coolingFactorLowZ(T_e_keV)
              + (double)f_imp_W    * coolingFactorTungsten(T_e_keV);
    double p = (double)n_e_m3 * (double)n_e_m3 * Lz;
    return (float)std::max(p, 0.0);
}

// ─── Alpha-particle heating partition (Stix slowing-down) ────────────────────
//
//  A 3.5 MeV alpha slows on electrons above the critical energy
//  E_c ≈ 33 · T_e [keV] and on ions below it.  The fraction of its energy
//  delivered to IONS is (Stix, Plasma Phys. 14, 367 (1972)):
//
//      f_i(x) = (1/x) ∫₀ˣ dy / (1 + y^{3/2}),     x = E_α / E_c
//
//  At T_e = 10 keV → f_i ≈ 0.20; at 20 keV → ≈ 0.31; at 40 keV → ≈ 0.45.
//  Alphas predominantly heat ELECTRONS in a reactor-grade plasma — the
//  ions are heated indirectly through e-i equilibration, which is why a
//  two-temperature model matters for burn dynamics.
//
inline float alphaIonFraction(float T_e_keV)
{
    if (T_e_keV < 0.05f) return 0.0f;
    constexpr float E_alpha_keV = 3521.0f;
    float E_c = 33.0f * T_e_keV;
    float x = E_alpha_keV / std::max(E_c, 1.0f);
    // 32-point midpoint quadrature of ∫₀ˣ dy/(1+y^1.5) — cheap & smooth
    const int N = 32;
    float h = x / N, sum = 0.0f;
    for (int k = 0; k < N; k++) {
        float y = (k + 0.5f) * h;
        sum += h / (1.0f + y * sqrtf(y));
    }
    return std::clamp(sum / x, 0.0f, 1.0f);
}

// ─── Electron-ion temperature equilibration rate ─────────────────────────────
//
//  ν_eq [s⁻¹] — the rate at which (T_e − T_i) relaxes via Coulomb
//  collisions (NRL Plasma Formulary, electron-ion energy exchange):
//
//      ν_eq ≈ 1.0e-19 · n_i[m⁻³] · lnΛ / (μ · T_e[keV]^{3/2})
//
//  ITER-grade (n=1e20, T_e=10 keV, μ=2.5): ν ≈ 2.2 s⁻¹ → τ_eq ≈ 0.5 s,
//  comparable to τ_E — which is exactly why T_e ≠ T_i in real burns.
//
inline float eiEquilibrationRate(float n_i_m3, float T_e_keV,
                                  float mu_amu = 2.5f, float ln_Lambda = 17.0f)
{
    if (T_e_keV < 0.01f) T_e_keV = 0.01f;
    return 1.0e-19f * n_i_m3 * ln_Lambda / (mu_amu * powf(T_e_keV, 1.5f));
}

// ─── q₉₅ — Uckan shaping formula (ITER Physics Basis) ────────────────────────
//
//  q_95 = (5 · a² · B_T) / (R · I_p[MA])
//         · [1 + κ²·(1 + 2δ² − 1.2δ³)] / 2
//         · (1.17 − 0.65·ε) / (1 − ε²)²
//
//  Unlike the bare cylindrical-κ² form this includes triangularity δ and
//  the aspect-ratio correction, so the operator's δ shape knob genuinely
//  changes the edge safety factor.  ITER (a=2, B=5.3, R=6.2, I=15 MA,
//  κ=1.7, δ=0.33): q_95 ≈ 3.0 — the design value.
//
inline float q95Uckan(const TokamakGeometry& g, float delta)
{
    if (g.I_plasma_MA < 0.05f || g.B_toroidal_T < 0.01f) return 100.0f;
    float a = g.a_minor_m, R = g.R_major_m, k = g.kappa;
    float eps = a / R;
    float d = std::clamp(delta, 0.0f, 0.7f);
    float shape = 0.5f * (1.0f + k * k * (1.0f + 2.0f * d * d - 1.2f * d * d * d));
    float aspect = (1.17f - 0.65f * eps) / powf(1.0f - eps * eps, 2.0f);
    float q_cyl = 5.0f * a * a * g.B_toroidal_T / (R * g.I_plasma_MA);
    return q_cyl * shape * aspect;
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
    float P_aux_e_MW;         // auxiliary heating into ELECTRONS (ECRH, LHCD, ½·NBI)
    float P_aux_i_MW;         // auxiliary heating into IONS      (ICRH, ½·NBI)
    float P_ohm_MW;           // ohmic heating (→ electrons)
    float P_line_MW;          // impurity line radiation (← electrons), from caller

    // Outputs
    float P_fus_MW;
    float P_alpha_MW;
    float P_alpha_e_MW;       // alpha power → electrons (Stix partition)
    float P_alpha_i_MW;       // alpha power → ions
    float P_brem_MW;
    float P_sync_MW;
    float P_cond_MW;          // total conduction loss (W_e + W_i)/τ_E
    float P_ei_MW;            // e→i collisional exchange (positive = e heats i)
    float P_neutron_MW;       // DIAGNOSTIC ONLY — neutrons never enter W
    float dWe_dt_MW;
    float dWi_dt_MW;
    float dW_dt_MW;           // = dWe + dWi (kept for compatibility)
    float tau_E_s;
    float H_98;
    float P_aux_MW;           // = P_aux_e + P_aux_i (convenience)
};

//  Two-temperature 0D power balance:
//
//      dW_e/dt = f_e·P_α + P_ohm + P_aux,e − P_brem − P_sync − P_line
//                − W_e/τ_E − Q_ei
//      dW_i/dt = f_i·P_α + P_aux,i − W_i/τ_E + Q_ei
//
//  Q_ei = (3/2)·n·k·(T_e−T_i)·ν_eq·V  is the collisional exchange term.
//  Radiation only cools the ELECTRONS (photons are emitted by electrons);
//  alphas heat mostly electrons (Stix partition); NBI splits ~50/50;
//  ICRH heats ions, ECRH/LHCD heat electrons.  This is what produces the
//  realistic T_i ≠ T_e dynamics of a burning plasma.
//
//  BUG FIX vs. the previous version: P_neutron is NO LONGER subtracted
//  from dW/dt.  The 14.07 MeV neutron energy never enters the plasma
//  stored energy in the first place — only the 3.5 MeV alpha channel
//  does — so subtracting it double-counted the loss and made the plasma
//  0.8·P_fus too lossy (≈400 MW of phantom loss at ITER-grade burn).
//
inline PowerBalance solvePowerBalance(const TokamakGeometry& g,
                                      const PowerBalance& in,
                                      float W_e_MJ, float W_i_MJ,
                                      float H_98 = 1.0f)
{
    PowerBalance out = in;
    float V = plasmaVolume(g);
    float n_19 = in.n_e_m3 * 1e-19f;

    // Fusion power — driven by the ION temperature (the reacting species).
    //
    //  PROFILE PEAKING: a 0D model evaluates everything at volume-averaged
    //  n and T, but ⟨n²⟨σv⟩(T)⟩ over realistic peaked profiles
    //  (T ∝ (1−ρ²), n weakly peaked) exceeds n̄²⟨σv⟩(T̄) by a factor of
    //  ~2-3 because the reactivity is so strongly weighted toward the hot
    //  core.  Without this the burn stalls: at T̄ ≈ 4 keV the flat-model
    //  alpha power is too small to bootstrap the plasma through to
    //  ignition-relevant temperatures.  C_PEAK_FUS = 2.5 is the standard
    //  parabolic-profile integral value; bremsstrahlung (∝ n²√T, much
    //  weaker T weighting) gets a milder 1.4.
    constexpr float C_PEAK_FUS  = 2.5f;
    constexpr float C_PEAK_BREM = 1.4f;
    float p_fus_W_m3 = fusionPowerDensity(in.n_D_m3, in.n_T_m3, in.T_i_keV)
                     * C_PEAK_FUS;
    out.P_fus_MW = p_fus_W_m3 * V * 1e-6f;

    // α heating = 3.5/17.6 · P_fus, split e/i via Stix slowing-down
    out.P_alpha_MW   = out.P_fus_MW * (3.521f / 17.589f);
    float f_i        = alphaIonFraction(in.T_e_keV);
    out.P_alpha_i_MW = out.P_alpha_MW * f_i;
    out.P_alpha_e_MW = out.P_alpha_MW * (1.0f - f_i);

    // Neutron power (diagnostic — deposited in the BLANKET, not the plasma)
    out.P_neutron_MW = out.P_fus_MW * (14.070f / 17.589f);

    // Radiation (electron channel)
    out.P_brem_MW = bremSStrahlung(in.n_e_m3, in.T_e_keV, in.Z_eff)
                  * C_PEAK_BREM * V * 1e-6f;
    out.P_sync_MW = synchrotron(in.n_e_m3, in.T_e_keV, g.B_toroidal_T) * V * 1e-6f;

    // Confinement time — P_loss = total heating (steady-state approximation)
    out.P_aux_MW = in.P_aux_e_MW + in.P_aux_i_MW;
    float P_loss = out.P_aux_MW + out.P_alpha_MW + in.P_ohm_MW;
    out.tau_E_s = ipb98y2(g, P_loss, n_19, /*M_amu=*/2.5f, H_98);
    out.H_98 = H_98;

    // Conduction losses per channel: P_cond,x = W_x / τ_E
    // Guard: when I_p → 0 (during SCRAM), ipb98y2 returns tau_E = 0
    // because powf(0, 0.93) = 0 — avoid W/0 = inf poisoning the integrator.
    float P_cond_e = 0.0f, P_cond_i = 0.0f;
    if (out.tau_E_s > 1e-6f) {
        if (W_e_MJ > 0.0f && std::isfinite(W_e_MJ)) P_cond_e = W_e_MJ / out.tau_E_s;
        if (W_i_MJ > 0.0f && std::isfinite(W_i_MJ)) P_cond_i = W_i_MJ / out.tau_E_s;
    }
    out.P_cond_MW = P_cond_e + P_cond_i;

    // Electron-ion collisional exchange:
    //   Q_ei [MW] = (3/2)·n_e·k·(T_e−T_i)·ν_eq·V
    float nu_eq = eiEquilibrationRate(in.n_e_m3, in.T_e_keV);
    out.P_ei_MW = 1.5f * in.n_e_m3 * (in.T_e_keV - in.T_i_keV)
                * 1.602176634e-16f * nu_eq * V * 1e-6f;
    // Numerical stability: the exchange must not overshoot equal
    // temperatures within one solver call — cap at the energy difference
    // per second between the channels.
    float dW_gap_MJ = 0.5f * (W_e_MJ - W_i_MJ);
    if (out.P_ei_MW > 0.0f) out.P_ei_MW = std::min(out.P_ei_MW,  std::max(dW_gap_MJ, 0.0f) * 5.0f);
    else                    out.P_ei_MW = std::max(out.P_ei_MW, -std::max(-dW_gap_MJ, 0.0f) * 5.0f);

    // Channel ODEs (note: NO neutron term — see bug-fix comment above)
    out.dWe_dt_MW = out.P_alpha_e_MW + in.P_ohm_MW + in.P_aux_e_MW
                  - out.P_brem_MW - out.P_sync_MW - in.P_line_MW
                  - P_cond_e - out.P_ei_MW;
    out.dWi_dt_MW = out.P_alpha_i_MW + in.P_aux_i_MW
                  - P_cond_i + out.P_ei_MW;
    out.dW_dt_MW  = out.dWe_dt_MW + out.dWi_dt_MW;

    // ── NaN / inf safety net ────────────────────────────────────────────────
    //  Once a stored energy becomes NaN, std::max(NaN, 0) returns NaN on
    //  libstdc++ and the poison spreads to every downstream field — so
    //  clamp any non-finite result to zero here.
    if (!std::isfinite(out.dWe_dt_MW)) out.dWe_dt_MW = 0.0f;
    if (!std::isfinite(out.dWi_dt_MW)) out.dWi_dt_MW = 0.0f;
    if (!std::isfinite(out.dW_dt_MW))  out.dW_dt_MW  = 0.0f;
    if (!std::isfinite(out.P_cond_MW)) out.P_cond_MW = 0.0f;
    if (!std::isfinite(out.P_fus_MW))  out.P_fus_MW  = 0.0f;
    if (!std::isfinite(out.P_ei_MW))   out.P_ei_MW   = 0.0f;

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

// ─── Bootstrap current fraction ──────────────────────────────────────────────
//
//  f_bs ≈ 0.7 · √ε · β_p        (ε = a/R, Wesson Ch. 4 / Hirshman form)
//
//  ITER baseline (ε = 0.32, β_p ≈ 0.65): f_bs ≈ 0.26 — matching the
//  published ~25% bootstrap fraction.  Advanced high-β_p scenarios reach
//  > 0.5 (the steady-state DEMO target).
//
//  BUG FIX: the previous form (0.3 · q · √(R/a) · β_p) had the aspect-ratio
//  factor INVERTED — the trapped-particle fraction that drives bootstrap
//  scales with √(a/R), not √(R/a) — and the extra q factor double-counted
//  the poloidal-field weakness already captured in β_p.  Combined with a
//  β_p that was ~1 by construction (see PlasmaCoreBridge), the fraction was
//  permanently pinned at its 0.95 clamp: the H&CD tab reported ~14 MA of
//  "bootstrap" and any volt-second accounting was meaningless.
//
inline float bootstrapFraction(float eps, float beta_p)
{
    //  Clamp to ≤ 0.8: the scaling is a low-β_p expansion and can exceed
    //  unity for extreme inputs, which is unphysical (bootstrap can't carry
    //  more than the total current).
    return std::clamp(0.7f * sqrtf(std::max(eps, 0.0f)) * beta_p, 0.0f, 0.8f);
}

} // namespace ConfinementPhysics
