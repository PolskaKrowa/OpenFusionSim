//
// src/Magnets/MagnetPhysics.cpp
//

#include "MagnetPhysics.h"
#include <cmath>
#include <algorithm>

namespace MagnetPhysics {

static constexpr float PI  = 3.14159265358979f;
static constexpr float MU0 = 1.25663706e-6f;

// ─── toroidalFieldOnAxis ───────────────────────────────────────────────────────
float toroidalFieldOnAxis(int N_turns, float I_A, float R_major_m)
{
    // B_axis = μ₀ * N * I / (2π * R)
    if (R_major_m < 1e-6f) return 0.0f;
    return MU0 * static_cast<float>(N_turns) * I_A / (2.0f * PI * R_major_m);
}

// ─── poloidalFieldCoilForces ──────────────────────────────────────────────────
//
//  Mutual inductance between two coaxial circular loops (Neumann formula):
//    M = μ₀ * sqrt(R1*R2) * ((2/k - k) * K(k) - (2/k) * E(k))
//  where k² = 4R1R2 / ((R1+R2)² + d²)
//  K, E = complete elliptic integrals of the first and second kind.
//
//  We use the AGM approximation for K and E accurate to ~1e-6.
//
static float ellipticK(float k2)
{
    // Arithmetic-Geometric Mean method for K(k)
    float kp2 = 1.0f - k2;
    float a   = 1.0f, b = sqrtf(kp2), c;
    for (int i = 0; i < 10; i++) {
        c = 0.5f * (a + b);
        b = sqrtf(a * b);
        a = c;
    }
    return PI / (2.0f * a);
}

static float ellipticE(float k2)
{
    // Via the relation E(k) = (1 - k²/2) * K(k) -- Gauss AGM not trivial for E,
    // use the Carlson series approximation instead:
    float a0 = 1.0f, g0 = sqrtf(1.0f - k2), s = 0.5f;
    float sum = 1.0f - k2 * 0.5f;
    float power = 0.25f;
    float a = a0, g = g0;
    for (int i = 0; i < 10; i++) {
        float a_new = 0.5f * (a + g);
        float g_new = sqrtf(a * g);
        float c2    = (a - g) * (a - g);
        sum   -= power * c2;
        power *= 2.0f;
        a = a_new; g = g_new;
    }
    return PI * sum / (2.0f * a);
}

CoilForceResult poloidalFieldCoilForces(float I1_A, float I2_A,
                                         float R1_m, float R2_m,
                                         float d_m)
{
    CoilForceResult res{};

    float denom = (R1_m + R2_m) * (R1_m + R2_m) + d_m * d_m;
    float k2    = 4.0f * R1_m * R2_m / denom;
    k2          = std::clamp(k2, 0.0f, 0.9999f);

    float K = ellipticK(k2);
    float E = ellipticE(k2);
    float k  = sqrtf(k2);

    // Neumann mutual inductance
    float M = MU0 * sqrtf(R1_m * R2_m)
            * ((2.0f / k - k) * K - (2.0f / k) * E);
    res.mutual_H = M;

    // Axial force: F_z = I1 * I2 * dM/dd
    // dM/dd analytically involves derivatives of K,E with respect to d.
    // Use numerical difference: δd = d * 1e-4
    float dd = std::max(d_m * 1e-4f, 1e-6f);
    float denom2 = (R1_m + R2_m)*(R1_m + R2_m) + (d_m + dd)*(d_m + dd);
    float k2b    = 4.0f * R1_m * R2_m / denom2;
    k2b          = std::clamp(k2b, 0.0f, 0.9999f);
    float kb     = sqrtf(k2b);
    float M2 = MU0 * sqrtf(R1_m * R2_m)
             * ((2.0f / kb - kb) * ellipticK(k2b) - (2.0f / kb) * ellipticE(k2b));

    float dM_dd = (M2 - M) / dd;
    res.F_axial_N  = I1_A * I2_A * dM_dd;   // [N], negative = attractive

    // Radial (centering) force from hoop stress — approximate as:
    // F_radial ≈ I * B_at_coil * L_coil  — not derived from mutual inductance
    // Here: F_r ≈ μ₀ * I1² * R1 * ln(8R1/a_wire - 2)  (self-force, hoop tension)
    // Simplified: use the interaction field from coil 2 at the position of coil 1
    float B2_at_R1 = MU0 * I2_A * R2_m * R2_m
                   / (2.0f * powf(R2_m * R2_m + d_m * d_m, 1.5f)); // on-axis approx
    float circumference_1 = 2.0f * PI * R1_m;
    res.F_radial_N = I1_A * B2_at_R1 * circumference_1;  // Lorentz [N]

    return res;
}

// ─── quenchPropagation ────────────────────────────────────────────────────────

// Material property table at 4.5 K, 12 T (nominal operating point)
struct SCProperties {
    float T_critical_K;     // critical temperature at zero field [K]
    float B_critical_T;     // upper critical field at 0 K [T]
    float rho_n_ohm_m;      // normal-state resistivity (stabiliser) [Ω·m]
    float Cv_J_m3K;         // volumetric heat capacity [J/(m³·K)]
    float k_thermal_W_mK;   // thermal conductivity in normal state [W/(m·K)]
    float J_c0_Am2;         // critical current density at 4.5 K, 12 T [A/m²]
};

static SCProperties getMaterialProps(SCMaterial mat, float B_T)
{
    (void)B_T;
    switch (mat) {
        case SCMaterial::Nb3Sn:
            return { 18.0f, 26.0f, 1.5e-9f, 3000.0f, 200.0f, 700e6f };
        case SCMaterial::NbTi:
            return {  9.6f, 14.5f, 1.2e-9f, 2500.0f, 150.0f, 500e6f };
        case SCMaterial::REBCO:
            return { 92.0f, 100.f, 0.8e-9f, 2000.0f, 300.0f, 3000e6f };
        default:
            return { 18.0f, 26.0f, 1.5e-9f, 3000.0f, 200.0f, 700e6f };
    }
}

QuenchResult quenchPropagation(float T_K,
                               float B_T,
                               float J_Am2,
                               SCMaterial material,
                               float strand_length_m,
                               float detect_threshold_V)
{
    QuenchResult res{};
    auto p = getMaterialProps(material, B_T);

    // Current-sharing temperature: T_cs = T_op + (T_c - T_op) * (1 - J/J_c)
    // Simplified linear model
    float J_c = p.J_c0_Am2 * (1.0f - B_T / p.B_critical_T)
                           * (1.0f - T_K   / p.T_critical_K);
    J_c = std::max(J_c, 0.0f);

    float T_cs = T_K + (p.T_critical_K - T_K) * (1.0f - J_Am2 / (J_c + 1.0f));
    T_cs       = std::clamp(T_cs, T_K, p.T_critical_K);

    res.quench_confirmed = (T_K >= T_cs) || (J_Am2 > J_c);

    if (!res.quench_confirmed) {
        res.v_normal_zone_ms  = 0.0f;
        res.T_hotspot_K       = T_K;
        res.t_detect_s        = 1e9f; // not quenching
        res.energy_deposited_J= 0.0f;
        return res;
    }

    // Normal-zone propagation velocity (adiabatic model):
    //   v = J * sqrt(ρ_n * k) / (C_v * (T_cs - T_op))
    float delta_T = T_cs - T_K + 0.1f; // avoid /0
    float v = J_Am2 * sqrtf(p.rho_n_ohm_m * p.k_thermal_W_mK)
            / (p.Cv_J_m3K * delta_T);
    res.v_normal_zone_ms = v;

    // Detection time: voltage = J * ρ_n * v * t_detect ≥ threshold
    // V = ρ_n * J * (zone length) = ρ_n * J * v * t
    if (J_Am2 * p.rho_n_ohm_m * v > 1e-20f) {
        res.t_detect_s = detect_threshold_V
                       / (J_Am2 * p.rho_n_ohm_m * v);
    } else {
        res.t_detect_s = 1.0f;
    }

    // Hot-spot temperature (adiabatic): MIIT integral
    // T_hs ~ T_cs + (J² * ρ_n * t_detect) / C_v
    res.T_hotspot_K = T_cs + (J_Am2 * J_Am2 * p.rho_n_ohm_m * res.t_detect_s)
                           / p.Cv_J_m3K;

    // Energy deposited in normal zone during detection window
    // E = J² * ρ_n * V_wire * t_detect  (Joule heating)
    float A_strand = J_Am2 > 0 ? 1e-6f : 0.0f; // assume 1 mm² strand cross-section
    float V_wire   = A_strand * v * res.t_detect_s;
    res.energy_deposited_J = J_Am2 * J_Am2 * p.rho_n_ohm_m * V_wire * res.t_detect_s;

    return res;
}

// ─── cryostatHeatLeak ─────────────────────────────────────────────────────────
CryostatHeatLeak cryostatHeatLeak(float A_m2,
                                   float k_insulation_Wm_K,
                                   float thickness_m,
                                   float T_warm_K,
                                   float T_cold_K,
                                   int   n_MLI_layers)
{
    CryostatHeatLeak res{};

    // Conduction through insulation bulk: Q = k * A * ΔT / L
    if (thickness_m > 1e-6f)
        res.Q_conduction_W = k_insulation_Wm_K * A_m2
                           * (T_warm_K - T_cold_K) / thickness_m;

    // MLI radiation suppression: each layer reduces radiation by a factor
    // Effective emissivity: ε_eff ≈ ε_single / (N+1)   (N = number of layers)
    // Radiation without MLI: Q_rad = σ * ε * A * (T_warm^4 - T_cold^4)
    constexpr float SIGMA = 5.67e-8f; // Stefan-Boltzmann [W/(m²·K⁴)]
    constexpr float EPS   = 0.05f;    // emissivity of polished Al MLI foil
    float eps_eff = (n_MLI_layers > 0)
                  ? EPS / (n_MLI_layers + 1)
                  : EPS * 0.9f; // bare surface
    float T4_diff = (T_warm_K * T_warm_K * T_warm_K * T_warm_K)
                  - (T_cold_K * T_cold_K * T_cold_K * T_cold_K);
    res.Q_radiation_W = SIGMA * eps_eff * A_m2 * T4_diff;

    res.Q_total_W = res.Q_conduction_W + res.Q_radiation_W;
    return res;
}

// ─── inductiveStoredEnergy ────────────────────────────────────────────────────
float inductiveStoredEnergy(float L_H, float I_A)
{
    return 0.5f * L_H * I_A * I_A;
}

float plasmaInductance(float R_major_m, float a_minor_m, float li, float kappa)
{
    // Shafranov formula:
    // L_p = μ₀ * R * (ln(8R / a_eff) - 2 + li/2)
    // where a_eff = a * sqrt(kappa)  (effective minor radius for elongated plasma)
    float a_eff = a_minor_m * sqrtf(kappa);
    if (a_eff < 1e-6f) return 0.0f;
    float arg = 8.0f * R_major_m / a_eff;
    return MU0 * R_major_m * (logf(arg) - 2.0f + li * 0.5f);
}

} // namespace MagnetPhysics