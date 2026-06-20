//
// src/Helium/HeliumCoolingPhysics.cpp
//

#include "HeliumCoolingPhysics.h"
#include <cmath>
#include <algorithm>

namespace HeliumCoolingPhysics {

// ─── heliumProperties ─────────────────────────────────────────────────────────
//
//  Property fits from NIST/REFPROP data for He-4 (near-ideal gas).
//  cp is nearly constant over 300–1000 K for He (monoatomic: cp = 5R/2M).
//
HeliumProps heliumProperties(float T_K, float P_MPa)
{
    HeliumProps p{};

    // He-4: M = 4.0026 g/mol, monoatomic ideal gas
    constexpr float R_specific = 2077.0f;   // specific gas constant [J/(kg·K)]
    constexpr float cp         = 5193.0f;   // cp = (5/2) * R_specific [J/(kg·K)] (constant)

    // Ideal gas density: ρ = P / (R * T)
    p.rho_kg_m3 = (P_MPa * 1e6f) / (R_specific * T_K);
    p.cp_J_kgK  = cp;

    // Viscosity: power-law fit μ(T) ≈ μ_ref * (T/T_ref)^0.67  (Sutherland-like)
    // μ_ref = 19.9e-6 Pa·s at T_ref = 300 K
    p.mu_Pa_s = 19.9e-6f * powf(T_K / 300.0f, 0.67f);

    // Thermal conductivity: k(T) ≈ k_ref * (T/T_ref)^0.67
    // k_ref = 0.152 W/(m·K) at 300 K
    p.k_W_mK = 0.152f * powf(T_K / 300.0f, 0.67f);

    // Prandtl: Pr = μ * cp / k  (should be ~0.67 for monoatomic gas)
    p.Pr = p.mu_Pa_s * p.cp_J_kgK / p.k_W_mK;

    return p;
}

// ─── heliumNusselt ────────────────────────────────────────────────────────────
float heliumNusselt(float Re, float Pr, GeometryType geometry, float L_over_D)
{
    // Gnielinski correlation (valid Re > 3000, Pr 0.5–2000)
    Re = std::max(Re, 2300.0f); // Gnielinski requires turbulent regime

    // Darcy friction factor from Petukhov equation:
    //   f = (0.790 * ln(Re) - 1.64)^-2
    float lnRe = logf(Re);
    float f = 1.0f / ((0.790f * lnRe - 1.64f) * (0.790f * lnRe - 1.64f));

    float f_8     = f / 8.0f;
    float Pr_term = powf(Pr, 2.0f / 3.0f) - 1.0f;
    float Nu = f_8 * (Re - 1000.0f) * Pr
             / (1.0f + 12.7f * sqrtf(f_8) * Pr_term);

    // Thermal entry-length correction (Hausen factor):
    //   Nu_corrected = Nu * (1 + (D/L)^0.7)  (approximate)
    if (L_over_D > 0.0f) {
        Nu *= (1.0f + powf(1.0f / L_over_D, 0.7f));
    }

    // Geometry correction for non-circular channels
    switch (geometry) {
        case GeometryType::Annulus:
            Nu *= 0.86f;  // annular channels: ~14 % lower than pipe
            break;
        case GeometryType::RectangularChannel:
            Nu *= 0.93f;  // aspect-ratio-dependent; 0.93 for AR~2
            break;
        case GeometryType::Pipe:
        default:
            break;
    }

    return Nu;
}

// ─── firstWallHeatFlux ────────────────────────────────────────────────────────
WallHeatFluxResult firstWallHeatFlux(float P_radiation_MW,
                                      float P_particle_MW,
                                      float A_wall_m2)
{
    WallHeatFluxResult res{};
    if (A_wall_m2 < 0.1f) return res;

    res.q_radiation_MW_m2 = P_radiation_MW / A_wall_m2;
    res.q_particle_MW_m2  = P_particle_MW  / A_wall_m2;
    res.q_total_MW_m2     = res.q_radiation_MW_m2 + res.q_particle_MW_m2;

    // ITER design limit: 2 MW/m² (normal operation), 20 MW/m² (ELMs, transients)
    res.exceeds_limit = (res.q_total_MW_m2 > 2.0f);
    return res;
}

// ─── divertorHeatLoad ─────────────────────────────────────────────────────────
DivertorHeatResult divertorHeatLoad(float P_scrape_off_MW,
                                     float R_major_m,
                                     float lambda_q_m,
                                     float f_expansion)
{
    DivertorHeatResult res{};
    if (lambda_q_m < 1e-5f || f_expansion < 0.1f) return res;

    constexpr float PI = 3.14159265f;

    // Wetted area: annular strip of width λ_q * f_expansion at divertor target
    float lambda_target = lambda_q_m * f_expansion;
    res.wetted_area_m2  = 2.0f * PI * R_major_m * lambda_target;  // single-null

    // Peak heat flux assumes an exponential decay profile across λ_q:
    // q_peak = P_SOL / (2π * R * λ_q * f_expansion * S_factor)
    // S_factor accounts for tilted target (typically 1.5–2 for ITER)
    constexpr float S_tilt = 1.5f;
    if (res.wetted_area_m2 > 0.0f) {
        res.q_peak_MW_m2    = P_scrape_off_MW / (res.wetted_area_m2 * S_tilt);
        res.q_average_MW_m2 = P_scrape_off_MW /  res.wetted_area_m2;
    }

    // W damage threshold: recrystallisation starts above 10 MW/m² (steady)
    res.tile_damage_risk = (res.q_peak_MW_m2 > 10.0f);
    return res;
}

// ─── thermalStressInWall ──────────────────────────────────────────────────────
ThermalStressResult thermalStressInWall(float alpha_K,
                                         float E_Pa,
                                         float deltaT_K,
                                         float nu,
                                         float T_surface_K)
{
    ThermalStressResult res{};

    // Biaxial thermoelastic stress (plane-stress, constrained edges):
    //   σ = α * E * ΔT / (1 − ν)
    res.sigma_MPa = (alpha_K * E_Pa * deltaT_K / (1.0f - nu)) * 1e-6f;

    // Tungsten yield strength (temperature-dependent):
    // σ_y(T) ≈ 550 MPa at 300 K, drops to ~100 MPa at 1800 K
    // Linear fit: σ_y = 550 - 0.25 * (T - 300)  [MPa]
    res.sigma_yield_MPa = std::max(550.0f - 0.25f * (T_surface_K - 300.0f), 50.0f);

    res.safety_factor = (res.sigma_MPa > 1e-6f)
                      ? res.sigma_yield_MPa / res.sigma_MPa
                      : 1e9f;
    res.yield_risk   = (res.safety_factor < 1.0f);
    res.fatigue_risk = (res.safety_factor < 1.5f);
    return res;
}

// ─── tungstenTemperatureProfile ───────────────────────────────────────────────
TungstenProfile tungstenTemperatureProfile(float q_surface_MW_m2,
                                            float thickness_m,
                                            float T_coolant_K)
{
    TungstenProfile prof{};

    // W thermal conductivity: k(T) ≈ 174 - 0.082*T  [W/(m·K)]
    // (valid ~300–2000 K; use mean temperature iteration)
    float q_W_m2 = q_surface_MW_m2 * 1e6f;

    // Iterative solution for T_surface (since k depends on T):
    // Start with a guess, then refine
    float T_back  = T_coolant_K;
    float T_surf  = T_back;
    for (int iter = 0; iter < 10; iter++) {
        float T_mean = 0.5f * (T_surf + T_back);
        float k_mean = std::max(174.0f - 0.082f * T_mean, 60.0f); // clamp to 60 W/mK minimum
        T_surf = T_back + q_W_m2 * thickness_m / k_mean;
    }

    prof.T_back_K      = T_back;
    prof.T_surface_K   = T_surf;
    prof.T_midplane_K  = 0.5f * (T_surf + T_back);

    float k_mid        = std::max(174.0f - 0.082f * prof.T_midplane_K, 60.0f);
    float T_diff       = prof.T_surface_K - prof.T_back_K;
    prof.q_conducted_MW_m2 = (thickness_m > 1e-6f)
                           ? k_mid * T_diff / thickness_m * 1e-6f : 0.0f;

    // ITER W armour limits
    prof.above_recrystallisation = (prof.T_surface_K > 1700.0f);
    prof.approaching_melt        = (prof.T_surface_K > 3400.0f);

    return prof;
}

} // namespace HeliumCoolingPhysics