#pragma once

//
// src/Helium/HeliumCoolingPhysics.h
// First-wall and divertor thermal physics using helium coolant.
//
//  All functions from helium_cooling.md spec:
//    - Gnielinski Nusselt correlation for turbulent He
//    - Helium thermophysical properties (near-ideal gas)
//    - First-wall surface heat flux
//    - Divertor peak heat load
//    - Thermal stress in tungsten armour
//    - 1D temperature profile through tungsten tile
//

namespace HeliumCoolingPhysics {

// ─── Helium Properties ─────────────────────────────────────────────────────────
//
//  Valid range: 300 K – 1000 K, 1 MPa – 20 MPa.
//  He is near-ideal; deviations from ideal gas <1 % below 100 MPa.
//
struct HeliumProps {
    float cp_J_kgK;    // specific heat at constant pressure [J/(kg·K)]
    float mu_Pa_s;     // dynamic viscosity [Pa·s]
    float k_W_mK;      // thermal conductivity [W/(m·K)]
    float rho_kg_m3;   // density [kg/m³]
    float Pr;          // Prandtl number  μ·cp / k
};

HeliumProps heliumProperties(float T_K, float P_MPa);

// ─── heliumNusselt ────────────────────────────────────────────────────────────
//
//  Gnielinski correlation for developed turbulent pipe flow:
//    Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))
//  Valid: 3000 < Re < 5e6, 0.5 < Pr < 2000.
//
//  For annular or non-circular geometry, Re and Nu are based on hydraulic diameter.
//
enum class GeometryType { Pipe, Annulus, RectangularChannel };

float heliumNusselt(float Re,              // Reynolds number
                    float Pr,              // Prandtl number
                    GeometryType geometry,
                    float L_over_D = 50.0f); // entry correction (L/D ratio)

// ─── firstWallHeatFlux ────────────────────────────────────────────────────────
//
//  Surface heat flux on the first wall from plasma radiation and charge-exchange:
//    q_surface = P_radiation / A_wall   [W/m²]
//
//  ITER first wall: ~0.5 MW/m² (radiation) + ~0.5 MW/m² (particle)
//
struct WallHeatFluxResult {
    float q_total_MW_m2;     // total surface heat flux [MW/m²]
    float q_radiation_MW_m2; // radiation component
    float q_particle_MW_m2;  // charge-exchange neutral particle component
    bool  exceeds_limit;     // true if > 2 MW/m² (ITER design limit)
};

WallHeatFluxResult firstWallHeatFlux(float P_radiation_MW,     // total radiated power [MW]
                                      float P_particle_MW,      // particle heat load [MW]
                                      float A_wall_m2);         // first-wall area [m²]

// ─── divertorHeatLoad ─────────────────────────────────────────────────────────
//
//  Power flowing through the scrape-off layer (SOL) concentrates onto the
//  divertor target with a characteristic flux width λ_q (Eich scaling).
//
//  Peak heat flux:  q_peak = P_SOL / (2π * R * λ_q * f_expansion)
//  Typical λ_q ~ 1–5 mm; f_expansion ~ 4–10 (magnetic flux expansion at target).
//
struct DivertorHeatResult {
    float q_peak_MW_m2;      // peak heat flux on target [MW/m²]
    float q_average_MW_m2;   // average over wetted area
    float wetted_area_m2;    // total wetted area
    bool  tile_damage_risk;  // true if > 10 MW/m² (W recrystallisation)
};

DivertorHeatResult divertorHeatLoad(float P_scrape_off_MW,  // power in SOL [MW]
                                     float R_major_m,        // major radius [m]
                                     float lambda_q_m,       // SOL width [m] (typ: 0.002)
                                     float f_expansion);     // flux tube expansion

// ─── thermalStressInWall ──────────────────────────────────────────────────────
//
//  Thermoelastic stress in a constrained plate under a through-thickness ΔT:
//    σ_thermal = α * E * ΔT / (1 − ν)
//
//  For tungsten armour tile (biaxial stress state):
//    α ≈ 4.5 × 10⁻⁶ K⁻¹ (CTE)
//    E ≈ 411 GPa (Young's modulus)
//    ν ≈ 0.28 (Poisson's ratio)
//
struct ThermalStressResult {
    float sigma_MPa;         // thermoelastic stress [MPa]
    float sigma_yield_MPa;   // material yield strength at T [MPa]
    float safety_factor;     // sigma_yield / sigma (> 1.5 = acceptable)
    bool  yield_risk;        // true if safety_factor < 1.0
    bool  fatigue_risk;      // true if safety_factor < 1.5 (cyclic fatigue)
};

ThermalStressResult thermalStressInWall(float alpha_K,      // CTE [1/K]
                                         float E_Pa,         // Young's modulus [Pa]
                                         float deltaT_K,     // through-thickness ΔT [K]
                                         float nu,           // Poisson's ratio
                                         float T_surface_K); // surface temperature [K]

// ─── tungstenTemperatureProfile ───────────────────────────────────────────────
//
//  1D steady-state heat conduction through a tungsten armour tile:
//    q = -k(T) * dT/dx  →  T(x) = T_back + (q/k) * (thickness - x)
//
//  k(T) for tungsten: k ≈ 174 - 0.082 * T  [W/(m·K)]  (valid 300–2000 K)
//
struct TungstenProfile {
    float T_surface_K;   // front-face (plasma-side) temperature [K]
    float T_midplane_K;  // mid-thickness temperature [K]
    float T_back_K;      // back-face (coolant-side) temperature [K]
    float q_conducted_MW_m2; // conducted power [MW/m²] (should ≈ q_surface)
    bool  above_recrystallisation; // T_surface > 1700 K → W recrystallises
    bool  approaching_melt;        // T_surface > 3400 K → approaching melting point
};

TungstenProfile tungstenTemperatureProfile(float q_surface_MW_m2,  // incident heat flux [MW/m²]
                                            float thickness_m,       // armour thickness [m]
                                            float T_coolant_K);      // coolant temperature [K]

} // namespace HeliumCoolingPhysics