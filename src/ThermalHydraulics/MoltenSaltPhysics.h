#pragma once

//
// src/ThermalHydraulics/MoltenSaltPhysics.h
// Thermal-hydraulics of the molten-salt (FLiBe/FLiNaK) breeding blanket.
//
//  Functions from Molten_Salt_Blanket_Thermal_Hydraulics.md spec:
//    Heat Transfer:
//      heatGenerationInBlanket     — volumetric heating from neutrons [W/m³]
//      dittusBoelterNusselt        — Nu for forced convection
//      moltenSaltViscosity         — FLiBe/FLiNaK μ(T, composition)
//      moltenSaltThermalConductivity — k(T)
//      moltenSaltHeatCapacity      — cp(T, composition)
//      bulkTemperatureRise         — ΔT = q / (ṁ * cp)
//    Flow & Pressure:
//      pressureDrop                — Darcy-Weisbach
//      magnetohydrodynamicDrag     — MHD pressure drop (Hartmann flow)
//      HartmannNumber              — Ha = BL√(σ/μ)
//      naturalConvectionFlow       — passive decay heat removal
//

namespace MoltenSaltPhysics {

// ─── Salt composition types ────────────────────────────────────────────────────
enum class SaltType {
    FLiBe,    // LiF-BeF₂  (2:1 molar): ITER HCLL reference
    FLiNaK,   // LiF-NaF-KF eutectic: higher operating temperature
    FLiBe_Li7 // Li-7 enriched FLiBe (reduced tritium co-generation)
};

// ─── heatGenerationInBlanket ──────────────────────────────────────────────────
//
//  Volumetric nuclear heating in the blanket salt [W/m³]:
//    q_dot = neutron_flux * Σ_eff * E_dep
//  where Σ_eff is the effective macroscopic cross-section and E_dep is the
//  energy deposited per reaction (including gamma heating from (n,γ)).
//
//  For FLiBe: dominant reactions are Li-6(n,t)α  +  Be(n,2n)  +  F(n,γ)
//  Typical volumetric heating: 5–20 MW/m³ near the first wall.
//
struct BlanketHeatGenResult {
    float q_dot_W_m3;            // total volumetric heating [W/m³]
    float q_Li6_reaction_W_m3;   // contribution from Li-6 breeding reaction
    float q_gamma_W_m3;          // gamma heating component
    float total_power_MW;        // integral over blanket volume [MW]
};

BlanketHeatGenResult heatGenerationInBlanket(float neutron_flux_m2s,     // 14 MeV neutron flux [n/m²/s]
                                              float Li6_density_m3,        // Li-6 number density [m^-3]
                                              float blanket_volume_m3,     // total blanket volume [m³]
                                              SaltType salt);

// ─── dittusBoelterNusselt ─────────────────────────────────────────────────────
//
//  Dittus-Boelter correlation for fully developed turbulent pipe flow:
//    Nu = 0.023 * Re^0.8 * Pr^n
//  where n = 0.4 for heating (T_wall > T_fluid), 0.3 for cooling.
//  Valid: Re > 10,000, 0.6 < Pr < 160, L/D > 10.
//
float dittusBoelterNusselt(float Re,      // Reynolds number
                            float Pr,     // Prandtl number
                            bool heating);// true = fluid is being heated

// ─── moltenSaltViscosity ──────────────────────────────────────────────────────
//
//  FLiBe viscosity from Cantor (1973) Arrhenius fit:
//    μ(T) = A * exp(B / T)
//  FLiBe: A = 1.16e-4 Pa·s, B = 3755 K  (valid 733–1073 K)
//  FLiNaK: A = 4.0e-5 Pa·s, B = 4170 K  (valid 730–1570 K)
//
float moltenSaltViscosity(float T_K, SaltType salt);

// ─── moltenSaltThermalConductivity ────────────────────────────────────────────
//
//  FLiBe k(T) — polynomial fit to Williams (2006) data:
//    k = 0.629697 + 5.0e-4 * T   [W/(m·K)]  valid 733–1073 K
//  FLiNaK: approximately constant k ≈ 0.6 W/(m·K) over operating range.
//
float moltenSaltThermalConductivity(float T_K, SaltType salt);

// ─── moltenSaltHeatCapacity ───────────────────────────────────────────────────
//
//  FLiBe cp (Sohal 2010):
//    cp = 2415.78 J/(kg·K)  (approximately constant 700–1000 K, within 2 %)
//  FLiNaK: cp ≈ 1905 J/(kg·K)
//
float moltenSaltHeatCapacity(float T_K, SaltType salt);

// ─── bulkTemperatureRise ──────────────────────────────────────────────────────
//
//  ΔT_bulk = Q / (ṁ * cp)
//  where Q = total heat added [W], ṁ = mass flow rate [kg/s].
//  Returns temperature rise [K].
//
float bulkTemperatureRise(float q_dot_W,      // total heat input [W]
                           float mdot_kg_s,    // mass flow rate [kg/s]
                           float cp_J_kgK);    // specific heat [J/(kg·K)]

// ─── pressureDrop ─────────────────────────────────────────────────────────────
//
//  Darcy-Weisbach pressure drop along a pipe:
//    ΔP = f_D * (L/D) * (ρ * v² / 2)
//  where f_D = Darcy friction factor (Colebrook-White for turbulent flow).
//
struct PressureDropResult {
    float dP_Pa;             // total pressure drop [Pa]
    float f_darcy;           // Darcy-Weisbach friction factor
    float v_m_s;             // mean velocity [m/s]
    float Re;                // Reynolds number
    float pump_power_W;      // required pumping power = Q_vol * ΔP [W]
};

PressureDropResult pressureDrop(float L_m,          // pipe length [m]
                                 float D_m,          // hydraulic diameter [m]
                                 float rho_kg_m3,    // fluid density [kg/m³]
                                 float v_m_s,        // mean velocity [m/s]
                                 float f_darcy,      // Darcy friction factor (pass 0 to compute)
                                 float mu_Pa_s);     // dynamic viscosity [Pa·s]

// ─── magnetohydrodynamicDrag ──────────────────────────────────────────────────
//
//  MHD pressure drop for conducting fluid in a transverse magnetic field:
//    ΔP_MHD = σ * v * B² * L  (for Ha >> 1, Hartmann flow in rectangular duct)
//
//  The Hartmann number Ha determines the flow regime:
//    Ha < 1:  negligible MHD effect
//    Ha ~ 10: transition regime
//    Ha > 100: fully MHD-dominated, flat core profile + thin Hartmann layers
//
struct MHDResult {
    float dP_MHD_Pa;         // additional MHD pressure drop [Pa]
    float Ha;                // Hartmann number
    float velocity_profile;  // 1.0 = flat (MHD-dominated), 0.0 = parabolic (no MHD)
    float effective_k_W_mK;  // effective k (MHD suppresses turbulence → lower k)
};

MHDResult magnetohydrodynamicDrag(float sigma_S_m,    // electrical conductivity [S/m]
                                   float B_T,          // transverse magnetic field [T]
                                   float v_m_s,        // mean flow velocity [m/s]
                                   float L_char_m,     // characteristic length (duct half-width) [m]
                                   float mu_Pa_s,      // dynamic viscosity [Pa·s]
                                   float rho_kg_m3);   // density [kg/m³]

// ─── HartmannNumber ───────────────────────────────────────────────────────────
//
//  Ha = B * L * sqrt(σ / μ)
//  Critical value:  Ha > 1 → MHD effects significant.
//
float HartmannNumber(float B_T,            // magnetic field [T]
                     float L_m,            // characteristic length [m]
                     float sigma_S_m,      // electrical conductivity [S/m]
                     float mu_Pa_s);       // dynamic viscosity [Pa·s]

// ─── naturalConvectionFlow ────────────────────────────────────────────────────
//
//  Passive decay heat removal via buoyancy-driven natural convection.
//  Rayleigh-number correlation for vertical cavity:
//    Ra = g * β * ΔT * L³ / (ν * α)
//  Nu = 0.15 * Ra^(1/3)  (turbulent natural convection, Ra > 10^9)
//
//  Returns volumetric flow rate [m³/s] driven by natural convection.
//
struct NatConvResult {
    float Ra;                // Rayleigh number
    float Nu;                // Nusselt number
    float Q_removed_W;       // heat removed by natural convection [W]
    float flow_velocity_ms;  // approximate buoyancy velocity [m/s]
    bool  adequate_cooling;  // true if Q_removed ≥ decay_power target
};

NatConvResult naturalConvectionFlow(float rho_kg_m3,       // salt density [kg/m³]
                                     float beta_K,          // volumetric expansion coeff [1/K]
                                     float deltaT_K,        // driving temperature diff [K]
                                     float L_m,             // vertical height of loop [m]
                                     float k_W_mK,          // thermal conductivity [W/(m·K)]
                                     float cp_J_kgK,        // specific heat [J/(kg·K)]
                                     float mu_Pa_s,         // viscosity [Pa·s]
                                     float decay_power_W);  // target removal power [W]

} // namespace MoltenSaltPhysics