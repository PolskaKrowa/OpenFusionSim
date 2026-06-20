#pragma once

//
// src/SteamPower/SteamCyclePhysics.h
// Rankine cycle and heat-exchanger physics.
//
//  Functions from steam_power_cycle.md spec:
//    heatExchangerEffectiveness  — ε-NTU method
//    steamProperties             — IAPWS-97 approximation (h, s, ρ, cp)
//    rankineEfficiency           — Carnot-bounded actual cycle efficiency
//    turbineWork                 — W = ṁ · η · (h_in − h_out_s)
//    condenserHeatRejection      — Q_rejected to cold sink
//    feedwaterPumpWork           — parasitic pumping power
//    steamGeneratorDesign        — size HX from duty / LMTD / U
//

namespace SteamCyclePhysics {

// ─── steamProperties ─────────────────────────────────────────────────────────
//
//  IAPWS-97 region approximations (Region 1: liquid, Region 2: vapour).
//  Valid: 273–1073 K, 0.001–100 MPa.
//  Full IAPWS-97 requires ~50 coefficients per region; this is a 4-term
//  Helmholtz energy approximation accurate to ~0.5 % for the steam cycle range.
//
struct SteamState {
    float h_J_kg;       // specific enthalpy [J/kg]
    float s_J_kgK;      // specific entropy [J/(kg·K)]
    float rho_kg_m3;    // density [kg/m³]
    float cp_J_kgK;     // specific heat at constant pressure [J/(kg·K)]
    float mu_Pa_s;      // dynamic viscosity [Pa·s]
    float k_W_mK;       // thermal conductivity [W/(m·K)]
    float Pr;           // Prandtl number
    bool  superheated;  // true if T > T_sat(P)
    float x_quality;    // steam quality (0=sat liquid, 1=sat vapour, >1=superheated)
};

SteamState steamProperties(float T_K, float P_MPa);

// Saturation temperature as a function of pressure (Antoine-type fit to IAPWS)
float saturationTemperature(float P_MPa);
float saturationPressure(float T_K);

// ─── heatExchangerEffectiveness ───────────────────────────────────────────────
//
//  ε-NTU method for heat exchanger sizing and performance.
//  Returns effectiveness ε (0–1) given NTU and heat capacity rate ratio.
//
//  ε = Q_actual / Q_max   where   Q_max = C_min * (T_hot_in - T_cold_in)
//
enum class HXConfig { CounterFlow, ParallelFlow, CrossFlow_UnmixedBoth, Shell1Pass2 };

struct HXResult {
    float epsilon;          // effectiveness [0–1]
    float Q_W;              // actual heat transferred [W]
    float T_hot_out_K;      // hot-side outlet temperature [K]
    float T_cold_out_K;     // cold-side outlet temperature [K]
    float NTU;              // number of transfer units
};

HXResult heatExchangerEffectiveness(float NTU,
                                     float C_hot_W_K,    // hot-side heat capacity rate ṁ·cp [W/K]
                                     float C_cold_W_K,   // cold-side heat capacity rate [W/K]
                                     float T_hot_in_K,
                                     float T_cold_in_K,
                                     HXConfig config);

// ─── rankineEfficiency ────────────────────────────────────────────────────────
//
//  Actual Rankine cycle thermal efficiency, bounded by Carnot:
//    η_Carnot = 1 - T_cold / T_hot
//    η_Rankine = η_Carnot * η_internal  (accounts for irreversibilities)
//
//  η_internal ≈ 0.85–0.92 for modern supercritical steam plants.
//  Also accounts for feedwater pump parasitic.
//
struct RankineResult {
    float eta_Carnot;        // Carnot limit
    float eta_Rankine;       // actual cycle efficiency
    float eta_internal;      // internal efficiency (turbine + pump losses)
    float specific_work_J_kg;// net work per kg of working fluid [J/kg]
};

RankineResult rankineEfficiency(float T_hot_K,          // turbine inlet temperature [K]
                                 float T_cold_K,         // condenser temperature [K]
                                 float eta_turbine,      // isentropic turbine efficiency
                                 float eta_pump);        // isentropic pump efficiency

// ─── turbineWork ──────────────────────────────────────────────────────────────
//
//  Actual turbine shaft work for an isentropic efficiency η:
//    W_actual = ṁ · η · (h_in − h_out_isentropic)
//
//  h_out_isentropic is determined from isentropic expansion to P_out
//  (requires knowing s_in and finding T at s_in, P_out — computed here).
//
struct TurbineResult {
    float W_shaft_W;         // shaft power output [W]
    float h_out_actual_J_kg; // actual outlet enthalpy [J/kg]
    float T_out_K;           // outlet temperature [K]
    float eta_used;          // isentropic efficiency used
};

TurbineResult turbineWork(float mdot_kg_s,    // mass flow rate [kg/s]
                           float h_in_J_kg,   // inlet enthalpy [J/kg]
                           float s_in_J_kgK,  // inlet entropy [J/(kg·K)]
                           float P_out_MPa,   // outlet pressure [MPa]
                           float eta_isentropic);

// ─── condenserHeatRejection ───────────────────────────────────────────────────
//
//  Heat rejected to the cold sink (cooling tower or sea):
//    Q_rejected = ṁ · (h_in − h_out)
//  where h_out = saturated liquid enthalpy at condenser pressure.
//
struct CondenserResult {
    float Q_rejected_W;      // heat rejection rate [W]
    float T_sat_K;           // saturation temperature at P_cond [K]
    float mdot_coolant_kg_s; // required coolant flow [kg/s]
    float LMTD_K;            // log-mean temperature difference [K]
};

CondenserResult condenserHeatRejection(float mdot_steam_kg_s,  // steam flow [kg/s]
                                        float h_in_J_kg,        // turbine exhaust enthalpy [J/kg]
                                        float P_cond_MPa,       // condenser pressure [MPa]
                                        float T_coolant_in_K,   // cooling water inlet [K]
                                        float T_coolant_out_K); // cooling water outlet [K]

// ─── feedwaterPumpWork ────────────────────────────────────────────────────────
//
//  Isentropic pump work (nearly incompressible liquid):
//    W_pump = ṁ · ΔP / (ρ · η_pump)
//
struct PumpResult {
    float W_pump_W;         // pump shaft power [W]
    float h_out_J_kg;       // feedwater outlet enthalpy [J/kg]
    float eta_used;
};

PumpResult feedwaterPumpWork(float mdot_kg_s,   // mass flow rate [kg/s]
                              float dP_Pa,       // pressure rise [Pa]
                              float rho_kg_m3,   // liquid density [kg/m³]
                              float h_in_J_kg,   // inlet enthalpy [J/kg]
                              float eta_pump);   // isentropic pump efficiency

// ─── steamGeneratorDesign ─────────────────────────────────────────────────────
//
//  Size the primary-to-secondary heat exchanger:
//    Q = U * A * LMTD   →   A = Q / (U * LMTD)
//
//  LMTD = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)   (counter-flow)
//
struct SteamGenDesign {
    float A_m2;              // required heat transfer area [m²]
    float LMTD_K;            // log-mean temperature difference [K]
    float NTU;               // number of transfer units
    float U_W_m2K;           // overall heat transfer coefficient used
    float Q_duty_W;          // thermal duty [W]
};

SteamGenDesign steamGeneratorDesign(float Q_duty_W,      // required thermal duty [W]
                                     float LMTD_K,        // log-mean temperature difference [K]
                                     float U_W_m2K,       // overall U [W/(m²·K)]
                                     float C_hot_W_K,     // hot-side heat capacity rate [W/K]
                                     float C_cold_W_K);   // cold-side heat capacity rate [W/K]

} // namespace SteamCyclePhysics