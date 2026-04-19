//
// src/SteamPower/SteamPower.cpp
// Rankine cycle power conversion system.
// Reads blanket thermal output from ReactorState, writes electrical outputs.
//

#include "SteamPower.h"
#include "SteamCyclePhysics.h"
#include <algorithm>
#include <cmath>

SteamPowerSystem::SteamPowerSystem(const SteamConfig& cfg)
    : cfg_(cfg)
{}

void SteamPowerSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    updateSteamGenerator(state, dt);
    updateTurbine        (state, dt);
    updateGenerator      (state, dt);

    // Write outputs
    state.gross_electric_MW  = gross_electric_MW_;
    state.parasitic_load_MW  = cfg_.pump_load_MW + cfg_.magnet_load_MW + cfg_.aux_load_MW;
    state.net_electric_MW    = net_electric_MW_;
    state.turbine_rpm        = turbine_rpm_;
    state.steam_pressure_MPa = steam_pressure_MPa_;
    state.steam_temp_K       = steam_temp_K_;
    state.turbine_trip       = (turbine_rpm_ > cfg_.nominal_rpm * 1.1f)   // overspeed
                             || (steam_pressure_MPa_ > cfg_.nominal_pressure_MPa * 1.15f); // overpressure

    // Scientific Q factor
    float P_aux_MW = 50.0f; // NBI + ICRH auxiliary heating (fixed assumption)
    state.Q_scientific = (P_aux_MW > 0.1f)
                       ? state.fusion_power_MW / P_aux_MW : 0.0f;
}

void SteamPowerSystem::updateSteamGenerator(ReactorState& state, float dt)
{
    // Thermal input = blanket heat + divertor heat (via separate cooling loop)
    thermal_input_MW_ = state.blanket_heat_MW;

    if (thermal_input_MW_ < 0.1f) {
        steam_pressure_MPa_ = std::max(steam_pressure_MPa_ - 0.5f * dt, 0.1f);
        steam_temp_K_       = std::max(steam_temp_K_       - 50.0f * dt, 300.0f);
        steam_flow_kg_s_    = std::max(steam_flow_kg_s_    - 100.0f * dt, 0.0f);
        return;
    }

    // Use ε-NTU to find actual steam conditions from primary coolant
    float C_primary  = state.coolant_flow_kg_s * 5200.0f; // FLiBe: cp ≈ 5200 J/(kg·K) (placeholder — real value from MoltenSaltPhysics)
    float C_secondary= cfg_.nominal_flow_kg_s  * 2100.0f; // steam: cp ≈ 2100 J/(kg·K) at ~600°C

    // HX area → NTU = UA / C_min
    float UA   = 150e6f;  // overall UA [W/K]  (designed for nominal 1500 MW duty)
    float C_min= std::min(C_primary, C_secondary);
    float NTU  = (C_min > 1.0f) ? UA / C_min : 0.0f;

    auto hx = SteamCyclePhysics::heatExchangerEffectiveness(
        NTU, C_primary, C_secondary,
        state.coolant_outlet_temp_K,   // primary hot side in
        500.0f,                        // secondary feedwater in [K]
        SteamCyclePhysics::HXConfig::CounterFlow);

    // Steam generator outlet conditions
    steam_temp_K_       = hx.T_cold_out_K;
    steam_pressure_MPa_ = std::clamp(
        cfg_.nominal_pressure_MPa * (thermal_input_MW_ / 1500.0f), // scale with heat input
        0.1f, cfg_.nominal_pressure_MPa * 1.05f);

    // Mass flow driven by pressure differential (turbine demand)
    steam_flow_kg_s_ = std::min(
        (thermal_input_MW_ * 1e6f) /
        std::max(SteamCyclePhysics::steamProperties(steam_temp_K_, steam_pressure_MPa_).h_J_kg
               - SteamCyclePhysics::steamProperties(500.0f, 0.02f).h_J_kg, 1.0f),
        cfg_.nominal_flow_kg_s);
    steam_flow_kg_s_ = std::max(steam_flow_kg_s_, 0.0f);

    (void)dt;
}

void SteamPowerSystem::updateTurbine(ReactorState& state, float dt)
{
    if (steam_flow_kg_s_ < 1.0f || steam_pressure_MPa_ < 0.2f) {
        turbine_rpm_ = std::max(turbine_rpm_ - 50.0f * dt, 0.0f);
        return;
    }

    auto steam_in = SteamCyclePhysics::steamProperties(steam_temp_K_, steam_pressure_MPa_);

    // Turbine expansion to condenser pressure (~3.5 kPa, typical LP condenser)
    constexpr float P_cond_MPa = 0.0035f;
    auto turb = SteamCyclePhysics::turbineWork(
        steam_flow_kg_s_,
        steam_in.h_J_kg,
        steam_in.s_J_kgK,
        P_cond_MPa,
        cfg_.turbine_efficiency);

    gross_electric_MW_ = turb.W_shaft_W * cfg_.generator_efficiency * 1e-6f;

    // Feedwater pump parasitic (pump condensate back to boiler pressure)
    float rho_fw = SteamCyclePhysics::steamProperties(
        SteamCyclePhysics::saturationTemperature(P_cond_MPa), P_cond_MPa).rho_kg_m3;
    float h_fw   = SteamCyclePhysics::steamProperties(
        SteamCyclePhysics::saturationTemperature(P_cond_MPa), P_cond_MPa).h_J_kg;

    auto pump = SteamCyclePhysics::feedwaterPumpWork(
        steam_flow_kg_s_,
        (steam_pressure_MPa_ - P_cond_MPa) * 1e6f,
        rho_fw, h_fw, 0.82f); // η_pump = 0.82

    float pump_parasitic_MW = pump.W_pump_W * 1e-6f;

    // Rotor dynamics: J * dω/dt = T_turbine - T_generator - T_friction
    float omega_nom  = cfg_.nominal_rpm * 2.0f * 3.14159f / 60.0f; // rad/s
    float omega      = turbine_rpm_     * 2.0f * 3.14159f / 60.0f;
    float T_turb     = (omega > 1.0f) ? turb.W_shaft_W / omega : 0.0f;
    float T_load     = (omega > 1.0f) ? gross_electric_MW_ * 1e6f / omega : 0.0f;
    float T_friction = 0.001f * cfg_.turbine_inertia * omega;

    float alpha = (T_turb - T_load - T_friction) / cfg_.turbine_inertia;  // [rad/s²]
    float omega_new = omega + alpha * dt;
    turbine_rpm_    = std::max(omega_new * 60.0f / (2.0f * 3.14159f), 0.0f);

    // Condenser heat rejection (informational)
    SteamCyclePhysics::condenserHeatRejection(
        steam_flow_kg_s_, turb.h_out_actual_J_kg,
        P_cond_MPa, 290.0f, 305.0f);

    (void)pump_parasitic_MW;
    (void)state;
}

void SteamPowerSystem::updateGenerator(ReactorState& state, float dt)
{
    // Net electrical output
    float parasitics = cfg_.pump_load_MW + cfg_.magnet_load_MW + cfg_.aux_load_MW;
    net_electric_MW_ = gross_electric_MW_ - parasitics;
    net_electric_MW_ = std::max(net_electric_MW_, -parasitics); // can be negative at low power

    // Grid frequency deviation (50 Hz nominal → 3000 RPM for 2-pole)
    // Used as a diagnostic but not fed back here (grid control is external)
    float rpm_error = turbine_rpm_ - cfg_.nominal_rpm;
    float freq_Hz   = turbine_rpm_ / 60.0f; // for 1-pole pair
    (void)rpm_error; (void)freq_Hz; (void)state; (void)dt;
}