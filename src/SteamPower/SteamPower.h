#pragma once

//
// src/SteamPower/SteamPower.h
// Rankine cycle power conversion: steam generator → turbine → condenser → generator.
//

#include "ReactorState.h"
#include "SimTime.h"

struct SteamConfig {
    // Steam conditions (supercritical water, DEMO-class)
    float nominal_pressure_MPa = 18.0f;   // live steam pressure [MPa]
    float nominal_temp_K       = 820.0f;  // live steam temperature [K]
    float nominal_flow_kg_s    = 3000.0f; // steam mass flow [kg/s]

    // Turbine
    float turbine_efficiency   = 0.88f;   // isentropic efficiency
    float generator_efficiency = 0.985f;  // electrical efficiency
    float nominal_rpm          = 3000.0f; // 50 Hz grid (or 3600 for 60 Hz)
    float turbine_inertia      = 5e5f;    // rotor moment of inertia [kg·m^2]

    // Parasitic loads
    float pump_load_MW         = 30.0f;   // coolant pumps
    float magnet_load_MW       = 50.0f;   // superconducting magnet cryo + power
    float aux_load_MW          = 20.0f;   // control, HVAC, misc
};

class SteamPowerSystem {
public:
    explicit SteamPowerSystem(const SteamConfig& cfg);

    void update(ReactorState& state, const SimTime& t);

    float netElectricMW() const { return net_electric_MW_; }

private:
    void updateSteamGenerator(ReactorState& state, float dt);
    void updateTurbine        (ReactorState& state, float dt);
    void updateGenerator      (ReactorState& state, float dt);

    // Approximate Rankine cycle thermal efficiency from steam conditions
    float rankineEfficiency(float T_hot_K, float T_cold_K) const;

    SteamConfig cfg_;

    float steam_pressure_MPa_  = 0.1f;
    float steam_temp_K_        = 300.0f;
    float steam_flow_kg_s_     = 0.0f;
    float turbine_rpm_         = 0.0f;
    float gross_electric_MW_   = 0.0f;
    float net_electric_MW_     = 0.0f;
    float thermal_input_MW_    = 0.0f;
};