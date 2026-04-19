//
// src/TurbineSystem/TurbineSystem.cpp
//
#include "TurbineSystem.h"

TurbineSystem::TurbineSystem()
    : units_{TurbineUnitController{1}, TurbineUnitController{2},
             TurbineUnitController{3}, TurbineUnitController{4}}
{}

void TurbineSystem::update(ReactorState& state, const SimTime& t,
                            float grid_frequency_Hz,
                            const float salt_heat_MW[4])
{
    float total_gross = 0.f;
    float total_parasitic = 0.f;

    for (int i = 0; i < 4; i++) {
        // Feed heat from molten salt to each steam generator
        units_[i].s.sg_heat_input_MW = salt_heat_MW[i];
        // Step physics
        units_[i].update(t.dt_s, grid_frequency_Hz);

        // Aggregate power
        total_gross     += units_[i].netPowerMW();
        // FW pump parasitic loads
        for (auto& p : units_[i].s.fw_pump) total_parasitic += p.power_MW;
        // Condensate pump ~ 2 MW each unit when running
        total_parasitic += units_[i].s.condenser.condensate_pump_speed * 2.f;
    }

    // Write summary to ReactorState (TurbineSystem overrides SteamPower legacy)
    state.gross_electric_MW = total_gross;
    // Use the first online unit's steam conditions for legacy display
    for (int i = 0; i < 4; i++) {
        if (units_[i].s.state == TurbineState::Online) {
            state.turbine_rpm         = units_[i].s.rpm;
            state.steam_pressure_MPa  = units_[i].s.sg_pressure_MPa;
            state.steam_temp_K        = units_[i].s.sg_steam_temp_K;
            break;
        }
    }
    (void)total_parasitic;
}

float TurbineSystem::totalPowerMW() const
{
    float sum = 0.f;
    for (int i = 0; i < 4; i++) sum += units_[i].netPowerMW();
    return sum;
}