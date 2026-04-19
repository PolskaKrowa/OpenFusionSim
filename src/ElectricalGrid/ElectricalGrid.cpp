//
// src/ElectricalGrid/ElectricalGrid.cpp
//
#include "ElectricalGrid.h"
#include "TurbineSystem/TurbineSystem.h"
#include <cmath>
#include <algorithm>
#include <numeric>

static constexpr float PI       = 3.14159265f;
static constexpr float FREQ_NOM = 50.f;

ElectricalGridSystem::ElectricalGridSystem()
{
    // Define on-site load groups
    site_loads_ = {{
        {"Plasma Control Systems",   2.5f,  true, true},
        {"Cryoplant (magnets)",      35.f,  true, true},
        {"He Reactor Cooling",       4.0f,  true, true},
        {"Molten Salt Pumps",        8.0f,  true, true},
        {"TF Coil Power Supplies",  12.0f,  true, true},
        {"FW Pumps (all units)",     8.0f,  true, true},
        {"Condenser CAR/SJAE",       3.0f,  true, false},
        {"HVAC & Lighting",          2.0f,  true, false},
        {"Tritium Plant",            5.0f,  true, false},
        {"Water Treatment",          1.5f,  true, false},
        {"Control Building",         0.8f,  true, true},
        {"Spare/Test Load",          0.0f,  false,false},
    }};
}

void ElectricalGridSystem::update(ReactorState& state, const SimTime& t,
                                   const TurbineSystem& turbines)
{
    float dt = t.dt_s;

    // Collect generator bus data from turbine units
    float total_gen = 0.f;
    for (int i = 0; i < 4; i++) {
        const auto& u = turbines.unit(i).s;
        grid_.bus[i].frequency_Hz   = u.generator.frequency_Hz;
        grid_.bus[i].phase_rad      = u.generator.phase_rad;
        grid_.bus[i].power_MW       = u.generator.power_MW;
        grid_.bus[i].breaker_closed = u.generator.breaker_closed;
        if (u.generator.breaker_closed) total_gen += u.generator.power_MW;
    }

    updateSiteLoads();
    updateFrequency(dt, total_gen);
    updateSynchCheck(turbines);
    updateExportImport();

    // Write summary to ReactorState
    state.grid_frequency_Hz = grid_.frequency_Hz;
    state.grid_connected    = (total_gen > 0.f);
    state.gross_electric_MW = total_gen;
    state.net_electric_MW   = total_gen - grid_.total_site_load_MW;
}

void ElectricalGridSystem::updateFrequency(float dt, float total_gen_MW)
{
    // If loss of offsite power, frequency is purely our generators vs site load
    // Otherwise, infinite bus holds frequency near 50 Hz with droop
    float load = grid_.total_site_load_MW;

    if (grid_.loss_of_offsite_power) {
        // Island mode: swing equation — only our generator inertia
        float H_total = 6.f * 4.f; // 4 generators × H=6 s each (rough)
        float P_base  = 4800.f;    // MW
        float df_dt   = (total_gen_MW - load) * FREQ_NOM / (2.f * H_total * P_base);
        grid_.frequency_Hz += df_dt * dt;
        grid_.frequency_Hz  = std::clamp(grid_.frequency_Hz, 40.f, 55.f);
    } else {
        // Grid-connected: external bus pulls frequency toward 50 Hz
        float df_dt = (FREQ_NOM - grid_.frequency_Hz) * 2.f;
        grid_.frequency_Hz += df_dt * dt;
        grid_.frequency_Hz  = std::clamp(grid_.frequency_Hz, 48.f, 52.f);
    }

    grid_.frequency_deviation_Hz = grid_.frequency_Hz - FREQ_NOM;
    grid_.underfrequency_alarm   = (grid_.frequency_Hz < 49.f);
    grid_.overfrequency_alarm    = (grid_.frequency_Hz > 51.f);
    grid_.underfrequency_trip    = (grid_.frequency_Hz < 47.5f);
}

void ElectricalGridSystem::updateSynchCheck(const TurbineSystem& turbines)
{
    // Check synchronisation conditions for each generator
    for (int i = 0; i < 4; i++) {
        const auto& u = turbines.unit(i).s;
        float df    = std::abs(u.generator.frequency_Hz - grid_.frequency_Hz);
        // Phase difference (wrap)
        float dp    = std::fmodf(std::abs(u.generator.phase_rad
                                  - grid_.frequency_Hz * 2.f * PI * 0.f), 2.f * PI);
        dp = std::min(dp, 2.f * PI - dp);
        // Synch OK: Δf < 0.15 Hz and Δφ < 15°
        grid_.bus[i].synch_ok = (df < 0.15f && dp < 0.26f &&
                                  u.rpm > 2900.f && u.state == TurbineState::Synchronizing);
    }
}

void ElectricalGridSystem::updateSiteLoads()
{
    float total = 0.f;
    for (auto& L : site_loads_) {
        if (L.energised) total += L.load_MW;
    }
    grid_.total_site_load_MW = total;
}

void ElectricalGridSystem::updateExportImport()
{
    float net = grid_.total_generation_MW - grid_.total_site_load_MW;
    grid_.export_MW = std::max(net, 0.f);
    grid_.import_MW = std::max(-net, 0.f);
}

void ElectricalGridSystem::setLoadEnergised(int idx, bool on)
{
    if (idx >= 0 && idx < (int)site_loads_.size())
        if (!site_loads_[idx].essential || on) // can't de-energise essential loads
            site_loads_[idx].energised = on;
}

bool ElectricalGridSystem::canSync(int gen_idx) const
{
    return grid_.bus[gen_idx].synch_ok;
}