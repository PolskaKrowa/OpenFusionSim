//
// src/MoltenSalt/MoltenSaltSystem.cpp
//
#include "MoltenSaltSystem.h"
#include <cmath>
#include <algorithm>

MoltenSaltSystem::MoltenSaltSystem()
{
    s_.hot_tank.temp_K   = 860.f;
    s_.cold_tank.temp_K  = 610.f;
    s_.hot_tank.level_m  = 8.f;
    s_.cold_tank.level_m = 8.f;

    // Pre-configure pumps
    for (auto& p : s_.hotleg)      { p.auto_start = true; }
    for (auto& p : s_.coldleg)     { p.auto_start = true; }
    for (auto& p : s_.blanket_circ){ p.auto_start = true; }
}

// ─── Hotleg pump flow computation ─────────────────────────────────────────────
void MoltenSaltSystem::updateHotlegPumps(float dt)
{
    for (int i = 0; i < 4; i++) {
        auto& p = s_.hotleg[i];
        if (p.running && !p.trip) {
            // Two pumps per group; each at full speed delivers 1500 kg/s
            p.flow_kg_s = p.speed_frac * 1500.f;
            p.power_MW  = p.flow_kg_s * 0.002f; // ~2 kJ/kg lift at pump head
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
    }
    (void)dt;
}

void MoltenSaltSystem::updateColdlegPumps(float dt)
{
    for (int i = 0; i < 4; i++) {
        auto& p = s_.coldleg[i];
        if (p.running && !p.trip) {
            p.flow_kg_s = p.speed_frac * 1400.f;
            p.power_MW  = p.flow_kg_s * 0.0015f;
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
    }
    (void)dt;
}

void MoltenSaltSystem::updateBlanketCirc(float dt, float blanket_heat_MW)
{
    for (int i = 0; i < 2; i++) {
        auto& p = s_.blanket_circ[i];
        if (p.running && !p.trip) {
            p.flow_kg_s = p.speed_frac * 3000.f; // total blanket circulation
            p.power_MW  = p.flow_kg_s * 0.003f;
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
    }

    // Blanket heats salt: ΔT = Q / (ṁ * cp)
    float total_flow = s_.blanket_circ[0].flow_kg_s + s_.blanket_circ[1].flow_kg_s;
    if (total_flow > 0.f) {
        float dT = blanket_heat_MW * 1e6f / (total_flow * CP_SALT);
        // Salt enters blanket from cold tank, exits to hot tank
        float T_in    = s_.cold_tank.temp_K;
        float T_out   = T_in + dT;
        // Hot tank receives heated salt
        float dHotLevel = total_flow * dt / (s_.hot_tank.salt_rho * 250.f); // 250 m² floor
        s_.hot_tank.level_m  += dHotLevel;
        s_.cold_tank.level_m -= dHotLevel;
        // Hot tank temperature: mix of arriving hot salt with existing
        float m_tank   = s_.hot_tank.level_m * 250.f * s_.hot_tank.salt_rho;
        float m_arrive = total_flow * dt;
        if (m_tank > 1.f)
            s_.hot_tank.temp_K = (s_.hot_tank.temp_K * m_tank + T_out * m_arrive)
                               / (m_tank + m_arrive);
    }
    (void)dt;
}

// ─── Distribution and SG heat calculation ────────────────────────────────────
void MoltenSaltSystem::updateDistribution(float dt)
{
    // Group 1&2: hotleg pumps 0 and 1
    float flow_12 = s_.dist_12_enabled
                  ? (s_.hotleg[0].flow_kg_s + s_.hotleg[1].flow_kg_s) : 0.f;
    // Group 3&4: hotleg pumps 2 and 3
    float flow_34 = s_.dist_34_enabled
                  ? (s_.hotleg[2].flow_kg_s + s_.hotleg[3].flow_kg_s) : 0.f;

    // Split within each group
    float flow[4];
    flow[0] = flow_12 * s_.dist_1_frac;
    flow[1] = flow_12 * (1.f - s_.dist_1_frac);
    flow[2] = flow_34 * s_.dist_3_frac;
    flow[3] = flow_34 * (1.f - s_.dist_3_frac);

    // Each SG: hot salt cools from hot_tank temp to some outlet temp
    // Outlet temperature set by SG duty (simplified: fixed approach to steam temp)
    for (int i = 0; i < 4; i++) {
        s_.sg_salt_inlet_K[i] = s_.hot_tank.temp_K;
        // SG outlet: cooled by ~200 K when at full flow (approach to steam saturation)
        float T_steam  = 560.f; // K ~18 MPa sat temp approximation
        float approach = 30.f;  // minimum approach temp
        float T_out    = (flow[i] > 10.f)
                       ? std::max(T_steam + approach, s_.hot_tank.temp_K - 220.f)
                       : s_.hot_tank.temp_K; // no cooling if no flow
        s_.sg_salt_outlet_K[i] = T_out;

        // Heat transferred [MW] = ṁ * cp * ΔT
        float dT = s_.hot_tank.temp_K - T_out;
        s_.sg_heat_MW[i] = (flow[i] * CP_SALT * dT) * 1e-6f;

        // Return cooled salt to cold tank via coldleg pumps
        // (simplified: coldleg flow matches hotleg flow per SG)
        float coldleg_flow = std::min(s_.coldleg[i].flow_kg_s, flow[i]);
        float returned     = coldleg_flow * dt;
        float dColdLevel   = returned / (s_.cold_tank.salt_rho * 250.f);
        s_.cold_tank.level_m  += dColdLevel;
        s_.hot_tank.level_m   -= returned / (s_.hot_tank.salt_rho * 250.f);

        // Cold tank temperature: mix of returning cool salt
        float m_tank    = s_.cold_tank.level_m * 250.f * s_.cold_tank.salt_rho;
        float m_return  = returned * s_.cold_tank.salt_rho;
        if (m_tank > 1.f)
            s_.cold_tank.temp_K = (s_.cold_tank.temp_K * m_tank + T_out * m_return)
                                / (m_tank + m_return);
    }

    // Clamp tank levels
    s_.hot_tank.level_m  = std::clamp(s_.hot_tank.level_m,  0.f, 16.f);
    s_.cold_tank.level_m = std::clamp(s_.cold_tank.level_m, 0.f, 16.f);
    (void)dt;
}

// ─── Tank alarms ─────────────────────────────────────────────────────────────
void MoltenSaltSystem::checkAlarms()
{
    auto& ht = s_.hot_tank;
    auto& ct = s_.cold_tank;

    ht.hi_level_alarm = (ht.level_m > 14.f);
    ht.lo_level_alarm = (ht.level_m < 2.f);
    ht.hi_temp_alarm  = (ht.temp_K  > 950.f);  // approaching material limits
    ht.lo_temp_alarm  = (ht.temp_K  < T_FREEZE + 50.f); // freeze risk

    ct.hi_level_alarm = (ct.level_m > 14.f);
    ct.lo_level_alarm = (ct.level_m < 2.f);
    ct.lo_temp_alarm  = (ct.temp_K  < T_FREEZE + 50.f);
}

// ─── Main update ─────────────────────────────────────────────────────────────
void MoltenSaltSystem::update(ReactorState& state, const SimTime& t,
                               float blanket_heat_MW)
{
    float dt = t.dt_s;
    updateBlanketCirc(dt, blanket_heat_MW);
    updateHotlegPumps(dt);
    updateColdlegPumps(dt);
    updateDistribution(dt);
    checkAlarms();

    // Write to ReactorState summary
    state.hot_tank_temp_K    = s_.hot_tank.temp_K;
    state.cold_tank_temp_K   = s_.cold_tank.temp_K;
    state.hot_tank_level_m   = s_.hot_tank.level_m;
    state.cold_tank_level_m  = s_.cold_tank.level_m;

    float total_flow = 0.f;
    for (auto& p : s_.hotleg) total_flow += p.flow_kg_s;
    state.salt_flow_total_kg_s = total_flow;
}