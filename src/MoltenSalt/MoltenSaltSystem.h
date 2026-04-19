#pragma once
//
// src/MoltenSalt/MoltenSaltSystem.h
// FLiBe molten salt secondary loop:
//
//   Blanket (heat source)
//       ↓
//   [Hot Tank] ← hot salt arrives from blanket
//       ↓
//   [Hotleg Pumps 1A/1B] → [1&2 Distribution Valve] → split → SG1, SG2
//   [Hotleg Pumps 2A/2B] → [3&4 Distribution Valve] → split → SG3, SG4
//       ↓ (cooled salt leaves SGs)
//   [Coldleg Pumps 1–4] → [Cold Tank]
//       ↓
//   [Blanket Circulation Pumps] → back to blanket
//

#include "ReactorState.h"
#include "SimTime.h"

struct SaltPump {
    bool  running        = false;
    float speed_frac     = 0.f;     // 0–1
    float flow_kg_s      = 0.f;     // computed
    float power_MW       = 0.f;
    bool  trip           = false;
    bool  auto_start     = false;   // auto-start on demand signal
};

struct SaltTank {
    float temp_K         = 700.f;
    float level_m        = 8.f;     // normal 8 m; max 16 m
    float capacity_m3    = 2000.f;  // volume per metre of level (~250 m² floor)
    float salt_rho       = 1940.f;  // FLiBe density [kg/m³]
    bool  hi_level_alarm = false;
    bool  lo_level_alarm = false;
    bool  hi_temp_alarm  = false;
    bool  lo_temp_alarm  = false;   // FLiBe freezes at ~460°C (733 K) — must stay hot
};

struct MoltenSaltState {
    // Tanks
    SaltTank hot_tank;
    SaltTank cold_tank;

    // Hotleg pumps: two pairs (for 1&2 group and 3&4 group)
    SaltPump hotleg[4];   // [0]=1A, [1]=1B (feed group 12), [2]=2A, [3]=2B (feed group 34)

    // Coldleg return pumps: one per SG
    SaltPump coldleg[4];

    // Blanket circulation pumps (primary, interface with ThermalHydraulics)
    SaltPump blanket_circ[2];

    // Distribution valves (fraction sent to turbine 1 vs 2 within group, etc.)
    bool  dist_12_enabled   = true;    // enable flow to turbines 1 & 2
    bool  dist_34_enabled   = true;    // enable flow to turbines 3 & 4
    float dist_1_frac       = 0.5f;   // 0–1: fraction of 1&2 group flow to T1 (rest to T2)
    float dist_3_frac       = 0.5f;   // fraction of 3&4 group flow to T3 (rest to T4)

    // Computed: heat delivered to each SG [MW]
    float sg_heat_MW[4]     = {};

    // Salt temperatures at SG inlets and outlets
    float sg_salt_inlet_K[4]  = {850.f,850.f,850.f,850.f};
    float sg_salt_outlet_K[4] = {600.f,600.f,600.f,600.f};
};

class MoltenSaltSystem {
public:
    MoltenSaltSystem();

    // blanket_heat_MW = heat produced in blanket this tick (from ThermalHydraulics)
    void update(ReactorState& state, const SimTime& t, float blanket_heat_MW);

    const MoltenSaltState& saltState() const { return s_; }
    MoltenSaltState&       saltState()       { return s_; }

    // sg_heat_MW[4] array to feed to TurbineSystem
    const float* sgHeatMW() const { return s_.sg_heat_MW; }

private:
    void updateHotlegPumps(float dt);
    void updateColdlegPumps(float dt);
    void updateBlanketCirc(float dt, float blanket_heat_MW);
    void updateDistribution(float dt);
    void updateTankTemps(float dt);
    void checkAlarms();

    MoltenSaltState s_;
    static constexpr float CP_SALT = 2415.f;   // FLiBe cp [J/(kg·K)]
    static constexpr float T_FREEZE = 733.f;    // FLiBe freeze point [K]
};