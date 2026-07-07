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
    //  temp_K and level_m are DOUBLE, not float: with a ~1.2e6 kg salt
    //  inventory and 1 ms ticks, the per-tick temperature increment
    //  (~9e-5 K) is smaller than the float ULP at ~860 K (~1e-4 K), so a
    //  float accumulator silently stops integrating — the hot tank froze
    //  at whatever temperature it first reached and "heating from core
    //  heat" appeared completely broken.
    double temp_K         = 700.0;
    double level_m        = 8.0;    // normal 8 m; max 16 m
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

    // Blanket circulation outlet temperature (hot leg into the hot tank) [K]
    float blanket_outlet_K  = 700.f;

    // Salt temperatures at SG inlets and outlets
    float sg_salt_inlet_K[4]  = {850.f,850.f,850.f,850.f};
    float sg_salt_outlet_K[4] = {800.f,800.f,800.f,800.f};
};

class MoltenSaltSystem {
public:
    MoltenSaltSystem();

    // blanket_heat_MW    = heat produced in blanket this tick (from ThermalHydraulics)
    // sg_steam_sat_K[4]  = saturation temperature of each SG's steam side [K]
    //                      (sets the pinch-point limit on how cold the salt
    //                      leaving that SG can get)
    // sg_demand_factor[4]= 0–1, how much of each SG's potential heat output
    //                      the secondary side can currently carry away
    //                      (TurbineUnitController::sgDemandFactor(), one tick
    //                      delayed — see main loop ordering)
    //
    // This is the salt/steam coupling: a steam generator is a heat
    // exchanger, and heat transfer is limited by BOTH sides. If the
    // turbine/bypass/relief on the steam side can't absorb heat (e.g.
    // turbine offline), the salt passes through largely unchanged and the
    // primary loop heats up instead of the SG cooling at a fixed rate
    // regardless of demand.
    void update(ReactorState& state, const SimTime& t, float blanket_heat_MW,
                const float sg_steam_sat_K[4], const float sg_demand_factor[4]);

    // Cold-restart: rebuild the salt state in place (tanks, pumps, alarms,
    // SG inlet/outlet temps) but keep the operator's pump enable/speed
    // settings — those are part of the plant configuration, not transient
    // state.  Required for RESET — COLD RESTART to actually deliver a
    // clean slate.
    void reset();

    const MoltenSaltState& saltState() const { return s_; }
    MoltenSaltState&       saltState()       { return s_; }

    // sg_heat_MW[4] array to feed to TurbineSystem
    const float* sgHeatMW() const { return s_.sg_heat_MW; }

private:
    void updateHotlegPumps(float dt);
    void updateColdlegPumps(float dt);
    void updateBlanketCirc(float dt, float blanket_heat_MW);
    void updateDistribution(float dt, const float sg_steam_sat_K[4],
                             const float sg_demand_factor[4]);
    void updateTankTemps(float dt);
    void checkAlarms();

    MoltenSaltState s_;
    static constexpr float CP_SALT = 2415.f;   // FLiBe cp [J/(kg·K)]
    static constexpr float T_FREEZE = 733.f;    // FLiBe freeze point [K]
    static constexpr float BCP_CAPACITY_KG_S    = 800.f;   // per blanket circ pump
    static constexpr float TANK_FLOOR_M2        = 80.f;    // salt tank floor area
    static constexpr float HOTLEG_CAPACITY_KG_S = 1500.f;  // per hotleg pump
    static constexpr float COLDLEG_CAPACITY_KG_S= 1400.f;  // per coldleg pump
};