//
// src/MoltenSalt/MoltenSaltSystem.cpp
//
#include "MoltenSaltSystem.h"
#include <cmath>
#include <algorithm>

MoltenSaltSystem::MoltenSaltSystem()
{
    //  FLiBe-compatible temperatures: the salt freezes at 733 K, so the
    //  COLD leg must run ~800 K (the old 610 K cold tank was 120 K below
    //  the freezing point — permanently frozen, with the lo-temp alarm
    //  latched on from the first frame).
    s_.hot_tank.temp_K   = 850.f;
    s_.cold_tank.temp_K  = 800.f;
    s_.hot_tank.level_m  = 8.f;
    s_.cold_tank.level_m = 8.f;

    // Pre-configure pumps
    for (auto& p : s_.hotleg)      { p.auto_start = true; }
    for (auto& p : s_.coldleg)     { p.auto_start = true; }
    for (auto& p : s_.blanket_circ){ p.auto_start = true; }
}

void MoltenSaltSystem::reset()
{
    // Preserve pump auto_start (it's a config flag) and any operator-set
    // running/speed — the operator shouldn't have to re-enable every pump
    // after pressing RESET.  Everything else (tank temps, levels, alarms,
    // SG heat duty) goes back to the construction defaults.
    auto save_pump = [](SaltPump& p) {
        bool was_running    = p.running;
        float was_speed     = p.speed_frac;
        bool was_auto       = p.auto_start;
        p = SaltPump{};
        p.running    = was_running;
        p.speed_frac = was_speed;
        p.auto_start = was_auto;
    };
    for (auto& p : s_.hotleg)       save_pump(p);
    for (auto& p : s_.coldleg)      save_pump(p);
    for (auto& p : s_.blanket_circ) save_pump(p);

    //  FLiBe-compatible temperatures: the salt freezes at 733 K, so the
    //  COLD leg must run ~800 K (the old 610 K cold tank was 120 K below
    //  the freezing point — permanently frozen, with the lo-temp alarm
    //  latched on from the first frame).
    s_.hot_tank.temp_K   = 850.f;
    s_.cold_tank.temp_K  = 800.f;
    s_.hot_tank.level_m  = 8.f;
    s_.cold_tank.level_m = 8.f;
    s_.hot_tank.hi_level_alarm  = false;
    s_.hot_tank.lo_level_alarm  = false;
    s_.hot_tank.hi_temp_alarm   = false;
    s_.hot_tank.lo_temp_alarm   = false;
    s_.cold_tank.hi_level_alarm = false;
    s_.cold_tank.lo_level_alarm = false;
    s_.cold_tank.lo_temp_alarm  = false;

    for (int i = 0; i < 4; i++) {
        s_.sg_heat_MW[i]       = 0.f;
        s_.sg_salt_inlet_K[i]  = 850.f;
        s_.sg_salt_outlet_K[i] = 800.f;
    }
}

// ─── Hotleg pump flow computation ─────────────────────────────────────────────
void MoltenSaltSystem::updateHotlegPumps(float dt)
{
    for (int i = 0; i < 4; i++) {
        auto& p = s_.hotleg[i];
        // ── AUTO: mass-balance against the blanket circulation ─────────────
        //  In steady state the flow leaving the hot tank (to the SGs) must
        //  match the flow arriving from the blanket, or the tank drains —
        //  the old auto ran the hotlegs flat out at 3000 kg/s against a
        //  ~1300 kg/s blanket input, letting the SGs "extract" 3+ GW from
        //  a 460 MW core by silently draining stored salt.  Each group
        //  (two pumps) targets half the blanket circ flow.
        if (p.auto_start && !p.trip) {
            bool group_enabled = (i < 2) ? s_.dist_12_enabled : s_.dist_34_enabled;
            float bcp_total = s_.blanket_circ[0].flow_kg_s + s_.blanket_circ[1].flow_kg_s;
            bool demand = group_enabled
                        && bcp_total > 10.f
                        && s_.hot_tank.level_m > 1.0f
                        && s_.hot_tank.temp_K  > T_FREEZE + 60.f;
            if (demand) {
                p.running    = true;
                // group target = bcp_total/2, split over two pumps
                p.speed_frac = std::clamp((bcp_total * 0.25f) / HOTLEG_CAPACITY_KG_S,
                                          0.02f, 1.f);
            } else {
                p.running = false;
            }
        }
        if (p.running && !p.trip) {
            // Two pumps per group; each at full speed delivers 1500 kg/s
            p.flow_kg_s = p.speed_frac * HOTLEG_CAPACITY_KG_S;
            p.power_MW  = p.flow_kg_s * 0.002f; // ~2 kJ/kg lift at pump head
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
        // ── Suction limit: an emptying hot tank starves the pumps ──────────
        //  Manual over-pumping now drains the tank and loses flow instead
        //  of conjuring salt (and heat) from nowhere.
        float avail = (float)std::clamp(s_.hot_tank.level_m, 0.0, 1.0);
        p.flow_kg_s *= avail;
    }
    (void)dt;
}

void MoltenSaltSystem::updateColdlegPumps(float dt)
{
    for (int i = 0; i < 4; i++) {
        auto& p = s_.coldleg[i];
        // ── AUTO: match the salt flow arriving at this SG so the return
        //    leg keeps pace with the supply (level-neutral operation).
        if (p.auto_start && !p.trip) {
            float supply = 0.f;
            if (i < 2) supply = (s_.dist_12_enabled ? (s_.hotleg[0].flow_kg_s + s_.hotleg[1].flow_kg_s) : 0.f)
                              * (i == 0 ? s_.dist_1_frac : (1.f - s_.dist_1_frac));
            else       supply = (s_.dist_34_enabled ? (s_.hotleg[2].flow_kg_s + s_.hotleg[3].flow_kg_s) : 0.f)
                              * (i == 2 ? s_.dist_3_frac : (1.f - s_.dist_3_frac));
            if (supply > 10.f) {
                p.running    = true;
                p.speed_frac = std::clamp(supply / COLDLEG_CAPACITY_KG_S, 0.05f, 1.f);
            } else {
                p.running = false;
            }
        }
        if (p.running && !p.trip) {
            p.flow_kg_s = p.speed_frac * COLDLEG_CAPACITY_KG_S;
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
    // ── AUTO mode: hold ~900 K blanket outlet temperature ──────────────────
    //  A real salt plant modulates primary flow to hold the hot-leg
    //  temperature: too much flow and the salt barely warms (the old bug —
    //  6000 kg/s of pump capacity gave a ~28 K ΔT, so "hot" salt arrived
    //  COLDER than the 860 K tank and running the reactor cooled it);
    //  too little flow and the salt overheats toward material limits.
    //  In AUTO the required total flow is ṁ = Q / (cp · (T_target − T_cold)),
    //  split across the two pumps.
    constexpr float T_OUTLET_TARGET_K = 950.f;
    bool any_auto = s_.blanket_circ[0].auto_start || s_.blanket_circ[1].auto_start;
    if (any_auto && blanket_heat_MW > 2.f) {
        float dT_target  = (float)std::max((double)T_OUTLET_TARGET_K - s_.cold_tank.temp_K, 30.0);
        float flow_needed = blanket_heat_MW * 1e6f / (CP_SALT * dT_target);
        int n_auto = (s_.blanket_circ[0].auto_start ? 1 : 0)
                   + (s_.blanket_circ[1].auto_start ? 1 : 0);
        for (auto& p : s_.blanket_circ) {
            if (!p.auto_start || p.trip) continue;
            p.running    = true;
            p.speed_frac = std::clamp(flow_needed / (n_auto * BCP_CAPACITY_KG_S),
                                      0.05f, 1.f);
        }
    } else if (any_auto) {
        // No blanket heat — auto pumps idle at low speed to keep salt moving
        // (FLiBe must not stagnate and freeze in the lines)
        for (auto& p : s_.blanket_circ)
            if (p.auto_start && !p.trip) { p.running = true; p.speed_frac = 0.05f; }
    }

    for (int i = 0; i < 2; i++) {
        auto& p = s_.blanket_circ[i];
        if (p.running && !p.trip) {
            //  Capacity: 500 kg/s per pump (1000 kg/s total).  At full
            //  blanket power (~460 MW) full flow gives ΔT ≈ 190 K; the
            //  design point for a 900 K hot leg from a 610 K cold tank
            //  (ΔT = 290 K) is ~660 kg/s — i.e. ~65% pump speed.  Throttle
            //  the pumps to trade outlet temperature against flow, exactly
            //  like a real salt plant.
            p.flow_kg_s = p.speed_frac * BCP_CAPACITY_KG_S;
            p.power_MW  = p.flow_kg_s * 0.003f;
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
    }

    // Blanket heats salt: ΔT = Q / (ṁ * cp)
    float total_flow = s_.blanket_circ[0].flow_kg_s + s_.blanket_circ[1].flow_kg_s;
    if (total_flow > 1.f) {
        float dT = blanket_heat_MW * 1e6f / (total_flow * CP_SALT);
        // Salt enters blanket from cold tank, exits to hot tank
        float T_in    = (float)s_.cold_tank.temp_K;
        float T_out   = std::min(T_in + dT, 1300.f);  // FLiBe limit clamp
        s_.blanket_outlet_K = T_out;
        // Hot tank receives heated salt
        float dHotLevel = total_flow * dt / (s_.hot_tank.salt_rho * TANK_FLOOR_M2); // tank floor area
        s_.hot_tank.level_m  += dHotLevel;
        s_.cold_tank.level_m -= dHotLevel;
        // Hot tank temperature: mix of arriving hot salt with existing
        float m_tank   = s_.hot_tank.level_m * TANK_FLOOR_M2 * s_.hot_tank.salt_rho;
        float m_arrive = total_flow * dt;
        if (m_tank > 1.f)
            s_.hot_tank.temp_K = (s_.hot_tank.temp_K * m_tank + T_out * m_arrive)
                               / (m_tank + m_arrive);
    } else {
        s_.blanket_outlet_K = (float)s_.cold_tank.temp_K;
    }

    // ── Tank heat losses ────────────────────────────────────────────────────
    //  Even well-insulated salt tanks lose heat (~1 MW-class for tanks of
    //  this size).  Without this, a hot tank stayed at 860 K forever with
    //  the reactor cold — now trace-heating/keeping the plant running
    //  actually matters, and the freeze alarm is a genuine threat during
    //  long outages.
    auto tank_loss = [&](SaltTank& tk) {
        float m = tk.level_m * TANK_FLOOR_M2 * tk.salt_rho;
        if (m < 1.f) return;
        constexpr float LOSS_W_PER_K = 3000.f;   // UA to ambient
        float Q_loss = LOSS_W_PER_K * (tk.temp_K - 300.f);
        tk.temp_K -= Q_loss * dt / (m * CP_SALT);
    };
    tank_loss(s_.hot_tank);
    tank_loss(s_.cold_tank);
}

// ─── Distribution and SG heat calculation ────────────────────────────────────
void MoltenSaltSystem::updateDistribution(float dt, const float sg_steam_sat_K[4],
                                           const float sg_demand_factor[4])
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

    // Each SG is a heat exchanger: heat transfer is limited by BOTH sides.
    //  - Salt side sets an upper bound on how much heat *could* be given up:
    //    the salt can't be cooled below the steam-side saturation
    //    temperature plus a pinch-point approach (2nd law).
    //  - Steam side sets how much of that potential heat can actually be
    //    *carried away* right now (sgDemandFactor: turbine load, bypass,
    //    relief venting). If the secondary side can't take steam, the salt
    //    passes through nearly unchanged and the primary loop heats up
    //    instead — this is what makes turbine availability matter for the
    //    molten salt loop, rather than the SG always cooling by a fixed ΔT.
    constexpr float APPROACH_K = 30.f; // minimum pinch-point approach temp

    for (int i = 0; i < 4; i++) {
        s_.sg_salt_inlet_K[i] = (float)s_.hot_tank.temp_K;

        //  Anti-freeze floor: FLiBe leaving the SG must stay comfortably
        //  above its 733 K freezing point regardless of how cold the steam
        //  side is — a real salt SG is designed so the salt-side outlet
        //  cannot approach freezing (frozen salt in a heat exchanger is a
        //  plant-killing event).  Without this floor the old model happily
        //  cooled the salt to T_sat + 30 ≈ 600 K — solid FLiBe.
        float T_min   = std::max(sg_steam_sat_K[i] + APPROACH_K, T_FREEZE + 40.f);
        float demand  = std::clamp(sg_demand_factor[i], 0.f, 1.f);

        float T_out;
        if (flow[i] > 10.f) {
            float coolable = (float)std::max(s_.hot_tank.temp_K - (double)T_min, 0.0);
            T_out = (float)s_.hot_tank.temp_K - demand * coolable;
        } else {
            T_out = (float)s_.hot_tank.temp_K; // no flow, no cooling
        }
        s_.sg_salt_outlet_K[i] = T_out;

        // Heat transferred [MW] = ṁ * cp * ΔT
        float dT = (float)s_.hot_tank.temp_K - T_out;
        s_.sg_heat_MW[i] = (flow[i] * CP_SALT * dT) * 1e-6f;

        // Return cooled salt to cold tank via coldleg pumps
        // (simplified: coldleg flow matches hotleg flow per SG)
        float coldleg_flow = std::min(s_.coldleg[i].flow_kg_s, flow[i]);
        float returned     = coldleg_flow * dt;
        float dColdLevel   = returned / (s_.cold_tank.salt_rho * TANK_FLOOR_M2);
        s_.cold_tank.level_m  += dColdLevel;
        s_.hot_tank.level_m   -= returned / (s_.hot_tank.salt_rho * TANK_FLOOR_M2);

        // Cold tank temperature: mix of returning cool salt
        //  m_tank is the mass of salt currently in the cold tank [kg].
        //  m_return is the mass of cooled salt returning from this SG [kg].
        //  The old code had `m_return = returned * salt_rho` which
        //  multiplied kg by kg/m³ = kg²/m³ (wrong units).  `returned` is
        //  already in kg (coldleg_flow [kg/s] × dt [s]), so m_return = returned.
        float m_tank    = s_.cold_tank.level_m * TANK_FLOOR_M2 * s_.cold_tank.salt_rho;
        float m_return  = returned;  // already in kg
        if (m_tank > 1.f)
            s_.cold_tank.temp_K = (s_.cold_tank.temp_K * m_tank + T_out * m_return)
                                / (m_tank + m_return);
    }

    // Clamp tank levels
    s_.hot_tank.level_m  = std::clamp(s_.hot_tank.level_m,  0.0, 16.0);
    s_.cold_tank.level_m = std::clamp(s_.cold_tank.level_m, 0.0, 16.0);
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
    ht.lo_temp_alarm  = (ht.temp_K  < T_FREEZE + 30.f); // freeze risk

    ct.hi_level_alarm = (ct.level_m > 14.f);
    ct.lo_level_alarm = (ct.level_m < 2.f);
    ct.lo_temp_alarm  = (ct.temp_K  < T_FREEZE + 30.f);
}

// ─── Main update ─────────────────────────────────────────────────────────────
void MoltenSaltSystem::update(ReactorState& state, const SimTime& t,
                               float blanket_heat_MW,
                               const float sg_steam_sat_K[4],
                               const float sg_demand_factor[4])
{
    float dt = t.dt_s;
    updateBlanketCirc(dt, blanket_heat_MW);
    updateHotlegPumps(dt);
    updateColdlegPumps(dt);
    updateDistribution(dt, sg_steam_sat_K, sg_demand_factor);
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