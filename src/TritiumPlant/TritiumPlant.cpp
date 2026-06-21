//
// src/TritiumPlant/TritiumPlant.cpp
//  Tritium extraction, isotope separation, and detritiation.
//

#include "TritiumPlant.h"
#include <algorithm>
#include <cmath>

TritiumPlantSystem::TritiumPlantSystem(const TritiumPlantConfig& cfg)
    : cfg_(cfg)
    , T_in_plant_g_(cfg.initial_T_in_plant_g)
    , T_in_storage_g_(0.f)
{}

void TritiumPlantSystem::reset()
{
    T_in_plant_g_ = cfg_.initial_T_in_plant_g;
    T_in_storage_g_ = 0.f;
    tes_on_ = false;
    tes_rate_g_s_ = 0.f;
    detritiation_on_ = false;
    detritiation_flow_m3_s_ = 0.f;
}

void TritiumPlantSystem::startDetritiation(float flow_m3_s)
{
    detritiation_on_ = true;
    detritiation_flow_m3_s_ = (flow_m3_s < 0.f) ? cfg_.detritiation_max_m3_s
                                                : std::min(flow_m3_s, cfg_.detritiation_max_m3_s);
}

void TritiumPlantSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // ── TES: extract T from blanket → fuel store ────────────────────────────
    //  Extraction rate scales with blanket coolant flow (more Li in flow →
    //  more T extracted).  T_in_plant_g decreases, fuel_T_inventory_g
    //  increases (the Fuel module reads this back via state.tritium_recovery_rate_g_s).
    if (tes_on_ && T_in_plant_g_ > 0.f) {
        float flow_factor = std::clamp(state.coolant_flow_kg_s / 18000.f, 0.f, 1.f);
        tes_rate_g_s_ = cfg_.tes_rate_g_per_s_per_kg_s
                       * state.coolant_flow_kg_s * flow_factor;
        float extracted = std::min(tes_rate_g_s_ * dt, T_in_plant_g_);
        T_in_plant_g_ -= extracted;
        // The Fuel module will pick this up via state.tritium_recovery_rate_g_s
        // and add it to T_inventory_g_ on its next update.
        state.tritium_recovery_rate_g_s = tes_rate_g_s_;
    } else {
        tes_rate_g_s_ = 0.f;
        state.tritium_recovery_rate_g_s = 0.f;
    }

    // ── Tritium breeding: TBR > 1 produces net new T ────────────────────────
    //  The blanket breeds T via Li-6(n,t)He-4.  ThermalHydraulics computes
    //  TBR; we just convert that to a T production rate and add to plant.
    if (state.plasma_status != PlasmaStatus::Cold && state.fusion_power_MW > 0.f) {
        // T consumed by fusion: 1 T per D-T reaction
        constexpr float E_fusion_J = 17.59e6f * 1.602176634e-19f;
        float T_consumed_g_s = state.fusion_power_MW * 1e6f / E_fusion_J
                              * 5.00735588e-27f * 1e3f;  // T mass in g
        // T bred: TBR * T_consumed
        float T_bred_g_s = state.tbr_current * T_consumed_g_s;
        T_in_plant_g_ += T_bred_g_s * dt;
        // Cap at storage capacity
        T_in_plant_g_ = std::min(T_in_plant_g_, cfg_.tes_max_g * 2.f);
    }

    // ── Detritiation ─────────────────────────────────────────────────────────
    //  When running, removes tritium from the building air.  This is mainly
    //  a safety system — it doesn't directly affect the plant inventory,
    //  but it does keep the building habitable after a tritium release.
    state.detritiation_flow_m3_s = detritiation_on_ ? detritiation_flow_m3_s_ : 0.f;

    // Write outputs
    state.tritium_in_plant_g = T_in_plant_g_;
    state.tritium_recovery_on = tes_on_;
}
