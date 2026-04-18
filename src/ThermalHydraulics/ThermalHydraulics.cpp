//
// src/ThermalHydraulics/ThermalHydraulics.cpp
//

#include "ThermalHydraulics.h"
#include <algorithm>
#include <cmath>

ThermalHydraulics::ThermalHydraulics(const ThermalHydraulicsConfig& cfg)
    : cfg_(cfg)
    , blanket_temp_K_(cfg.coolant == CoolantType::FLiBe ? 700.0f : 573.0f)
    , coolant_inlet_K_(blanket_temp_K_)
    , coolant_outlet_K_(blanket_temp_K_ + 200.0f)
{}

void ThermalHydraulics::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    updateFirstWall   (state, dt);
    updateBlanket     (state, dt);
    updateCoolantLoop (state, dt);

    // Write outputs
    state.coolant_inlet_temp_K  = coolant_inlet_K_;
    state.coolant_outlet_temp_K = coolant_outlet_K_;
    state.coolant_flow_kg_s     = flow_kg_s_;
    state.blanket_heat_MW       = thermal_power_MW_;
    state.tbr_current           = tbr_current_;
    state.first_wall_temp_K     = first_wall_temp_K_;
    state.thermal_runaway       = (blanket_temp_K_ > cfg_.max_coolant_temp_K);
    state.alarm_overtemp       |= state.thermal_runaway;
    state.alarm_loss_of_coolant = (flow_kg_s_ < 100.0f && thermal_power_MW_ > 10.0f);
}

void ThermalHydraulics::updateFirstWall(ReactorState& state, float dt)
{
    // First wall absorbs ~3 % of fusion power as surface heat flux (alpha + radiation)
    float P_surface_MW = state.fusion_power_MW * 0.03f + state.radiated_power_MW * 0.2f;
    float q_dot        = P_surface_MW * 1e6f / cfg_.first_wall_area_m2; // W/m^2

    // Tungsten armour: lumped thermal capacity
    constexpr float rho_W = 19300.0f, Cp_W = 134.0f, thick = 0.02f; // 20 mm armour
    float C_wall = rho_W * Cp_W * thick * cfg_.first_wall_area_m2;    // J/K

    // Coolant back-side removes heat (convective)
    float h_conv  = 50000.0f; // W/(m^2·K) — helium-cooled first wall
    float Q_in    = q_dot * cfg_.first_wall_area_m2;
    float Q_out   = h_conv * cfg_.first_wall_area_m2 * (first_wall_temp_K_ - coolant_inlet_K_);

    first_wall_temp_K_ += (Q_in - Q_out) / C_wall * dt;
    first_wall_temp_K_  = std::max(first_wall_temp_K_, coolant_inlet_K_);
    (void)state;
}

void ThermalHydraulics::updateBlanket(ReactorState& state, float dt)
{
    if (state.neutron_flux_m2s <= 0.0f) {
        thermal_power_MW_ = 0.0f;
        return;
    }

    // 14 MeV neutrons deposit most of their energy in the Li blanket
    // P_blanket ≈ (14.07/17.59) * P_fusion * energy_mult  (80 % + amplification)
    float P_neutron_MW = state.fusion_power_MW * (14.07f / 17.59f);
    thermal_power_MW_  = P_neutron_MW * cfg_.blanket_energy_mult;

    // Blanket temperature: energy in from neutrons, energy out to coolant
    float blanket_mass = cfg_.rho_kg_m3 * cfg_.blanket_volume_m3;
    float C_blanket    = blanket_mass * cfg_.Cp_J_kgK;

    float Q_in  = thermal_power_MW_ * 1e6f;
    float Q_out = cfg_.HX_UA_W_K * (blanket_temp_K_ - coolant_inlet_K_)
                * std::min(flow_kg_s_ / cfg_.max_flow_kg_s, 1.0f);

    blanket_temp_K_ += (Q_in - Q_out) / C_blanket * dt;
    blanket_temp_K_  = std::max(blanket_temp_K_, coolant_inlet_K_);

    // TBR: depends on Li-6 enrichment, blanket thickness, and neutron energy spectrum
    tbr_current_ = computeTBR(state.neutron_flux_m2s);
    (void)dt;
}

void ThermalHydraulics::updateCoolantLoop(ReactorState& state, float dt)
{
    // Coolant flow tracks setpoint
    float target_flow = state.sp_coolant_flow * cfg_.max_flow_kg_s;
    float ramp        = cfg_.max_flow_kg_s * 0.5f; // kg/s per second ramp
    flow_kg_s_ += std::clamp(target_flow - flow_kg_s_, -ramp * dt, ramp * dt);
    flow_kg_s_  = std::clamp(flow_kg_s_, 0.0f, cfg_.max_flow_kg_s);

    if (flow_kg_s_ < 1.0f) {
        // No flow — blanket just heats up (loss of flow accident)
        coolant_outlet_K_ = blanket_temp_K_;
        return;
    }

    // Outlet temperature: coolant carries away Q_out from blanket
    float Q_extracted = cfg_.HX_UA_W_K * (blanket_temp_K_ - coolant_inlet_K_)
                      * std::min(flow_kg_s_ / cfg_.max_flow_kg_s, 1.0f);
    float delta_T = Q_extracted / (flow_kg_s_ * cfg_.Cp_J_kgK);
    coolant_outlet_K_ = coolant_inlet_K_ + delta_T;

    // Simple steam-generator model: HX returns cooled fluid to inlet
    // (SteamPower reads outlet temp; it sets the turbine conditions)
    float T_secondary = 500.0f; // secondary-side saturation temp [K]
    float Q_to_steam  = cfg_.HX_UA_W_K * 0.5f * (coolant_outlet_K_ - T_secondary);
    coolant_inlet_K_  = coolant_outlet_K_ - Q_to_steam / (flow_kg_s_ * cfg_.Cp_J_kgK);
    coolant_inlet_K_  = std::max(coolant_inlet_K_, T_secondary + 10.0f);
    (void)dt;
    (void)state;
}

float ThermalHydraulics::computeTBR(float neutron_flux) const
{
    // Simplified TBR model: scales with Li-6 fraction and blanket coverage
    // Reference: TBR ≈ 1.0 for 7.5 % Li-6 (natural); scale linearly with enrichment
    constexpr float Li6_natural = 0.075f;
    float enrichment_factor = cfg_.Li6_enrichment / Li6_natural;

    // Flux dependence: TBR saturates above ~1e14 n/m^2/s
    float flux_factor = std::min(neutron_flux / 1e14f, 1.0f);
    float tbr = enrichment_factor * 1.02f * flux_factor;

    return std::min(tbr, 1.5f); // physical upper bound
}