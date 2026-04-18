//
// src/Fuel/Fuel.cpp
//

#include "Fuel.h"
#include <algorithm>
#include <cmath>

FuelSystem::FuelSystem(const FuelConfig& cfg)
    : cfg_(cfg)
    , D_inventory_g_(cfg.initial_D_g)
    , T_inventory_g_(cfg.initial_T_g)
{}

void FuelSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // Warn if tritium is low (game mechanic: T takes time to breed/import)
    state.alarm_low_tritium = (T_inventory_g_ < 10.0f);

    runGasPuffing    (state, dt);
    runPelletInjector(state, dt);
    accountBurnup    (state, dt);
    accountRecycling (state, dt);

    // Write inventory to state for UI / other modules
    state.fuel_D_inventory_g = D_inventory_g_;
    state.fuel_T_inventory_g = T_inventory_g_;
}

void FuelSystem::runGasPuffing(ReactorState& state, float dt)
{
    if (state.plasma_status == PlasmaStatus::Cold) return;

    // Gas puff rate driven by Control setpoint and density error
    float rate = state.sp_fuel_rate * cfg_.max_gas_rate_Pa_m3s;
    rate = std::clamp(rate, 0.0f, cfg_.max_gas_rate_Pa_m3s);

    // Convert Pa·m^3/s → grams/s: n = PV/kT at ~room temp, then × mass
    constexpr float T_room = 293.0f, k_B = 1.38e-23f;
    constexpr float m_D = 2.0f * 1.67e-27f * 1000.0f; // g per molecule
    constexpr float m_T = 3.0f * 1.67e-27f * 1000.0f;

    float n_per_s = rate / (k_B * T_room);
    float D_rate_g = n_per_s * m_D * state.D_T_ratio / (1.0f + state.D_T_ratio);
    float T_rate_g = n_per_s * m_T / (1.0f + state.D_T_ratio);

    float D_used = D_rate_g * dt;
    float T_used = T_rate_g * dt;

    // Clamp to available inventory
    D_used = std::min(D_used, D_inventory_g_);
    T_used = std::min(T_used, T_inventory_g_);
    D_inventory_g_ -= D_used;
    T_inventory_g_ -= T_used;

    state.fuel_injection_rate = rate;
}

void FuelSystem::runPelletInjector(ReactorState& state, float dt)
{
    if (state.plasma_status == PlasmaStatus::Cold) return;
    if (state.pellet_frequency_Hz <= 0.0f) return;

    float freq = std::min(state.pellet_frequency_Hz, cfg_.max_pellet_freq_Hz);
    float interval = 1.0f / freq;

    pellet_timer_s_ -= dt;
    if (pellet_timer_s_ <= 0.0f) {
        pellet_timer_s_ = interval;

        // Fire one pellet (50% D / 50% T by mass)
        float pellet_D = cfg_.pellet_mass_g * 0.5f;
        float pellet_T = cfg_.pellet_mass_g * 0.5f;

        // Deduct from inventory
        pellet_D = std::min(pellet_D, D_inventory_g_);
        pellet_T = std::min(pellet_T, T_inventory_g_);
        D_inventory_g_ -= pellet_D;
        T_inventory_g_ -= pellet_T;

        // Deposited mass increases plasma density (simple model)
        constexpr float V_plasma = 840.0f; // m^3
        constexpr float m_D_kg   = 2.0f * 1.67e-27f;
        float N_deposited = (pellet_D * 1e-3f / m_D_kg) * cfg_.pellet_efficiency;
        state.plasma_density_m3 += N_deposited / V_plasma;
    }
}

void FuelSystem::accountBurnup(ReactorState& state, float dt)
{
    // D-T fusion: 1 D + 1 T → He-4 + n per reaction
    float energy_MJ  = state.fusion_power_MW * dt;
    float D_burned_g = cfg_.D_consumption_g_per_MWs * energy_MJ;
    float T_burned_g = cfg_.T_consumption_g_per_MWs * energy_MJ;

    D_burned_g  = std::min(D_burned_g, D_inventory_g_);
    T_burned_g  = std::min(T_burned_g, T_inventory_g_);
    D_inventory_g_ -= D_burned_g;
    T_inventory_g_ -= T_burned_g;
    D_consumed_g_  += D_burned_g;
    T_consumed_g_  += T_burned_g;
}

void FuelSystem::accountRecycling(ReactorState& state, float dt)
{
    // Fuel that hits the wall partially returns to the plasma edge.
    // Net effect: reduce apparent consumption.  Wall saturation not modelled.
    float recycled_D = cfg_.recycling_coeff
                     * state.fuel_injection_rate * 0.5f * dt * 1e-6f; // crude
    D_inventory_g_ += recycled_D * 0.01f; // very small — mostly stays at wall
    (void)state;
}

void FuelSystem::resupplyDeuterium(float grams)
{
    D_inventory_g_ = std::min(D_inventory_g_ + grams, cfg_.max_D_inventory_g);
}

void FuelSystem::resupplyTritium(float grams)
{
    T_inventory_g_ = std::min(T_inventory_g_ + grams, cfg_.max_T_inventory_g);
}