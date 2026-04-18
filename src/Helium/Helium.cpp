//
// src/Helium/Helium.cpp
//

#include "Helium.h"
#include <algorithm>
#include <cmath>

HeliumSystem::HeliumSystem(const HeliumConfig& cfg)
    : cfg_(cfg), divertor_temp_K_(300.0f)
{}

void HeliumSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    updateAshAccumulation(state, dt);
    updatePumping        (state, dt);
    updateDivertorThermal(state, dt);

    state.helium_fraction   = he_fraction_;
    state.pump_throughput_Pa = pump_throughput_;
    state.divertor_temp_K    = divertor_temp_K_;
    state.divertor_overtemp  = (divertor_temp_K_ > cfg_.max_tile_temp_K);
    state.alarm_overtemp    |= state.divertor_overtemp;
}

void HeliumSystem::updateAshAccumulation(ReactorState& state, float dt)
{
    if (state.plasma_status != PlasmaStatus::Burning) {
        // No burn — He drains naturally
        he_fraction_ = std::max(0.0f, he_fraction_ - 0.01f * dt);
        return;
    }

    // He production rate proportional to fusion power
    // Each D-T reaction produces one alpha that thermalises in ~τ_E
    // He source rate (fraction per second) ≈ P_fusion / (E_alpha * N_plasma)
    constexpr float E_alpha_J = 3.52e6f * 1.602e-19f;
    constexpr float V_plasma  = 840.0f;
    float N_plasma = state.plasma_density_m3 * V_plasma;
    if (N_plasma < 1.0f) return;

    float alphas_per_s = state.alpha_power_MW * 1e6f / E_alpha_J;
    float source_frac  = alphas_per_s / N_plasma; // fraction/s added

    // He confinement time = τ_He_mult * τ_E (τ_E ~ 3 s for ITER)
    float tau_E = 3.0f; // energy confinement time [s] — simplification
    float tau_He = cfg_.He_confinement_mult * tau_E;
    float sink_frac = he_fraction_ / tau_He; // exhaust rate

    he_fraction_ += (source_frac - sink_frac) * dt;
    he_fraction_  = std::clamp(he_fraction_, 0.0f, cfg_.max_He_fraction * 2.0f);

    // Feed back: He dilution reduces effective D-T density and thus fusion power
    // PlasmaCore uses this fraction to weight its density
    // (ReactorState propagates it; PlasmaCoreBridge reads it)
}

void HeliumSystem::updatePumping(ReactorState& state, float dt)
{
    // Divertor pumping removes He + D/T from the scrape-off layer
    // Throughput ~ pump_speed * edge_pressure
    float edge_pressure_Pa = state.plasma_density_m3 * 1.38e-23f
                           * (state.plasma_temp_keV * 1e3f * 1.602e-19f / 1.38e-23f)
                           * 1e-4f; // very rough edge estimate

    pump_throughput_ = std::min(cfg_.pump_speed_m3s * edge_pressure_Pa,
                                cfg_.max_throughput_Pa_m3s);

    // He removed per second
    float He_removal_rate = pump_throughput_ / (state.plasma_density_m3 * 840.0f + 1.0f);
    he_fraction_ -= He_removal_rate * he_fraction_ * dt;
    he_fraction_  = std::max(0.0f, he_fraction_);
    (void)dt;
    (void)state;
}

void HeliumSystem::updateDivertorThermal(ReactorState& state, float dt)
{
    // Power to divertor ≈ P_rad + power not captured by blanket
    // Simplified: ~15 % of total fusion power goes to divertor
    float P_div = state.fusion_power_MW * 0.15f + state.radiated_power_MW * 0.5f;

    // Impurity seeding radiates some of this before hitting the tile
    // Real control: inject Ar/Ne to keep tile below limit
    if (divertor_temp_K_ > cfg_.max_tile_temp_K * 0.8f)
        seed_rate_ = std::min(seed_rate_ + 0.1f * dt, cfg_.max_seed_rate);
    else
        seed_rate_ = std::max(seed_rate_ - 0.05f * dt, 0.0f);

    float P_seeded_MW = seed_rate_ * P_div * 0.6f; // seed radiates up to 60 %
    float P_tile_MW   = P_div - P_seeded_MW;

    // Tile temperature: Q_dot = P_tile / A → dT/dt = P_tile / C_heat
    divertor_temp_K_ += (P_tile_MW * 1e6f / cfg_.tile_heat_capacity) * dt;

    // Water-cooled tiles: active cooling removes heat proportional to ΔT
    float T_coolant  = 500.0f; // K
    float UA_diver   = 2e6f;   // [W/K] heat transfer coefficient × area
    float Q_removed  = UA_diver * (divertor_temp_K_ - T_coolant) * dt
                     / cfg_.tile_heat_capacity;
    divertor_temp_K_ -= Q_removed;
    divertor_temp_K_  = std::max(divertor_temp_K_, T_coolant);

    state.divertor_power_MW = P_tile_MW;
    state.radiated_power_MW = P_seeded_MW + state.fusion_power_MW * 0.03f; // Bremss
}