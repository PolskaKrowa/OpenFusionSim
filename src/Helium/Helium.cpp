//
// src/Helium/Helium.cpp
//

#include "Helium.h"
#include <algorithm>
#include <cmath>

HeliumSystem::HeliumSystem(const HeliumConfig& cfg)
    : cfg_(cfg), divertor_temp_K_(300.0f)
{}

void HeliumSystem::reset()
{
    he_fraction_      = 0.0f;
    divertor_temp_K_  = 300.0f;
    seed_rate_        = 0.0f;
    pump_throughput_  = 0.0f;
}

void HeliumSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    updateAshAccumulation(state, dt);
    updatePumping        (state, dt);
    updateDivertorThermal(state, dt);

    // ── helium_fraction ownership ───────────────────────────────────────────
    //  The Fuel module's accountExhaust() writes state.helium_fraction based
    //  on a rigorous particle-count bookkeeping (N_He_ash / N_total).  We
    //  previously OVERWROTE that here with our own ODE-based he_fraction_,
    //  which was a separate, inconsistent model of the same quantity.  The
    //  PlasmaCoreBridge read whichever value was written last (Helium's, since
    //  it runs after Fuel in the main loop), ignoring Fuel's model entirely.
    //
    //  Now: Fuel owns state.helium_fraction.  Helium's he_fraction_ is still
    //  tracked internally (for the divertor pumping model) but is NOT written
    //  back to ReactorState.  This eliminates the double-write race.
    //
    //  (We DO still write pump_throughput_Pa and divertor_temp_K — those are
    //  Helium's own fields that no other module touches.)
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

    // He confinement time = τ_He_mult × τ_E.
    //  Use the actual energy confinement time from the PlasmaCoreBridge
    //  (state.tau_E_s, computed by IPB98(y,2)) rather than the old hardcoded
    //  3.0 s.  At full ITER power, τ_E is ~3.7 s; during ramp-up it's much
    //  lower (~1 s), which means He exhausts faster early in the discharge.
    float tau_E = (state.tau_E_s > 0.1f) ? state.tau_E_s : 1.0f;
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
    // Guard: if plasma state is non-finite (shouldn't happen after the
    // power-balance guards, but belt-and-suspenders), skip the update
    // rather than producing NaN in edge_pressure_Pa / pump_throughput_.
    if (!std::isfinite(state.plasma_density_m3) ||
        !std::isfinite(state.plasma_temp_keV)) {
        return;
    }

    // Divertor pumping removes He + D/T from the scrape-off layer.
    //  Throughput = pump_speed [m³/s] × edge_pressure [Pa]  →  [Pa·m³/s] = [W]
    float edge_pressure_Pa = state.plasma_density_m3 * 1.38e-23f
                           * (state.plasma_temp_keV * 1e3f * 1.602e-19f / 1.38e-23f)
                           * 1e-4f; // very rough edge estimate

    pump_throughput_ = std::min(cfg_.pump_speed_m3s * edge_pressure_Pa,
                                cfg_.max_throughput_Pa_m3s);

    // He removal rate (fraction per second).
    //  The volumetric pumping rate is pump_speed [m³/s].  The fraction of
    //  the plasma volume pumped per second is pump_speed / V_plasma.  He is
    //  removed at that rate × its fraction.  This gives units of 1/s, which
    //  is correct for a fraction-per-second decay.
    //
    //  The old code used `pump_throughput_ / (n_e × V)` which has units of
    //  (J/s) / particles = J·s⁻¹·particle⁻¹ — not a fraction-per-second
    //  rate.  Numerically it gave ~2.4e-22 (effectively zero), so He pumping
    //  did nothing.  The new formula gives a realistic ~0.24/s at full pump
    //  speed (200 m³/s / 840 m³), so He exhausts with τ ≈ 4 s.
    constexpr float V_plasma = 840.0f;  // m³
    float He_removal_rate = cfg_.pump_speed_m3s / V_plasma;  // 1/s
    he_fraction_ -= He_removal_rate * he_fraction_ * dt;
    he_fraction_  = std::max(0.0f, he_fraction_);
}

void HeliumSystem::updateDivertorThermal(ReactorState& state, float dt)
{
    // Power to divertor ≈ P_rad + power not captured by blanket
    // Realistic ITER-like split: ~5 % of total fusion power hits the divertor
    // (alpha + neutral beam shine-through + charge-exchange), and ~50 % of
    // the *radiated* power reaches the divertor (the other 50 % is
    // deposited on the first wall).
    //
    // IMPORTANT: do NOT overwrite state.radiated_power_MW here.  That field
    // is owned by PlasmaCoreBridge (brem + sync + line radiation, computed
    // from local n_e/T_e/B per cell).  Overwriting it with `P_seeded + 3 %
    // of P_fus` discards the bridge's value and creates a runaway feedback
    // loop: high radiation → hot divertor → high seed_rate → "radiated" goes
    // up → hotter divertor → ... which was the second half of the
    // spurious-overtemperature-SCRAM bug.
    float P_div = state.fusion_power_MW * 0.05f + state.radiated_power_MW * 0.5f;

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
    // Note: state.radiated_power_MW is NOT touched here.  The bridge owns it.
}