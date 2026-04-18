#pragma once

//
// src/Helium/Helium.h
// Helium ash exhaust and divertor system.
//
//  Helium ash (alpha particles that have slowed down) accumulates in the plasma
//  and dilutes the D-T fuel, reducing fusion power.  The divertor scrapes off
//  the edge plasma into a target plate and pumps it out.
//
//  Also models:
//    - Divertor heat load and tile temperature
//    - Impurity seeding (Ar/Ne) for radiative cooling of the divertor
//    - Pumping speed and throughput
//

#include "ReactorState.h"
#include "SimTime.h"

struct HeliumConfig {
    // He ash accumulation
    float He_confinement_mult   = 5.0f;   // τ_He / τ_E (helium stays ~5x longer)
    float max_He_fraction       = 0.10f;  // above this, burn degrades severely

    // Divertor geometry (ITER-like single-null)
    float divertor_area_m2      = 4.0f;   // target plate area [m^2]
    float max_tile_temp_K       = 2500.0f;// tungsten melting ~3695 K; alarm at 2500 K
    float tile_heat_capacity    = 8e4f;   // J/K (per divertor assembly)

    // Pumping
    float pump_speed_m3s        = 200.0f; // effective pumping speed [m^3/s]
    float max_throughput_Pa_m3s = 200.0f; // maximum throughput [Pa·m^3/s]

    // Impurity seeding (radiative power exhaust)
    float max_seed_rate         = 1.0f;   // normalised 0–1
};

class HeliumSystem {
public:
    explicit HeliumSystem(const HeliumConfig& cfg);

    void update(ReactorState& state, const SimTime& t);

    float heFraction() const { return he_fraction_; }

private:
    void  updateAshAccumulation(ReactorState& state, float dt);
    void  updatePumping        (ReactorState& state, float dt);
    void  updateDivertorThermal(ReactorState& state, float dt);
    float impurityRadiationMW  (float n_e, float Z_eff) const;

    HeliumConfig cfg_;

    float he_fraction_      = 0.0f;    // He ash / total ion number fraction
    float divertor_temp_K_  = 300.0f;
    float seed_rate_        = 0.0f;    // impurity seeding level (0–1)
    float pump_throughput_  = 0.0f;
};