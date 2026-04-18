#pragma once

//
// src/Fuel/Fuel.h
// D-T fuel management system.
//
//  Tracks:
//    - Gas puffing (continuous, low rate — edge fuelling)
//    - Pellet injection (high-frequency, deep fuelling — main density control)
//    - D-T inventory (separate stores; tritium is scarce and radioactive)
//    - Fuel burnup and recycling from the wall
//

#include "ReactorState.h"
#include "SimTime.h"
#include <deque>

struct FuelConfig {
    // Inventory limits
    float max_D_inventory_g   = 1000.0f;  // [g]
    float max_T_inventory_g   =  200.0f;  // [g] — limited by breeding & regulation
    float initial_D_g         =  800.0f;
    float initial_T_g         =  150.0f;

    // Gas puffing
    float max_gas_rate_Pa_m3s = 10.0f;    // max fuelling rate [Pa·m^3/s]
    float fuelling_efficiency  = 0.25f;   // fraction that enters plasma core

    // Pellet injector
    float max_pellet_freq_Hz  = 10.0f;    // max injection frequency [Hz]
    float pellet_mass_g       = 0.05e-3f; // mass per pellet [g]  (50 μg typical)
    float pellet_efficiency   = 0.80f;    // fraction deposited in core

    // Burnup
    float D_consumption_g_per_MWs = 3.7e-8f;  // g of D burned per MJ of fusion
    float T_consumption_g_per_MWs = 5.6e-8f;  // g of T burned per MJ of fusion

    // Wall recycling coefficient
    float recycling_coeff = 0.85f;  // fraction of particles hitting wall that return
};

class FuelSystem {
public:
    explicit FuelSystem(const FuelConfig& cfg);

    void update(ReactorState& state, const SimTime& t);

    float deuteriumInventory() const  { return D_inventory_g_; }
    float tritiumInventory()   const  { return T_inventory_g_; }

    // External resupply (called by game events / logistics system)
    void  resupplyDeuterium(float grams);
    void  resupplyTritium  (float grams);

private:
    void  runGasPuffing    (ReactorState& state, float dt);
    void  runPelletInjector(ReactorState& state, float dt);
    void  accountBurnup    (ReactorState& state, float dt);
    void  accountRecycling (ReactorState& state, float dt);

    FuelConfig cfg_;
    float D_inventory_g_;
    float T_inventory_g_;
    float pellet_timer_s_ = 0.0f;  // time until next pellet fires
    float D_consumed_g_   = 0.0f;  // cumulative
    float T_consumed_g_   = 0.0f;
};