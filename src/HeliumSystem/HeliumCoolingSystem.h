#pragma once
//
// src/HeliumSystem/HeliumCoolingSystem.h
// Two separate helium circuits:
//
//  A) Reactor structure cooling (NOT cryogenic)
//     - Helium at 8 MPa, ~700 K inlet / ~800 K outlet
//     - Cools first wall support structure and blanket manifolds
//     - 2 circulating pumps
//     - Dumps heat to secondary heat exchanger (→ air coolers or aux steam)
//
//  B) Magnet cryogenic circuit (4.5 K)
//     - Supercritical He at 3 MPa, ~4.5 K
//     - Cools all 18 TF coils + 6 CS coils
//     - 2 cryo-compressors (large, ~35 MW each)
//     - Intermediate 80 K thermal shields
//     - Cryostat: multi-layer vacuum vessel
//

#include "ReactorState.h"
#include "SimTime.h"

struct HeCircuitPump {
    bool  running       = false;
    float speed_frac    = 1.f;
    float flow_kg_s     = 0.f;
    float inlet_temp_K  = 0.f;
    float outlet_temp_K = 0.f;
    float pressure_MPa  = 0.f;
    float power_MW      = 0.f;
    bool  trip          = false;
};

// ─── Reactor He circuit (warm) ────────────────────────────────────────────────
struct ReactorHeCircuit {
    HeCircuitPump pump[2];
    float inlet_temp_K        = 700.f;  // cold He in from HX
    float outlet_temp_K       = 800.f;  // hot He returning from reactor
    float pressure_MPa        = 8.f;
    float heat_removed_MW     = 0.f;
    float aux_HX_duty_MW      = 0.f;    // heat rejected to air coolers
    bool  hi_outlet_temp_alarm = false;
    bool  lo_flow_alarm        = false;
};

// ─── Cryostat ─────────────────────────────────────────────────────────────────
struct Cryostat {
    float vacuum_pressure_Pa  = 1e-4f;  // cryostat insulation vacuum [Pa]
    float thermal_shield_K    = 80.f;   // 80 K nitrogen/He shield temp
    float outer_wall_K        = 293.f;  // room temperature outer vessel
    float heat_leak_W         = 0.f;    // computed static heat in-leak
    bool  vacuum_ok           = true;
    bool  roughing_pump_on    = false;
    bool  turbo_pump_on       = false;
    // Cryostat inner dimensions (roughly ITER-sized)
    float inner_diameter_m    = 11.f;
    float height_m            = 12.f;
    float surface_area_m2     = 400.f;  // cold mass surface area
};

// ─── Magnet cryo circuit (4.5 K) ─────────────────────────────────────────────
struct MagnetCryoCircuit {
    HeCircuitPump cryo_pump[2];         // cryo-compressors/circulators
    float cold_mass_temp_K    = 4.5f;   // TF + CS coil temperature
    float supply_temp_K       = 4.3f;   // He supply from cold box
    float return_temp_K       = 4.7f;   // warm He returning to cold box
    float pressure_MPa        = 0.3f;   // supercritical He pressure
    float total_heat_load_W   = 0.f;    // ohmic + radiation + eddy current
    float cryo_refrigerator_MW= 0.f;    // electrical power to cryo refrigerator
    bool  refrigerator_on     = false;
    bool  lo_temp_alarm       = false;  // <4.0 K → risk of heater overshoot
    bool  hi_temp_alarm       = false;  // >5.5 K → quench risk
    bool  lo_pressure_alarm   = false;
};

struct HeliumSystemState {
    ReactorHeCircuit reactor_circuit;
    MagnetCryoCircuit magnet_circuit;
    Cryostat         cryostat;
};

class HeliumCoolingSystem {
public:
    HeliumCoolingSystem();

    void update(ReactorState& state, const SimTime& t);

    const HeliumSystemState& heState() const { return s_; }
    HeliumSystemState&       heState()       { return s_; }

    // Operator commands
    void startReactorCooling();
    void stopReactorCooling();
    void startCryoplant();
    void stopCryoplant();
    void startCryostatPumping();

private:
    void updateReactorCircuit(float dt, float fusion_power_MW,
                               float first_wall_temp_K);
    void updateCryostat(float dt);
    void updateMagnetCircuit(float dt, float magnet_temp_K,
                              float stored_energy_GJ);

    HeliumSystemState s_;
};