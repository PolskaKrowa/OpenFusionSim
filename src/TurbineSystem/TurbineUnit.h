#pragma once
//
// src/TurbineSystem/TurbineUnit.h
// State and physics for one 1200 MWe turbine unit.
//
//  Each unit contains:
//    - Steam generator (SG) with 2 feedwater pumps
//    - Main Steam Isolation Valve (MSIV)
//    - Turbine bypass valve
//    - Steam relief valve
//    - 4 feedwater preheaters (togglable)
//    - Condenser (4 CAR pumps, SJAE, variable condensate pumps)
//    - Hotwell (makeup + drain)
//    - Generator with grid synchronisation logic
//

#include "SimTime.h"
#include <cmath>
#include <cstdint>
#include <algorithm>

// ─── Turbine unit operating state ─────────────────────────────────────────────
enum class TurbineState : uint8_t {
    Offline,       // cold, no steam, MSIV closed
    RollingUp,     // MSIV cracking, turbine spinning up
    Synchronizing, // ~3000 RPM, matching grid phase/frequency
    Online,        // breaker closed, exporting power
    Runback,       // active load reduction in progress
    Tripping,      // MSIV closing, coasting to stop
    Tripped,       // shutdown, awaiting reset
};

// ─── Feedwater pump ──────────────────────────────────────────────────────────
struct FeedwaterPump {
    bool  running          = false;
    float speed_frac       = 0.f;   // 0–1, normally 1.0 when running
    float flow_kg_s        = 0.f;   // computed
    float discharge_MPa    = 0.f;   // computed discharge pressure
    float power_MW         = 0.f;   // pump electrical consumption
    bool  trip             = false; // motor protection trip
};

// ─── Feedwater preheater ──────────────────────────────────────────────────────
struct Preheater {
    bool  enabled              = true;  // operator can toggle
    float extraction_frac      = 0.f;  // fraction of steam extracted from turbine
    float fw_outlet_temp_K     = 0.f;  // feedwater temperature leaving this stage
    float drain_level_m        = 0.5f; // condensate level in shell
    float heat_transferred_MW  = 0.f;
};

// ─── Condenser ────────────────────────────────────────────────────────────────
struct CondenserState {
    float pressure_kPa         = 6.0f;   // target ~5–6 kPa
    float temp_K               = 309.f;  // saturation temp at condenser pressure
    float air_fraction         = 0.001f; // non-condensable gas mole fraction
    bool  car_pump[4]          = {};     // Condenser Air Removal pumps
    bool  sjae_running         = false;  // Steam Jet Air Ejector
    float sjae_steam_kg_s      = 0.f;   // steam consumed by SJAE
    float condensate_pump_speed= 0.f;   // 0–1 variable condensate pump
    float condensate_flow_kg_s = 0.f;
    float cooling_water_temp_K = 295.f;  // cooling water (sea/tower) inlet
    float CW_flow_kg_s         = 8000.f; // cooling water flow
};

// ─── Hotwell ─────────────────────────────────────────────────────────────────
struct Hotwell {
    float level_m              = 1.5f;   // normal: 1.5 m (range 0.5–3.0 m)
    float area_m2              = 40.f;   // cross-section area
    bool  makeup_valve         = false;  // adds demineralised water
    bool  drain_valve          = false;  // removes excess
    float makeup_flow_kg_s     = 5.f;   // flow when valve open
    float drain_flow_kg_s      = 8.f;
    bool  lo_level_alarm       = false;
    bool  hi_level_alarm       = false;
};

// ─── Generator ────────────────────────────────────────────────────────────────
struct GeneratorState {
    float terminal_voltage_kV  = 24.f;   // stator terminal voltage
    float frequency_Hz         = 0.f;    // rotor electrical frequency
    float phase_rad            = 0.f;    // rotor phase angle
    float power_MW             = 0.f;    // active power output
    float reactive_MVAR        = 0.f;    // reactive power
    float field_current_A      = 800.f;  // excitation (controls voltage/reactive)
    bool  breaker_closed       = false;  // grid connection breaker
    bool  exciter_on           = false;  // AVR / excitation system
    bool  overspeed_trip       = false;  // >110% RPM
    bool  diff_trip            = false;  // differential protection
};

// ─── Full turbine unit ────────────────────────────────────────────────────────
struct TurbineUnit {
    int id = 0; // 1–4

    TurbineState state = TurbineState::Offline;

    // ── Steam generator ──────────────────────────────────────────────────────
    float sg_pressure_MPa       = 0.1f;  // steam drum pressure
    float sg_steam_temp_K       = 310.f; // saturated steam temperature
    float sg_level_m            = 7.5f;  // water level in drum (normal 7.5 m)
    float sg_steam_flow_kg_s    = 0.f;   // steam production rate
    float sg_heat_input_MW      = 0.f;   // heat from molten salt (written by MoltenSalt)
    FeedwaterPump fw_pump[2];

    // ── MSIV & steam path ────────────────────────────────────────────────────
    bool  msiv_open             = false;
    float msiv_position         = 0.f;   // 0=closed, 1=full open (ramps)
    float msiv_setpoint         = 1.f;   // operator demand 0–1 (throttle)
    bool  msiv_trip_latch       = false; // latching trip; reset required

    float bypass_valve_pos      = 0.f;   // steam bypass to condenser, 0=closed
    bool  relief_valve_open     = false;
    float relief_setpoint_MPa   = 20.f;  // opens when sg_pressure > this

    // ── Preheaters ────────────────────────────────────────────────────────────
    Preheater preheater[4];
    float fw_temp_after_ph_K    = 400.f; // feedwater temp entering SG after all PHs

    // ── Turbine shaft ─────────────────────────────────────────────────────────
    float rpm                   = 0.f;   // rotor speed [RPM], target 3000
    float shaft_power_MW        = 0.f;   // mechanical shaft power
    float steam_flow_to_turbine = 0.f;   // kg/s through HP stages
    float governor_demand       = 0.f;   // governor valve 0–1
    bool  overspeed_trip        = false;

    // ── Condenser ─────────────────────────────────────────────────────────────
    CondenserState condenser;

    // ── Hotwell ───────────────────────────────────────────────────────────────
    Hotwell hotwell;

    // ── Generator ─────────────────────────────────────────────────────────────
    GeneratorState generator;

    // ── Per-unit alarms ───────────────────────────────────────────────────────
    bool  alarm_hi_sg_pressure  = false;
    bool  alarm_lo_sg_level     = false;
    bool  alarm_lo_condenser_vac= false; // condenser pressure high
    bool  alarm_turb_trip       = false;
};

// ─── TurbineUnitController — physics + automatic logic for one unit ───────────
class TurbineUnitController {
public:
    TurbineUnit s; // all mutable state

    explicit TurbineUnitController(int id) { s.id = id; initPreheaters(); }

    // Call each simulation tick
    void update(float dt, float grid_frequency_Hz);

    // Operator commands
    void cmdStart();                   // initiate roll-up sequence
    void cmdStop();                    // controlled shutdown
    void cmdTrip();                    // immediate trip (MSIV close)
    void cmdReset();                   // reset trip latches, return to Offline
    void cmdCloseBreakerRequest();     // attempt grid synchronisation
    void cmdOpenBreaker();             // disconnect from grid

    float netPowerMW() const { return s.generator.breaker_closed ? s.generator.power_MW : 0.f; }

private:
    void updateSteamGenerator(float dt);
    void updateMSIVAndSteamPath(float dt);
    void updatePreheaters(float dt);
    void updateTurbineShaft(float dt);
    void updateCondenser(float dt);
    void updateHotwell(float dt);
    void updateGenerator(float dt, float grid_freq_Hz);
    void checkAlarms();
    void runStateMachine(float dt, float grid_freq_Hz);
    void initPreheaters();

    bool  breaker_close_requested_ = false;
    float sync_timer_s_            = 0.f;  // time in Synchronizing state
    float phase_error_integrated_  = 0.f;
};