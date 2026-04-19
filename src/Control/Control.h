#pragma once

//
// src/Control/Control.h
// Top-level plasma control system.
// The ONLY module permitted to write to ReactorState setpoint fields.
//
//  Implements:
//    - Plasma current PID (targets sp_plasma_current_MA via V_loop command)
//    - Electron temperature PID (targets sp_electron_temp_keV via aux heating)
//    - Density PID (targets sp_density_m3 via fuelling rate)
//    - Disruption prediction and avoidance
//    - Runaway electron monitoring
//    - SCRAM logic (triggers cmd_scram on any hard limit breach)
//

#include "ReactorState.h"
#include "SimTime.h"
#include "ControlPhysics.h"

struct ControlConfig {
    // PID gains for plasma current [setpoint in MA, output in V_loop volts]
    float Kp_current = 0.5f,  Ki_current = 0.05f,  Kd_current = 0.01f;
    // PID gains for electron temperature [setpoint keV, output aux heat MW]
    float Kp_temp    = 2.0f,  Ki_temp    = 0.1f,   Kd_temp    = 0.5f;
    // PID gains for density [setpoint m^-3, output fuel rate normalised 0-1]
    float Kp_dens    = 5e-21f,Ki_dens    = 5e-22f, Kd_dens    = 1e-21f;

    // Disruption avoidance limits
    float q_limit_low   = 2.0f;   // q < 2 → disruption risk
    float beta_limit    = 0.035f; // β_N > 3.5 % → ideal MHD limit
    float runaway_E_Vm  = 0.1f;   // alarm if E_loop > this [V/m]

    // SCRAM triggers
    float scram_magnet_temp_K  = 15.0f;   // SC quench margin
    float scram_first_wall_K   = 2000.0f;
    float scram_divertor_K     = 2200.0f;
    float scram_beta_hard      = 0.05f;   // hard beta limit
};

class ControlSystem {
public:
    explicit ControlSystem(const ControlConfig& cfg);

    // Main update — reads full ReactorState, writes setpoints and alarms.
    void update(ReactorState& state, const SimTime& t);

private:
    void runCurrentControl    (ReactorState& state, float dt);
    void runTemperatureControl(ReactorState& state, float dt);
    void runDensityControl    (ReactorState& state, float dt);
    void runDisruptionWatch   (ReactorState& state, float dt);
    void runRunawayMonitor    (ReactorState& state);
    void runScramLogic        (ReactorState& state);

    ControlConfig cfg_;

    ControlPhysics::PIDState pid_current_;
    ControlPhysics::PIDState pid_temp_;
    ControlPhysics::PIDState pid_density_;

    float V_loop_command_V_   = 0.0f;  // commanded loop voltage
    float aux_heat_command_MW_= 0.0f;  // commanded auxiliary heating
    bool  scram_latched_      = false;

    // Disruption predictor state
    float q_history_[8]   = {};
    int   q_hist_idx_      = 0;
    float dq_dt_           = 0.0f;    // rate of change of q
};