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
    //  Reduced Kd from 0.5 to 0.1 — the high derivative gain was causing
    //  oscillation when T_e changed rapidly (e.g. on H&CD enable/disable).
    float Kp_temp    = 2.0f,  Ki_temp    = 0.1f,   Kd_temp    = 0.1f;
    // PID gains for density [setpoint m^-3, output fuel rate normalised 0-1]
    float Kp_dens    = 5e-21f,Ki_dens    = 5e-22f, Kd_dens    = 1e-21f;

    // Disruption avoidance limits
    float q_limit_low   = 2.0f;   // q < 2 → disruption risk
    float beta_limit    = 2.5f;   // β_N > 2.5 → Troyon ideal MHD limit
    float runaway_E_Vm  = 0.1f;   // alarm if E_loop > this [V/m]

    // SCRAM triggers
    float scram_magnet_temp_K  = 15.0f;   // SC quench margin
    float scram_first_wall_K   = 2000.0f;
    float scram_divertor_K     = 2200.0f;
    float scram_beta_hard      = 0.06f;   // hard beta limit (6% — ITER design is ~3%)

    // Controlled shutdown (Rampdown mode) ramp rates
    //  These are MUCH slower than the SCRAM zero-everything path — the
    //  whole point of a controlled shutdown is to keep the plasma stable
    //  while you ramp it down.  These rates keep q_95 and β_N inside their
    //  safe envelopes throughout.
    float rampdown_dIp_dt   = 0.3f;     // MA/s  (vs 0.5 MA/s during startup)
    float rampdown_dTe_dt   = 0.5f;     // keV/s
    float rampdown_dne_dt   = 1.0e19f;  // m⁻³/s
    float rampdown_dPaux_dt = 5.0f;     // MW/s  (back out the H&CD systems gracefully)
};

class ControlSystem {
public:
    explicit ControlSystem(const ControlConfig& cfg);

    // Main update — reads full ReactorState, writes setpoints and alarms.
    void update(ReactorState& state, const SimTime& t);

    // ── Reset the latched SCRAM condition ────────────────────────────────────
    //
    //  Once SCRAM is triggered, scram_latched_ stays true forever (the SCRAM
    //  is "latched" — operators must explicitly reset it).  Call this method
    //  when the operator presses RESET — without it, the next update() call
    //  sees the stale latch and re-triggers cmd_scram, trapping the sim in
    //  a permanent SCRAM loop.
    //
    //  Also clears the PID integrators and q-history so the controller
    //  restarts cleanly after a cold restart.
    //
    void resetScramLatch();

private:
    void runCurrentControl    (ReactorState& state, float dt);
    void runTemperatureControl(ReactorState& state, float dt);
    void runDensityControl    (ReactorState& state, float dt);
    void runDisruptionWatch   (ReactorState& state, float dt);
    void runRunawayMonitor    (ReactorState& state, float dt);
    void runScramLogic        (ReactorState& state);
    void runRampdown          (ReactorState& state, float dt);
    void runShapeControl      (ReactorState& state, float dt);

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