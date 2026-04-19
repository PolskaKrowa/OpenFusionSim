//
// src/Control/Control.cpp
//

#include "Control.h"
#include "ControlPhysics.h"
#include <cmath>
#include <algorithm>
#include <numeric>

ControlSystem::ControlSystem(const ControlConfig& cfg)
    : cfg_(cfg)
{}

void ControlSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // SCRAM overrides everything — check first
    runScramLogic(state);
    if (scram_latched_) {
        state.cmd_scram          = true;
        state.sp_plasma_current_MA = 0.0f;
        state.sp_fuel_rate         = 0.0f;
        return;
    }

    runCurrentControl    (state, dt);
    runTemperatureControl(state, dt);
    runDensityControl    (state, dt);
    runDisruptionWatch   (state, dt);
    runRunawayMonitor    (state);
}

// ─── Plasma Current PID ────────────────────────────────────────────────────────
//  Output: loop voltage command → used by PlasmaCore to drive Ip
void ControlSystem::runCurrentControl(ReactorState& state, float dt)
{
    float measured = state.plasma_current_MA;
    float setpoint = state.sp_plasma_current_MA;

    // Disruption avoidance: clamp setpoint if q is too low
    if (state.q_safety < cfg_.q_limit_low * 1.1f)
        setpoint = std::min(setpoint, state.plasma_current_MA * 0.95f);

    V_loop_command_V_ = ControlPhysics::feedbackControlPID(
        setpoint, measured,
        cfg_.Kp_current, cfg_.Ki_current, cfg_.Kd_current,
        dt, pid_current_,
        /*output_min=*/ -5.0f,   // V
        /*output_max=*/  20.0f,  // V  (typical inductive startup: ~5–10 V)
        /*windup=*/      100.0f
    );

    // Compute bootstrap fraction and net required loop voltage
    float pressure_Pa = state.plasma_density_m3
                      * state.plasma_temp_keV * 1e3f * 1.602e-19f * 2.0f; // n(Ti+Te)
    float epsilon     = 2.0f / 6.2f; // a/R for ITER-class

    float I_bootstrap = ControlPhysics::bootstrapCurrent(
        pressure_Pa, state.B_toroidal_T,
        epsilon, state.q_safety,
        state.plasma_current_MA * 1e6f);

    float bootstrap_fraction = (state.plasma_current_MA > 0.1f)
        ? I_bootstrap / (state.plasma_current_MA * 1e6f) : 0.0f;

    // Reduce required Ohmic current by bootstrap fraction
    // (informational — PlasmaCore uses V_loop_command directly)
    (void)bootstrap_fraction;
}

// ─── Temperature PID ──────────────────────────────────────────────────────────
//  Output: auxiliary heating power command [MW] (NBI + ICRH)
void ControlSystem::runTemperatureControl(ReactorState& state, float dt)
{
    aux_heat_command_MW_ = ControlPhysics::feedbackControlPID(
        state.sp_electron_temp_keV,
        state.electron_temp_keV,
        cfg_.Kp_temp, cfg_.Ki_temp, cfg_.Kd_temp,
        dt, pid_temp_,
        /*output_min=*/ 0.0f,    // can't have negative heating
        /*output_max=*/ 100.0f,  // 100 MW max aux heating (NBI + ICRH)
        /*windup=*/     1000.0f
    );
    // PlasmaCore reads aux_heat from ReactorState — write it via a shared field
    // (extend ReactorState with sp_aux_heat_MW if needed)
}

// ─── Density PID ──────────────────────────────────────────────────────────────
//  Output: fuel injection rate setpoint [0–1 normalised]
void ControlSystem::runDensityControl(ReactorState& state, float dt)
{
    float rate = ControlPhysics::feedbackControlPID(
        state.sp_density_m3,
        state.plasma_density_m3,
        cfg_.Kp_dens, cfg_.Ki_dens, cfg_.Kd_dens,
        dt, pid_density_,
        /*output_min=*/ 0.0f,
        /*output_max=*/ 1.0f,
        /*windup=*/     10.0f
    );
    state.sp_fuel_rate = rate;
}

// ─── Disruption Watch ─────────────────────────────────────────────────────────
void ControlSystem::runDisruptionWatch(ReactorState& state, float dt)
{
    // Track q safety factor history
    q_history_[q_hist_idx_ % 8] = state.q_safety;
    q_hist_idx_++;

    // Estimate dq/dt by linear fit over history window
    if (q_hist_idx_ >= 8) {
        float q_old = q_history_[(q_hist_idx_ + 1) % 8];
        float q_new = q_history_[(q_hist_idx_    ) % 8];
        dq_dt_      = (q_new - q_old) / (8.0f * dt);
    }

    // Predict q in 100 ms — if projected below limit, pre-emptively reduce Ip
    float q_predicted = state.q_safety + dq_dt_ * 0.1f;
    if (q_predicted < cfg_.q_limit_low) {
        // Ramp down plasma current setpoint gently
        state.sp_plasma_current_MA = std::max(
            state.sp_plasma_current_MA - 0.1f * dt,
            0.5f  // never command below 0.5 MA (then do a controlled shutdown)
        );
        state.disruption_flag  = true;
        state.alarm_disruption = true;
    }

    // Beta limit check (normalised beta β_N = β * a * B / (μ₀ * Ip))
    const float a_m = 2.0f;
    float beta_N = state.beta * a_m * state.B_toroidal_T
                 / (1.257e-6f * state.plasma_current_MA * 1e6f + 1.0f);
    if (beta_N > cfg_.beta_limit) {
        // Reduce heating / density setpoint to drop pressure
        state.sp_electron_temp_keV *= 0.99f;
    }
}

// ─── Runaway Electron Monitor ─────────────────────────────────────────────────
void ControlSystem::runRunawayMonitor(ReactorState& state)
{
    // Loop electric field from Faraday's law: E_loop = V_loop / (2π * R)
    const float R = 6.2f; // major radius [m]
    float E_loop = fabsf(V_loop_command_V_) / (2.0f * 3.14159f * R);

    float T_eV = state.electron_temp_keV * 1000.0f;
    float coulomb_log = 15.0f; // typical for DT plasma at 10-20 keV

    ControlPhysics::RunawayResult re = ControlPhysics::runawayElectronThreshold(
        E_loop, state.plasma_density_m3,
        1.5f,   // Z_eff estimate
        T_eV, coulomb_log);

    if (re.runaway_risk && re.growth_rate_s > 1e6f) {
        // Significant runaway: increase density (higher collisionality) to suppress
        state.sp_density_m3 = std::min(state.sp_density_m3 * 1.01f, 1.5e20f);
    }
}

// ─── SCRAM Logic ──────────────────────────────────────────────────────────────
//  Latching — once triggered, requires an explicit reset (not modelled here)
void ControlSystem::runScramLogic(ReactorState& state)
{
    if (scram_latched_) return;

    bool trigger = false;

    // Hard limits — any one trips the SCRAM
    trigger |= (state.magnet_temp_K     > cfg_.scram_magnet_temp_K);
    trigger |= (state.quench_detected);
    trigger |= (state.first_wall_temp_K > cfg_.scram_first_wall_K);
    trigger |= (state.divertor_temp_K   > cfg_.scram_divertor_K);
    trigger |= (state.beta              > cfg_.scram_beta_hard);
    trigger |= (state.thermal_runaway);
    trigger |= (state.alarm_loss_of_coolant);

    // Disruption that has already crossed q < 2 (not just predicted)
    trigger |= (state.q_safety < 1.5f && state.plasma_current_MA > 1.0f);

    if (trigger) {
        scram_latched_  = true;
        state.cmd_scram = true;
        state.mode      = ReactorMode::Emergency;
    }
}