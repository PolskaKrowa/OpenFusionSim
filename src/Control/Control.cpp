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

void ControlSystem::resetScramLatch()
{
    scram_latched_  = false;
    // Clear PID integrators so a cold restart doesn't immediately re-trip
    // due to a large integrated error from before the SCRAM.
    pid_current_  = ControlPhysics::PIDState{};
    pid_temp_     = ControlPhysics::PIDState{};
    pid_density_  = ControlPhysics::PIDState{};
    // Reset q-history so the disruption predictor doesn't fire on stale data
    for (int i = 0; i < 8; i++) q_history_[i] = 0.0f;
    q_hist_idx_  = 0;
    dq_dt_       = 0.0f;
    // Reset actuator commands
    V_loop_command_V_    = 0.0f;
    aux_heat_command_MW_ = 0.0f;
}

void ControlSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // SCRAM overrides everything — check first
    runScramLogic(state);
    if (scram_latched_) {
        state.cmd_scram            = true;
        // Zero ALL plasma setpoints, not just current and fuel.  Without
        // zeroing T_e and density setpoints, the PlasmaCoreBridge keeps
        // ramping T_e toward 20 keV and ne toward 1e20 even as I_p → 0.
        // That leaves the plasma in a numerically degenerate state
        // (high T, high n, zero I) where tau_E → 0 and several formulas
        // in the power balance divide by zero, producing NaN that
        // propagates to every ReactorState field.
        state.sp_plasma_current_MA = 0.0f;
        state.sp_fuel_rate         = 0.0f;
        state.sp_electron_temp_keV = 0.0f;
        state.sp_density_m3        = 0.0f;
        // Also zero the H&CD setpoints so the PlasmaCoreBridge sees a
        // ramp-down of P_aux instead of holding it at the last commanded
        // value (which was a major contributor to the "fusion power drops
        // to 0 at ~30 MW" bug — the bridge kept P_aux at 50 MW while
        // I_p ramped to 0, forcing a limit cycle).
        state.hcd_nbi_setpoint_MW  = 0.0f;
        state.hcd_icrh_setpoint_MW = 0.0f;
        state.hcd_ecrh_setpoint_MW = 0.0f;
        state.hcd_lhcd_setpoint_MW = 0.0f;
        state.hcd_nbi_on  = false;
        state.hcd_icrh_on = false;
        state.hcd_ecrh_on = false;
        state.hcd_lhcd_on = false;
        return;
    }

    // ── Controlled shutdown (Rampdown mode) ────────────────────────────────
    //  Operator-initiated, gentle ramp-down of EVERYTHING (I_p, T_e, ne,
    //  P_aux).  This is the other half of the "fusion power drops to 0"
    //  fix: previously only I_p and fuel_rate were zeroed on Rampdown, so
    //  T_e and ne setpoints stayed high — the bridge kept the plasma
    //  hot/dense while I_p collapsed, τ_E collapsed with it, and the
    //  power balance entered a 2-tick limit cycle between 30 MW and 0.
    //  Now we ramp down ALL four knobs simultaneously, keeping the plasma
    //  in a stable operating point as it decays.
    if (state.mode == ReactorMode::Rampdown) {
        runRampdown(state, dt);
        // Don't run normal PIDs during rampdown — the operator has taken
        // manual control via the rampdown ramps.  But still run safety
        // monitors so a quench or runaway is caught.
        runDisruptionWatch(state, dt);
        runRunawayMonitor(state);
        return;
    }

    runCurrentControl    (state, dt);
    runTemperatureControl(state, dt);
    runDensityControl    (state, dt);
    runShapeControl      (state, dt);
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
//  Output: auxiliary heating power command [MW] (NBI + ICRH + ECRH + LHCD total)
//
//  The PID output is written to state.sp_aux_heat_MW, which is then distributed
//  across the four H&CD systems by the HCD module's update() based on which are
//  enabled and their per-system setpoint caps.  This replaces the old hardcoded
//  P_aux = 50 MW in PlasmaCoreBridge — the bridge now reads sp_aux_heat_MW.
void ControlSystem::runTemperatureControl(ReactorState& state, float dt)
{
    aux_heat_command_MW_ = ControlPhysics::feedbackControlPID(
        state.sp_electron_temp_keV,
        state.electron_temp_keV,
        cfg_.Kp_temp, cfg_.Ki_temp, cfg_.Kd_temp,
        dt, pid_temp_,
        /*output_min=*/ 0.0f,    // can't have negative heating
        /*output_max=*/ 100.0f,  // 100 MW max aux heating (NBI + ICRH + ECRH + LHCD)
        /*windup=*/     1000.0f
    );
    // Write the unified aux heat setpoint to ReactorState.  The HCD module
    // picks this up on the next tick and distributes it across the enabled
    // H&CD systems (NBI/ICRH/ECRH/LHCD).
    state.sp_aux_heat_MW = aux_heat_command_MW_;
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

    // Beta limit check (normalised beta β_N)
    //  β_N = β[%] · a[m] · B[T] / I_p[MA], with Troyon limit at β_N = 2.5.
    //  The PlasmaCoreBridge already computes state.beta_N correctly, so we
    //  just use that here.  The old code recomputed beta_N with a wrong
    //  formula (divided by μ₀ and used I_p in Amps instead of MA), which
    //  gave beta_N values ~1000x too small — so the beta_limit check never
    //  fired, and the disruption-watch T_e clamp was effectively disabled.
    if (state.beta_N > cfg_.beta_limit) {
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
    trigger |= (state.alarm_vessel_breach);

    // Disruption that has already crossed q < 1.5 (not just predicted).
    // Only trigger if the toroidal field is at least 80% of setpoint —
    // during magnet ramp-up, B is low and q is artificially low (q ∝ B),
    // which would spuriously SCRAM before the field is established.
    bool field_ready = state.B_toroidal_T > 0.8f * state.sp_B_toroidal_T
                     && state.sp_B_toroidal_T > 0.5f;
    trigger |= (state.q_safety < 1.5f && state.plasma_current_MA > 1.0f && field_ready);

    if (trigger) {
        scram_latched_  = true;
        state.cmd_scram = true;
        state.mode      = ReactorMode::Emergency;
    }
}

// ─── Controlled Shutdown (Rampdown mode) ──────────────────────────────────────
//  Ramps ALL operating-point knobs (I_p, T_e, ne, P_aux) down at safe rates
//  so the plasma decays smoothly.  This is the fix for the third user-reported
//  issue's underlying physics problem — the old "CONTROLLED SHUTDOWN" button
//  only zeroed I_p and fuel_rate, which left T_e and ne setpoints high.  With
//  I_p collapsing but T_e/ne held at 20 keV / 1e20 m⁻³, the plasma power
//  balance entered a 2-tick limit cycle (W re-init from stale T_e → fusion
//  power spike → W collapse → fusion power 0 → repeat) — exactly the
//  "fusion power drops to 0 at ~30 MW then bounces back" bug.
//
//  By ramping ALL four knobs simultaneously, the plasma stays in a stable
//  operating point as it decays, and the power balance stays well-behaved.
void ControlSystem::runRampdown(ReactorState& state, float dt)
{
    // Ramp I_p setpoint down at the safe rampdown rate.
    // PlasmaCoreBridge will track this on the next tick.
    state.sp_plasma_current_MA = std::max(
        state.sp_plasma_current_MA - cfg_.rampdown_dIp_dt * dt, 0.0f);

    // Ramp T_e setpoint down (this also drops the aux heat PID output,
    // since the error goes negative).
    state.sp_electron_temp_keV = std::max(
        state.sp_electron_temp_keV - cfg_.rampdown_dTe_dt * dt, 0.0f);

    // Ramp density setpoint down
    state.sp_density_m3 = std::max(
        state.sp_density_m3 - cfg_.rampdown_dne_dt * dt, 0.0f);

    // Ramp fuel rate down (no point puffing fuel if density setpoint is 0)
    state.sp_fuel_rate = std::max(
        state.sp_fuel_rate - 0.5f * dt, 0.0f);

    // Ramp H&CD setpoints down — this is the other half of the "fusion
    // power drops to 0" fix.  Even with T_e setpoint dropping, if P_aux
    // stayed at 50 MW the bridge would keep reheating the plasma.  Each
    // H&CD system ramps down independently at the safe rate.
    float dP = cfg_.rampdown_dPaux_dt * dt;
    state.hcd_nbi_setpoint_MW  = std::max(state.hcd_nbi_setpoint_MW  - dP, 0.0f);
    state.hcd_icrh_setpoint_MW = std::max(state.hcd_icrh_setpoint_MW - dP, 0.0f);
    state.hcd_ecrh_setpoint_MW = std::max(state.hcd_ecrh_setpoint_MW - dP, 0.0f);
    state.hcd_lhcd_setpoint_MW = std::max(state.hcd_lhcd_setpoint_MW - dP, 0.0f);
    if (state.hcd_nbi_setpoint_MW  < 0.1f) state.hcd_nbi_on  = false;
    if (state.hcd_icrh_setpoint_MW < 0.1f) state.hcd_icrh_on = false;
    if (state.hcd_ecrh_setpoint_MW < 0.1f) state.hcd_ecrh_on = false;
    if (state.hcd_lhcd_setpoint_MW < 0.1f) state.hcd_lhcd_on = false;

    // Update the unified aux heat setpoint (sum of the four systems)
    state.sp_aux_heat_MW = state.hcd_nbi_setpoint_MW + state.hcd_icrh_setpoint_MW
                         + state.hcd_ecrh_setpoint_MW + state.hcd_lhcd_setpoint_MW;

    // Once everything has ramped to ~0, transition to Startup so the
    // operator can re-initiate via the INITIATE PLASMA button (the button
    // is gated on plasma_status == Cold).  We don't auto-Cold the plasma
    // — the bridge's early-return guard will Quench it once I_p drops
    // below 0.01 MA, and the operator's INITIATE button (gated on Cold)
    // then becomes the path back to a new discharge.
    if (state.sp_plasma_current_MA < 0.05f
        && state.plasma_current_MA  < 0.05f
        && state.sp_density_m3      < 1e17f
        && state.sp_electron_temp_keV < 0.5f) {
        state.mode = ReactorMode::Startup;
        if (state.plasma_status == PlasmaStatus::Burning
         || state.plasma_status == PlasmaStatus::Initiating) {
            state.plasma_status = PlasmaStatus::Cold;
        }
    }
}

// ─── Plasma Shape Controller ──────────────────────────────────────────────────
//  Drives the PF coil currents to track the operator's elongation (κ) and
//  triangularity (δ) setpoints.  Outputs the actual achieved κ/δ (which lag
//  the setpoint on a ~1 s PF coil current ramp) into state.kappa/delta —
//  these feed back into the q_95 calculation in PlasmaCoreBridge.
//
//  This is the shape-control loop a real tokamak has (a separate controller
//  from the plasma current loop).  Without it the operator had no way to
//  influence the plasma cross-section, which is one of the most important
//  operational knobs on a real machine (elongation trades confinement vs
//  vertical stability).
void ControlSystem::runShapeControl(ReactorState& state, float dt)
{
    // Simple first-order lag toward the setpoints — a real shape controller
    // would solve an inverse equilibrium problem, but for gameplay the lag
    // is what matters (the operator sees κ take ~1 s to settle after a
    // setpoint change, which matches the PF coil current ramp).
    constexpr float shape_tau = 1.0f;  // s
    float alpha = 1.0f - expf(-dt / shape_tau);

    state.kappa += (state.sp_kappa - state.kappa) * alpha;
    state.delta += (state.sp_delta - state.delta) * alpha;

    // R_out and Z_x aren't operator-settable in this version — they're
    // determined by the PF coil geometry.  We just hold them at ITER nominal.
    state.R_out_m = 6.2f;
    state.Z_x_m   = -1.0f * (state.kappa / 1.7f);  // scale with elongation

    // Update the plasma volume from the new shape
    // V = 2π² · R · a² · κ  (with a = 2.0 m fixed for ITER-class)
    constexpr float PI = 3.14159265358979f;
    constexpr float a_minor = 2.0f;
    state.plasma_volume_m3 = 2.0f * PI * PI * state.R_out_m
                           * a_minor * a_minor * state.kappa;
}