//
// src/Magnets/Magnets.cpp
//

#include "Magnets.h"
#include <cmath>
#include <algorithm>

MagnetSystem::MagnetSystem(const MagnetConfig& cfg)
    : cfg_(cfg)
{
    for (auto& c : tf_coils_) { c.temp_K = cfg_.T_op_K; }
    for (auto& c : cs_coils_) { c.temp_K = cfg_.T_op_K; }
}

void MagnetSystem::update(ReactorState& state, const SimTime& t)
{
    if (state.cmd_scram && !dump_triggered_)
        triggerQuenchDump();

    updateTFCoils(state, t.dt_s);
    updateCSCoil (state, t.dt_s);
    runQPS       (state, t.dt_s);
    coolCryoplant(t.dt_s);

    // Write aggregate outputs to state
    state.B_toroidal_T    = fieldOnAxis(tf_coils_[0].current_kA);
    state.magnet_temp_K   = tf_coils_[0].temp_K;
    state.stored_energy_GJ = 0.5f * cfg_.TF_stored_GJ
                           * (tf_coils_[0].current_kA / 68.0f)   // 68 kA nominal
                           * (tf_coils_[0].current_kA / 68.0f);
    state.coil_current_kA = tf_coils_[0].current_kA;
    state.quench_detected = dump_triggered_;
    state.alarm_quench    = dump_triggered_;
}

void MagnetSystem::updateTFCoils(ReactorState& state, float dt)
{
    // Ramp TF coil current toward setpoint
    float target_kA = (state.sp_B_toroidal_T / cfg_.nominal_B_axis_T) * 68.0f;
    float max_ramp  = 2.0f; // kA/s — superconductor ramp rate limit

    for (auto& c : tf_coils_) {
        if (c.quenched) continue;
        if (dump_triggered_) {
            // L/R decay during dump
            c.current_kA *= expf(-dt / cfg_.current_decay_tau);
        } else {
            float delta = std::clamp(target_kA - c.current_kA,
                                     -max_ramp * dt, max_ramp * dt);
            c.current_kA += delta;
        }

        // Ohmic heating from joint resistance: Q = I^2 * R
        float R_joint = 0.5e-9f; // 0.5 nΩ per joint, typical
        float I = c.current_kA * 1e3f;
        float Q_W = I * I * R_joint;
        // Cryoplant removes heat; temp fluctuates slightly
        c.temp_K += (Q_W / 500.0f) * dt; // 500 J/K heat capacity (rough)
        c.temp_K -= (c.temp_K - cfg_.T_op_K) * 0.1f * dt; // cryo feedback
    }
}

void MagnetSystem::updateCSCoil(ReactorState& state, float dt)
{
    // CS solenoid drives plasma current via V_loop = -dΦ/dt
    // Simplified: ramp CS current to induce target Ip
    float target_kA = state.sp_plasma_current_MA * 8.0f; // ~8 kA-turns per MA Ip
    target_kA = std::clamp(target_kA, 0.0f, 45.0f);

    for (auto& c : cs_coils_) {
        if (c.quenched || dump_triggered_) { c.current_kA = 0; continue; }
        float delta = std::clamp(target_kA - c.current_kA, -5.0f * dt, 5.0f * dt);
        c.current_kA += delta;
        c.temp_K += (c.current_kA * c.current_kA * 1e-9f / 200.0f) * dt;
        c.temp_K -= (c.temp_K - cfg_.T_op_K) * 0.1f * dt;
    }

    // Approximate induced loop voltage
    float dI_dt = (cs_coils_[0].current_kA - target_kA) / dt;
    state.B_poloidal_T = cs_coils_[0].current_kA / 45.0f * 6.0f; // peak ~6 T at 45 kA
    (void)dI_dt;
}

void MagnetSystem::runQPS(ReactorState& state, float dt)
{
    // Quench detection: any coil exceeding T_critical or developing resistance
    for (auto& c : tf_coils_) {
        if (c.temp_K > cfg_.T_critical_K * 0.9f && !c.quenched) {
            c.quenched        = true;
            c.resistance_uOhm = 1000.0f; // transition to normal state
            dump_triggered_   = true;
        }
    }
    (void)dt;
    (void)state;
}

void MagnetSystem::coolCryoplant(float dt)
{
    // Very simple cryoplant model: drives all coil temps toward T_op
    // Real system: 60 kW of cryoplant cooling for ITER
    cryo_load_W_ = 0.0f;
    for (auto& c : tf_coils_)
        cryo_load_W_ += std::max(0.0f, c.temp_K - cfg_.T_op_K) * 500.0f;
    (void)dt;
}

float MagnetSystem::fieldOnAxis(float coil_current_kA) const
{
    // Linear scaling: 68 kA → 5.3 T on axis (ITER-class)
    return cfg_.nominal_B_axis_T * (coil_current_kA / 68.0f);
}

void MagnetSystem::triggerQuenchDump()
{
    dump_triggered_ = true;
}