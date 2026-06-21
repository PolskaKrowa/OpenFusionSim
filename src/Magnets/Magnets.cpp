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

    // Write aggregate outputs to state.
    //  All TF coils are in series (carry the same current), so reading
    //  tf_coils_[0] is representative of all 18.  This is documented in
    //  the MagnetConfig header.
    state.B_toroidal_T    = fieldOnAxis(tf_coils_[0].current_kA);
    state.magnet_temp_K   = tf_coils_[0].temp_K;
    state.stored_energy_GJ = 0.5f * cfg_.TF_stored_GJ
                           * (tf_coils_[0].current_kA / 68.0f)   // 68 kA nominal
                           * (tf_coils_[0].current_kA / 68.0f);
    state.coil_current_kA = tf_coils_[0].current_kA;

    // ── Central Solenoid summary fields (H1) ────────────────────────────────
    //  Sum the 6 CS module currents to get the total CS current, and estimate
    //  the loop voltage from the CS current ramp rate.  The loop voltage is
    //  what drives the plasma current via transformer action: V_loop = -dΦ/dt
    //  where Φ = L_CS × I_CS.  With L_CS ≈ 6 H (ITER-class) and I_CS in kA:
    //    V_loop = -L × dI/dt  [V]
    //  Previously these fields were declared in ReactorState.h but never
    //  written — the UI showed "CS Current: 0.0 kA" permanently.
    float cs_total_kA = 0.f;
    for (const auto& c : cs_coils_) cs_total_kA += c.current_kA;
    state.cs_current_kA = cs_total_kA;
    // Loop voltage: derivative of CS current × CS inductance.
    //  Track previous-tick CS current to compute dI/dt.
    constexpr float L_CS_H = 6.0f;  // ITER CS inductance [H]
    float dI_cs_dt = (cs_total_kA - cs_current_prev_kA_) / std::max(t.dt_s, 1e-6f);
    state.cs_loop_voltage_V = -L_CS_H * dI_cs_dt * 1e3f;  // kA→A: ×1e3
    cs_current_prev_kA_ = cs_total_kA;
    // PF coil current (simplified: assume PF mirrors CS for shape control)
    state.pf_current_kA = cs_total_kA * 0.5f;
    // CS flux remaining: ITER CS has ~90 Wb swing; decreases as CS ramps up
    state.cs_flux_remaining_Wb = std::max(0.f, 90.f - cs_total_kA * 2.0f);

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

    // Poloidal field from CS current (peak ~6 T at 45 kA).
    //  This is the only field written here — the loop voltage is now computed
    //  in update() from the CS current ramp rate (using cs_current_prev_kA_).
    //  The old code computed a dead `dI_dt` here and discarded it via (void).
    state.B_poloidal_T = cs_coils_[0].current_kA / 45.0f * 6.0f;
}

void MagnetSystem::runQPS(ReactorState& state, float dt)
{
    // Quench detection: any TF or CS coil exceeding T_critical.
    //  Previously only TF coils were checked — a CS quench was never
    //  detected even though CS coils carry significant current and can
    //  quench under abnormal conditions.
    for (auto& c : tf_coils_) {
        if (c.temp_K > cfg_.T_critical_K * 0.9f && !c.quenched) {
            c.quenched        = true;
            c.resistance_uOhm = 1000.0f; // transition to normal state
            dump_triggered_   = true;
        }
    }
    for (auto& c : cs_coils_) {
        if (c.temp_K > cfg_.T_critical_K * 0.9f && !c.quenched) {
            c.quenched        = true;
            c.resistance_uOhm = 1000.0f;
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

void MagnetSystem::clearQuenchLatch()
{
    // Clear the dump latch and per-coil quenched flags, but leave coil
    // currents and temperatures as-is.  This lets the operator resume
    // operation after a SCRAM without a full cold restart — the magnets
    // can start ramping back up from their current state.
    //
    // The quenched coils have their resistance zeroed so they can carry
    // supercurrent again.  (In reality, a quenched coil would need to be
    // inspected before re-energizing, but for gameplay purposes we allow
    // it — the operator is responsible for verifying the magnets are OK.)
    dump_triggered_ = false;
    for (auto& c : tf_coils_) {
        c.quenched       = false;
        c.resistance_uOhm = 0.0f;
    }
    for (auto& c : cs_coils_) {
        c.quenched       = false;
        c.resistance_uOhm = 0.0f;
    }
}

void MagnetSystem::reset()
{
    // Restore every coil to its as-built state and forget any quench/dump
    // history.  This is what the operator's RESET — COLD RESTART button
    // needs to call: the previous SCRAM latched dump_triggered_=true, and
    // that flag would otherwise keep state.quench_detected=true forever
    // (which runScramLogic reads as "magnets are quenched" and re-trips
    // the SCRAM on the next tick, trapping the sim in a soft-lock).
    for (auto& c : tf_coils_) {
        c.current_kA     = 0.0f;
        c.temp_K         = cfg_.T_op_K;
        c.quenched       = false;
        c.resistance_uOhm = 0.0f;
    }
    for (auto& c : cs_coils_) {
        c.current_kA     = 0.0f;
        c.temp_K         = cfg_.T_op_K;
        c.quenched       = false;
        c.resistance_uOhm = 0.0f;
    }
    dump_triggered_   = false;
    cryo_load_W_      = 0.0f;
    current_ramp_rate_ = 0.0f;
    cs_current_prev_kA_ = 0.0f;
}