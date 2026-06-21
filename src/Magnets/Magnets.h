#pragma once

//
// src/Magnets/Magnets.h
// Superconducting magnet system.
//
//  Subsystems modelled:
//    - TF (toroidal field) coils — set the main confinement field.
//    - PF (poloidal field) coils — shape and position the plasma.
//    - CS (central solenoid)    — drives the plasma current via transformer action.
//    - Quench detection and protection (QPS).
//    - Cryogenic plant (4.5 K helium cooling).
//

#include "ReactorState.h"
#include "SimTime.h"

struct MagnetConfig {
    // TF coil system (18 coils for ITER-class)
    float max_B_T          = 11.8f;   // peak field on conductor [T]
    float nominal_B_axis_T =  5.3f;   // field on plasma axis [T]
    float TF_stored_GJ     = 41.0f;   // stored magnetic energy [GJ]

    // Superconductor parameters (Nb3Sn)
    float T_critical_K     =  18.0f;  // critical temp [K]
    float T_op_K           =   4.5f;  // operating temp [K]
    float J_critical_Am2   =  500e6f; // critical current density [A/m^2]

    // Quench protection
    float quench_detect_dt_s = 5e-3f;    // detection time [s]
    float dump_resistor_ohm  = 0.3f;     // energy dump resistor [Ω]
    float current_decay_tau  = 10.0f;    // L/R time constant [s]
};

// ─── Individual coil state ─────────────────────────────────────────────────────
struct CoilState {
    float current_kA   = 0.0f;
    float temp_K       = 4.5f;
    bool  quenched     = false;
    float resistance_uOhm = 0.0f; // SC: ~0; normal: large
};

class MagnetSystem {
public:
    explicit MagnetSystem(const MagnetConfig& cfg);

    // Called each game tick.  Reads setpoints, writes B fields and quench flags.
    void update(ReactorState& state, const SimTime& t);

    // Manually trigger an emergency dump (called by Control on SCRAM)
    void triggerQuenchDump();

    // Clear the latched quench-dump flag WITHOUT resetting coil currents
    // or temperatures.  This is what the SCRAM popup's ACKNOWLEDGE button
    // needs to call: once the QPS has fired and the energy dump is
    // complete, the quench_detected flag stays latched forever (it's
    // only cleared by reset()).  But ACKNOWLEDGE is meant to let the
    // operator resume operation without a full cold restart — so we
    // need to clear just the quench latch, not the entire magnet state.
    //
    // Also clears the per-coil quenched flags so the coils can carry
    // current again.
    void clearQuenchLatch();

    // Cold-restart: clear all internal coil state, quench flags and the dump
    // latch.  Without this, the dump_triggered_ flag from a previous SCRAM
    // keeps state.quench_detected=true forever and the control system
    // re-trips SCRAM on the very next tick (classic soft-lock).
    void reset();

    const MagnetConfig& config() const { return cfg_; }

private:
    void  updateTFCoils  (ReactorState& state, float dt);
    void  updateCSCoil   (ReactorState& state, float dt);
    void  runQPS         (ReactorState& state, float dt);
    void  coolCryoplant  (float dt);
    float fieldOnAxis    (float coil_current_kA) const;

    MagnetConfig cfg_;

    // TF coil array (18 coils)
    static constexpr int N_TF = 18;
    CoilState tf_coils_[N_TF];

    // CS solenoid (6 modules)
    static constexpr int N_CS = 6;
    CoilState cs_coils_[N_CS];

    float cryo_load_W_       = 0.0f;  // cryoplant heat load [W]
    float current_ramp_rate_ = 0.0f;  // kA/s
    bool  dump_triggered_    = false;
    // Previous-tick total CS current (for loop voltage dI/dt calculation).
    //  Reset to 0 in reset() so the first tick after reset doesn't produce
    //  a huge dI/dt spike.
    float cs_current_prev_kA_ = 0.0f;
};