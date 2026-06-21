#pragma once
//
// src/VacuumVessel/VacuumVessel.h
// Vacuum vessel pumping and conditioning system.
//
//  The vacuum vessel must be pumped down from atmospheric pressure to a base
//  pressure of < 1e-4 Pa before plasma initiation.  This module models:
//
//    - Roughing pump (rotary vane): pulls from atmosphere (~101325 Pa) down
//      to ~1 Pa.  Speed: ~1000 m³/h = 0.28 m³/s.  Cannot operate below 1 Pa
//      (back-streaming oil).
//
//    - Turbo-molecular pump: pulls from 1 Pa down to 1e-5 Pa.  Cannot
//      operate above 10 Pa (bearing damage).  The operator must wait for
//      the roughing pump to drop the pressure before starting the turbo.
//
//    - Wall bakeout: heats the vessel wall to 423 K (150 °C) for ~24 h to
//      desorb water and other gases.  Without bakeout the outgassing load
//      exceeds the turbo pump capacity and the base pressure plateaus at
//      ~1e-3 Pa (too high for plasma initiation).
//
//    - Wall conditioning (boronization): deposits a thin boron film on the
//      first wall to getter oxygen and reduce impurity radiation.  Performed
//      periodically as the film erodes.
//
//  The PlasmaCoreBridge is interlocked: it won't transition out of Cold
//  status until vessel_pressure_Pa < 1e-3 Pa.  The UI's INITIATE PLASMA
//  button checks this before allowing initiation.
//

#include "ReactorState.h"
#include "SimTime.h"
#include <string>

struct VacuumConfig {
    // Roughing pump
    float roughing_speed_m3s   = 50.0f;
    float roughing_min_Pa      = 1.0f;     // can't operate below this
    float roughing_max_Pa      = 110000.f; // can't operate above this

    // Turbo pump
    float turbo_speed_m3s      = 100.0f;   // 14400 m³/h (large turbopump)
    float turbo_max_Pa         = 10.0f;    // bearing damage above this
    float turbo_min_Pa         = 1e-7f;    // ultimate vacuum

    // Vessel volume and surface area (ITER-class)
    float vessel_volume_m3     = 1400.f;   // vacuum vessel interior
    float vessel_surface_m2    = 800.f;    // inner wall surface

    // Wall outgassing (after bakeout)
    float outgas_rate_Pa_m3_s_m2 = 1e-7f;  // baked wall: 1e-7 Pa·m³/(s·m²)
    float outgas_rate_unbaked    = 1e-4f;  // unbaked wall: 1e-4 Pa·m³/(s·m²)
                                            // (water vapor dominates)

    // Bakeout
    float bakeout_temp_K       = 423.f;    // 150 °C
    //  Real tokamak wall bakeout takes ~24 h at 150 °C to desorb water vapor
    //  from the vessel walls.  For gameplay we compress this to 5 minutes —
    //  still long enough that the operator has to plan around it (you can't
    //  just spam-bake between every discharge), but short enough that it
    //  doesn't kill the pace of the game.  The physics is the same: the
    //  wall outgassing factor decays exponentially during bakeout, just on
    //  a faster timescale.
    float bakeout_duration_s   = 5.f * 60.f;  // 5 min (gameplay-compressed from 24 h)

    // Plasma-initiation interlock
    float initiation_max_Pa    = 1e-3f;    // must be below this for INITIATE

    // Vessel breach (loss-of-vacuum accident)
    float breach_max_Pa        = 1000.f;   // if pressure rises above this
                                            // without pumps, suspect breach
};

class VacuumVesselSystem {
public:
    explicit VacuumVesselSystem(const VacuumConfig& cfg = {});

    void update(ReactorState& state, const SimTime& t);

    // Cold-restart: pressure back to atmospheric, pumps off, no bakeout
    void reset();

    // Operator commands (also exposed via ReactorState flags)
    void startRoughing() { roughing_on_ = true; }
    void stopRoughing()  { roughing_on_ = false; }
    void startTurbo();
    void stopTurbo()     { turbo_on_ = false; }
    void startBakeout()  { bakeout_on_ = true; }
    void stopBakeout()   { bakeout_on_ = false; }
    void startBoronization();
    void triggerVesselBreach();
    void clearBreach();

    // Force the internal pressure to at least `Pa` (used by the Disruption
    // Mitigation module when MGI fires — the gas injection raises vessel
    // pressure to ~10 Pa).  Without this, MGI writes state.vessel_pressure_Pa
    // directly, but VacuumVessel::update overwrites it next tick with the
    // internal pressure_Pa_ (which wasn't updated), so the MGI pressure
    // effect was lost.
    void forcePressureRise(float Pa) { pressure_Pa_ = std::max(pressure_Pa_, Pa); }

    // Diagnostic accessors
    float pressure() const { return pressure_Pa_; }
    bool  vacuumOK() const { return pressure_Pa_ < cfg_.initiation_max_Pa; }
    bool  bakeoutComplete() const { return bakeout_complete_; }
    float bakeoutProgress() const;  // 0..1
    bool  breachDetected() const { return breach_detected_; }
    const std::string& breachReason() const { return breach_reason_; }

private:
    VacuumConfig cfg_;

    // Vessel state
    float pressure_Pa_         = 101325.f;  // atmospheric at start
    float wall_temp_K_         = 300.f;     // room temp
    float wall_outgas_factor_  = 1.0f;      // 1 = unbaked; decays toward 0 with bakeout

    // Pump states (mirrored to ReactorState)
    bool  roughing_on_         = false;
    bool  turbo_on_            = false;
    bool  bakeout_on_          = false;
    bool  bakeout_complete_    = false;
    float bakeout_elapsed_s_   = 0.f;

    // Boronization state
    float boronization_thickness_nm_ = 0.f;
    bool  boronization_in_progress_  = false;
    float boronization_elapsed_s_    = 0.f;
    constexpr static float boronization_duration_s = 4.f * 3600.f;  // 4 h

    // Vessel breach
    bool  breach_detected_     = false;
    std::string breach_reason_;
    // Previous-tick pressure for breach detection (dP/dt).  This was
    // previously a function-static in checkBreach, which meant it survived
    // across ReactorState resets — after a SCRAM+RESET, pressure_Pa_ jumps
    // back to 101325 Pa but prev_pressure_ still held the pre-reset value
    // (could be 1e-5 Pa), producing a spurious dP/dt of ~1e8 Pa/s and a
    // false breach alarm.  Now it's a member that reset() clears.
    float prev_pressure_Pa_    = 101325.f;

    // Internal helpers
    float effectiveOutgasRate() const;
    void  updateBakeout(float dt);
    void  updateBoronization(float dt);
    void  updatePumping(float dt, float plasma_load_Pa_m3_s);
    void  checkBreach(float dt);
};
