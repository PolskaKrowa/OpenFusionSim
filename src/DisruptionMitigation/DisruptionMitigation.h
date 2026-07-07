#pragma once
//
// src/DisruptionMitigation/DisruptionMitigation.h
// Disruption mitigation systems: Massive Gas Injection (MGI) and
// Shattered Pellet Injection (SPI).
//
//  When a disruption is detected or imminent, the operator can fire one of
//  two mitigation systems to reduce the damage:
//
//    MGI — Massive Gas Injection
//           A fast-acting valve opens, releasing a large quantity of noble
//           gas (Ne or Ar) into the vacuum vessel.  The gas radiatively
//           cools the plasma on a ~1 ms timescale, forcing a controlled
//           thermal quench rather than a localized one.  Reduces:
//             - Halo current forces (factor ~3 reduction)
//             - Runaway electron generation (collisional damping)
//             - First-wall heat flux (radiation distributes the load)
//           Downside: the vessel pressure rises to ~10 Pa, so the next
//           discharge can't start until the pumps recover the vacuum.
//
//    SPI — Shattered Pellet Injection
//           A large deuterium (or D-T) pellet is fired at high speed into
//           a shatter ring, producing a spray of solid fragments that
//           disperse through the plasma.  More effective than MGI at high
//           plasma current (ITER-class), and the deuterium doesn't pollute
//           the vessel like neon/argon does.  Used as the primary
//           mitigation on ITER.
//
//  Both systems require the operator to ARM them first (they're not
//  auto-fired — the operator makes the decision).  The DisruptionMitigation
//  module consumes the cmd_disrupt_mitigation flag from ReactorState and
//  fires whichever system is armed (MGI takes priority if both are armed).
//
//  The mitigation effect is fed back into the PlasmaCoreBridge via
//  state.mitigation_force_N (peak halo force reduction) and a brief
//  "mitigation in progress" window during which the plasma is forced into
//  Disrupting status.
//

#include "ReactorState.h"
#include "SimTime.h"
#include <string>

struct DMConfig {
    // MGI
    float mgi_gas_amount_mol   = 10.f;    // moles of Ne injected per shot
    float mgi_injection_time_s = 1e-3f;   // 1 ms valve open time
    float mgi_pressure_rise_Pa = 10.f;    // vessel pressure after MGI
    float mgi_halo_reduction   = 0.66f;   // 66% halo current reduction
    float mgi_thermal_quench_s = 1e-3f;   // 1 ms thermal quench

    // SPI
    float spi_pellet_mass_g     = 0.5f;   // 500 mg deuterium pellet
    float spi_pellet_velocity_ms= 250.f;  // m/s
    float spi_injection_time_s  = 2e-3f;  // 2 ms transit time
    float spi_halo_reduction    = 0.75f;  // 75% halo current reduction
    float spi_thermal_quench_s  = 2e-3f;  // 2 ms thermal quench

    // Recovery
    float post_mitigation_pumpdown_s = 30.f;  // 30 s to recover vacuum after MGI

    // Auto-fire threshold (operator can enable "auto-MGI on disruption")
    bool  auto_mgi_enabled = false;
    float auto_mgi_q_threshold = 1.8f;  // auto-fire if q_95 drops below this
};

class DisruptionMitigationSystem {
public:
    explicit DisruptionMitigationSystem(const DMConfig& cfg = {});

    void update(ReactorState& state, const SimTime& t);

    // Cold-restart: disarm both systems, clear all latches
    void reset();

    // Operator commands (also exposed via ReactorState flags)
    void armMGI()  { mgi_armed_ = true; }
    void disarmMGI() { mgi_armed_ = false; }
    void armSPI()  { spi_armed_ = true; }
    void disarmSPI() { spi_armed_ = false; }

    // Fire whichever system is armed (called when operator clicks FIRE or
    // when the auto-MGI logic detects a disruption threshold crossing).
    void fireMitigation(ReactorState& state);

    // Enable / disable auto-MGI
    void enableAutoMGI(bool on) { cfg_.auto_mgi_enabled = on; }
    bool autoMGIEnabled() const { return cfg_.auto_mgi_enabled; }

    // Diagnostic accessors
    bool  mgiArmed()  const { return mgi_armed_; }
    bool  spiArmed()  const { return spi_armed_; }
    bool  mgiFired()  const { return mgi_fired_; }
    bool  spiFired()  const { return spi_fired_; }
    bool  active()    const { return mitigation_active_; }
    float timeSinceLastFire_s() const { return time_since_fire_s_; }
    const std::string& lastFireType() const { return last_fire_type_; }

private:
    DMConfig cfg_;

    bool  mgi_armed_ = false;
    bool  spi_armed_ = false;
    bool  mgi_fired_ = false;  // latched — needs reset to fire again
    bool  spi_fired_ = false;

    bool  mitigation_active_ = false;
    float mitigation_remaining_s_ = 0.f;
    float time_since_fire_s_ = 1e9f;  // large = "never fired"
    std::string last_fire_type_ = "none";

    // State of the post-mitigation pumpdown
    float pumpdown_remaining_s_ = 0.f;
};
