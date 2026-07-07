//
// src/DisruptionMitigation/DisruptionMitigation.cpp
//  Implementation of MGI/SPI disruption mitigation systems.
//

#include "DisruptionMitigation.h"
#include <cmath>
#include <algorithm>

DisruptionMitigationSystem::DisruptionMitigationSystem(const DMConfig& cfg)
    : cfg_(cfg)
{}

void DisruptionMitigationSystem::reset()
{
    mgi_armed_ = false;
    spi_armed_ = false;
    mgi_fired_ = false;
    spi_fired_ = false;
    mitigation_active_ = false;
    mitigation_remaining_s_ = 0.f;
    time_since_fire_s_ = 1e9f;
    last_fire_type_ = "none";
    pumpdown_remaining_s_ = 0.f;
}

void DisruptionMitigationSystem::fireMitigation(ReactorState& state)
{
    if (mitigation_active_) return;  // already firing

    bool fire_mgi = mgi_armed_ && !mgi_fired_;
    bool fire_spi = spi_armed_ && !spi_fired_;

    if (fire_mgi) {
        mgi_fired_ = true;
        mitigation_active_ = true;
        mitigation_remaining_s_ = cfg_.mgi_thermal_quench_s;
        last_fire_type_ = "MGI";
        state.mgi_fired = true;
        state.mgi_pressure_Pa = cfg_.mgi_pressure_rise_Pa;
        state.mitigation_force_N = state.plasma_current_MA * 1e6f
                                 * cfg_.mgi_halo_reduction;
        // ── The actual mitigation physics ────────────────────────────────
        //  MGI works by flooding the plasma with a strong radiator (Ne/Ar):
        //  the injected impurity drives a radiative thermal quench through
        //  the line-radiation model in PlasmaCoreBridge.  Previously the
        //  module just flipped plasma_status for 1 ms and the bridge
        //  cleared the flag next tick — the plasma survived an MGI shot
        //  completely unharmed.
        state.dm_injected_impurity += 0.25f;   // massive Ne/Ar dose
        state.plasma_density_m3    *= 1.5f;    // injected neutrals fuel up
        state.dm_mitigated = true;             // gentle current quench
        // Request a pressure rise in the vessel.  VacuumVessel::update will
        // pick this up via forcePressureRise() and apply it to its internal
        // pressure_Pa_.  Previously we wrote state.vessel_pressure_Pa
        // directly here, but VacuumVessel::update overwrote it next tick.
        state.dm_pressure_rise_Pa = cfg_.mgi_pressure_rise_Pa;
        pumpdown_remaining_s_ = cfg_.post_mitigation_pumpdown_s;
        time_since_fire_s_ = 0.f;
    } else if (fire_spi) {
        spi_fired_ = true;
        mitigation_active_ = true;
        mitigation_remaining_s_ = cfg_.spi_thermal_quench_s;
        last_fire_type_ = "SPI";
        state.spi_fired = true;
        state.spi_pellet_mass_g = cfg_.spi_pellet_mass_g;
        state.mitigation_force_N = state.plasma_current_MA * 1e6f
                                 * cfg_.spi_halo_reduction;
        //  SPI: Ne-doped shattered D2 pellet — deeper penetration, better
        //  assimilation than MGI (denser radiator per mole injected) plus
        //  a large density rise from the deuterium fragments, which
        //  collisionally cools the plasma as well as radiating.
        state.dm_injected_impurity += 0.30f;
        state.plasma_density_m3    *= 2.0f;
        state.dm_mitigated = true;
        // SPI injects less gas than MGI (deuterium pellets, not noble gas)
        state.dm_pressure_rise_Pa = cfg_.mgi_pressure_rise_Pa * 0.3f;
        pumpdown_remaining_s_ = cfg_.post_mitigation_pumpdown_s * 0.3f;
        time_since_fire_s_ = 0.f;
    }
    //  NOTE: time_since_fire_s_ is only reset when a system actually fired
    //  (it used to reset unconditionally, so clicking FIRE with nothing
    //  armed lied about the last-fire time).
}

void DisruptionMitigationSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    time_since_fire_s_ += dt;

    // ── Arm state: the UI checkboxes are the source of truth ────────────────
    //  The operator arms/disarms via ImGui::Checkbox on state.mgi_armed /
    //  state.spi_armed.  The old code did the OPPOSITE — it overwrote the
    //  state flags from the (never-set) internals every tick, so ticking
    //  "Arm MGI" was reverted one frame later and neither system could
    //  ever be armed.
    mgi_armed_ = state.mgi_armed;
    spi_armed_ = state.spi_armed;
    // Fired latches flow the other way (module → UI lamps)
    state.mgi_fired = mgi_fired_;
    state.spi_fired = spi_fired_;

    // ── Auto-MGI enable follows the (persistent) UI flag ────────────────────
    cfg_.auto_mgi_enabled = state.dm_auto_mgi;

    // ── Auto-MGI: fire on detected disruption or q_95 threshold ─────────────
    //  This is the whole point of auto-mitigation: the machine-protection
    //  system reacts faster than any operator can.  Fires on either an
    //  actual disruption flag or the q_95 precursor threshold.
    if (cfg_.auto_mgi_enabled && mgi_armed_ && !mgi_fired_
        && state.plasma_current_MA > 1.0f
        && (state.disruption_flag
            || state.q_safety < cfg_.auto_mgi_q_threshold)) {
        fireMitigation(state);
    }

    // ── Operator-triggered mitigation (FIRE button) ─────────────────────────
    if (state.cmd_disrupt_mitigation) {
        fireMitigation(state);
        state.cmd_disrupt_mitigation = false;  // consume the command
    }

    // ── Active mitigation: injection window in progress ─────────────────────
    //  The thermal quench itself now happens through the radiation physics
    //  in PlasmaCoreBridge (the injected radiator collapses T_e); here we
    //  just track the injection window for the UI.
    if (mitigation_active_) {
        mitigation_remaining_s_ -= dt;
        if (mitigation_remaining_s_ <= 0.f) {
            mitigation_active_ = false;
        }
    }

    // ── Post-mitigation pumpdown ────────────────────────────────────────────
    if (pumpdown_remaining_s_ > 0.f) {
        pumpdown_remaining_s_ -= dt;
        // Pressure stays elevated during pumpdown (VacuumVessel module will
        // pump it back down — we just gate the next plasma-initiation on
        // state.vessel_vacuum_ok, which VacuumVessel updates)
    }

    // ── Reset fired latches on operator RESET ───────────────────────────────
    // (The operator's RESET — COLD RESTART path calls reset() on every
    // module, so the latches clear there.  No additional logic needed here.)
}
