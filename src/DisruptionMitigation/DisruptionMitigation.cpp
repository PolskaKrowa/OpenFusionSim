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
        // Force the vessel pressure up (MGI releases gas into vessel)
        state.vessel_pressure_Pa = std::max(state.vessel_pressure_Pa,
                                              cfg_.mgi_pressure_rise_Pa);
        pumpdown_remaining_s_ = cfg_.post_mitigation_pumpdown_s;
    } else if (fire_spi) {
        spi_fired_ = true;
        mitigation_active_ = true;
        mitigation_remaining_s_ = cfg_.spi_thermal_quench_s;
        last_fire_type_ = "SPI";
        state.spi_fired = true;
        state.spi_pellet_mass_g = cfg_.spi_pellet_mass_g;
        state.mitigation_force_N = state.plasma_current_MA * 1e6f
                                 * cfg_.spi_halo_reduction;
        pumpdown_remaining_s_ = cfg_.post_mitigation_pumpdown_s * 0.3f;  // less gas
    }
    time_since_fire_s_ = 0.f;
}

void DisruptionMitigationSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    time_since_fire_s_ += dt;

    // Mirror armed state to ReactorState for UI
    state.mgi_armed = mgi_armed_;
    state.spi_armed = spi_armed_;

    // ── Auto-MGI: fire if q_95 drops below threshold ────────────────────────
    if (cfg_.auto_mgi_enabled && mgi_armed_ && !mgi_fired_
        && state.q_safety < cfg_.auto_mgi_q_threshold
        && state.plasma_current_MA > 1.0f) {
        fireMitigation(state);
    }

    // ── Operator-triggered mitigation (FIRE button) ─────────────────────────
    if (state.cmd_disrupt_mitigation) {
        fireMitigation(state);
        state.cmd_disrupt_mitigation = false;  // consume the command
    }

    // ── Active mitigation: thermal quench in progress ───────────────────────
    if (mitigation_active_) {
        mitigation_remaining_s_ -= dt;
        if (mitigation_remaining_s_ <= 0.f) {
            mitigation_active_ = false;
            // Force plasma into disruption / quench — the mitigation has
            // done its job (radiatively cooled the plasma), now the plasma
            // is gone.  PlasmaCoreBridge will pick up that I_p should be
            // quenched on the next tick.
            state.plasma_status = PlasmaStatus::Disrupting;
            state.disruption_flag = true;
            state.disruption_cause = DisruptionCause::None;  // operator-initiated
        } else {
            // During the thermal quench, force plasma to disrupting
            state.plasma_status = PlasmaStatus::Disrupting;
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
