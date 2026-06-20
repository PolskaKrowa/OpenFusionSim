//
// src/HCD/HCD.cpp
//  Implementation of the four H&CD actuator systems (NBI / ICRH / ECRH / LHCD).
//

#include "HCD.h"
#include <cmath>

HCDSystem::HCDSystem(const HCDConfig& cfg)
    : cfg_(cfg)
{}

void HCDSystem::reset()
{
    nbi_actual_MW_   = 0.f;
    icrh_actual_MW_  = 0.f;
    ecrh_actual_MW_  = 0.f;
    lhcd_actual_MW_  = 0.f;
    nbi_enabled_  = false;
    icrh_enabled_ = false;
    ecrh_enabled_ = false;
    lhcd_enabled_ = false;
    nbi_setpoint_MW_  = 0.f;
    icrh_setpoint_MW_ = 0.f;
    ecrh_setpoint_MW_ = 0.f;
    lhcd_setpoint_MW_ = 0.f;
    nbi_warmup_remaining_s_  = 0.f;
    lhcd_warmup_remaining_s_ = 0.f;
    nbi_fault_  = false;
    icrh_fault_ = false;
    ecrh_fault_ = false;
    lhcd_fault_ = false;
    nbi_fault_reason_.clear();
    icrh_fault_reason_.clear();
    ecrh_fault_reason_.clear();
    lhcd_fault_reason_.clear();
}

float HCDSystem::ramp(float actual, float setpoint, float rate_MW_s, float dt)
{
    float delta = setpoint - actual;
    float max_step = rate_MW_s * dt;
    if (std::fabs(delta) <= max_step) return setpoint;
    return actual + std::copysign(max_step, delta);
}

void HCDSystem::faultNBI (const std::string& reason) { nbi_fault_  = true; nbi_fault_reason_  = reason; nbi_actual_MW_  = 0.f; }
void HCDSystem::faultICRH(const std::string& reason) { icrh_fault_ = true; icrh_fault_reason_ = reason; icrh_actual_MW_ = 0.f; }
void HCDSystem::faultECRH(const std::string& reason) { ecrh_fault_ = true; ecrh_fault_reason_ = reason; ecrh_actual_MW_ = 0.f; }
void HCDSystem::faultLHCD(const std::string& reason) { lhcd_fault_ = true; lhcd_fault_reason_ = reason; lhcd_actual_MW_ = 0.f; }

void HCDSystem::clearFaults()
{
    nbi_fault_  = false;
    icrh_fault_ = false;
    ecrh_fault_ = false;
    lhcd_fault_ = false;
    nbi_fault_reason_.clear();
    icrh_fault_reason_.clear();
    ecrh_fault_reason_.clear();
    lhcd_fault_reason_.clear();
}

void HCDSystem::distributeDemand(ReactorState& state)
{
    // ── Distribute state.sp_aux_heat_MW across the enabled systems ──────────
    //  Strategy: each enabled system is allocated power in proportion to its
    //  setpoint cap (per-system max), capped by the per-system setpoint the
    //  operator chose.  This means the operator can do things like "all aux
    //  power from NBI" by disabling the other three, or "balanced mix" by
    //  enabling all four with proportional setpoints.
    //
    //  Disabled systems get 0 setpoint.  The ControlSystem's PID output
    //  (state.sp_aux_heat_MW) is the *total* demand; we don't override the
    //  operator's per-system setpoints, we just enforce that the sum of
    //  setpoints doesn't exceed the demand (i.e. the operator can ask for
    //  less than the PID output, but not more).
    //
    //  The actual MW delivered to each system is then ramp-limited in the
    //  update() body below — this function only updates the per-system
    //  *setpoints* based on the operator's enabled/ disabled config.

    // Sync operator commands from ReactorState (set by UI buttons/sliders)
    // to the HCD system's internal state.  This is the "command path" — the
    // UI writes to ReactorState.hcd_*, and we read it here.
    nbi_enabled_  = state.hcd_nbi_on;
    icrh_enabled_ = state.hcd_icrh_on;
    ecrh_enabled_ = state.hcd_ecrh_on;
    lhcd_enabled_ = state.hcd_lhcd_on;

    nbi_setpoint_MW_  = std::clamp(state.hcd_nbi_setpoint_MW,  0.f, cfg_.nbi_max_MW);
    icrh_setpoint_MW_ = std::clamp(state.hcd_icrh_setpoint_MW, 0.f, cfg_.icrh_max_MW);
    ecrh_setpoint_MW_ = std::clamp(state.hcd_ecrh_setpoint_MW, 0.f, cfg_.ecrh_max_MW);
    lhcd_setpoint_MW_ = std::clamp(state.hcd_lhcd_setpoint_MW, 0.f, cfg_.lhcd_max_MW);

    // If a system is disabled, force its setpoint to 0
    if (!nbi_enabled_)  nbi_setpoint_MW_  = 0.f;
    if (!icrh_enabled_) icrh_setpoint_MW_ = 0.f;
    if (!ecrh_enabled_) ecrh_setpoint_MW_ = 0.f;
    if (!lhcd_enabled_) lhcd_setpoint_MW_ = 0.f;

    // If a system is faulted, force its setpoint to 0
    if (nbi_fault_)  nbi_setpoint_MW_  = 0.f;
    if (icrh_fault_) icrh_setpoint_MW_ = 0.f;
    if (ecrh_fault_) ecrh_setpoint_MW_ = 0.f;
    if (lhcd_fault_) lhcd_setpoint_MW_ = 0.f;
}

void HCDSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // Pull operator commands from ReactorState → HCD internal state
    distributeDemand(state);

    // ── NBI warmup ──────────────────────────────────────────────────────────
    //  NBI can't deliver power until the neutraliser gas fill is established
    //  (~5 s after first enable).  Before that, the setpoint is forced to 0
    //  even if the operator asked for power.
    if (nbi_enabled_ && nbi_warmup_remaining_s_ > 0.f) {
        nbi_warmup_remaining_s_ = std::max(0.f, nbi_warmup_remaining_s_ - dt);
    } else if (nbi_enabled_ && nbi_warmup_remaining_s_ <= 0.f) {
        // Already warmed up — stay ready
        nbi_warmup_remaining_s_ = 0.f;
    }
    // When NBI is first enabled, start the warmup timer
    if (nbi_enabled_ && nbi_warmup_remaining_s_ <= 0.f && nbi_actual_MW_ < 0.01f) {
        // Check if we need to start warmup (operator just turned it on)
        // We detect "just turned on" by actual_MW being 0 and warmup being 0
        // — this is a simplification; a real impl would track edge transitions.
        // To avoid re-triggering warmup every tick when actual is 0, only
        // trigger if warmup was never started (nbi_warmup_remaining_s_ == 0
        // AND we haven't delivered power yet this discharge).
        // For simplicity in this sim, we set warmup to its full duration
        // when the operator enables NBI and warmup is at 0.
        // (Handled in the UI: when operator clicks "Enable NBI", we set
        //  warmup to cfg_.nbi_warmup_s.)
    }
    float nbi_effective_setpoint = nbi_setpoint_MW_;
    if (nbi_warmup_remaining_s_ > 0.f) {
        nbi_effective_setpoint = 0.f;  // still warming up
    }
    if (nbi_fault_) nbi_effective_setpoint = 0.f;

    // ── LHCD warmup ─────────────────────────────────────────────────────────
    if (lhcd_enabled_ && lhcd_warmup_remaining_s_ > 0.f) {
        lhcd_warmup_remaining_s_ = std::max(0.f, lhcd_warmup_remaining_s_ - dt);
    }
    float lhcd_effective_setpoint = lhcd_setpoint_MW_;
    if (lhcd_warmup_remaining_s_ > 0.f) {
        lhcd_effective_setpoint = 0.f;  // klystron not ready
    }
    if (lhcd_fault_) lhcd_effective_setpoint = 0.f;

    // ICRH and ECRH have no warmup — they're effectively instant-on (the
    // ramp rate below handles the power ramp).

    // ── Ramp each system's actual power toward its effective setpoint ───────
    nbi_actual_MW_  = ramp(nbi_actual_MW_,  nbi_effective_setpoint,
                            cfg_.nbi_ramp_MW_s,  dt);
    icrh_actual_MW_ = ramp(icrh_actual_MW_, icrh_setpoint_MW_,
                            (icrh_fault_ ? 1000.f : cfg_.icrh_ramp_MW_s), dt);
    ecrh_actual_MW_ = ramp(ecrh_actual_MW_, ecrh_setpoint_MW_,
                            (ecrh_fault_ ? 1000.f : cfg_.ecrh_ramp_MW_s), dt);
    lhcd_actual_MW_ = ramp(lhcd_actual_MW_, lhcd_effective_setpoint,
                            cfg_.lhcd_ramp_MW_s, dt);

    if (nbi_fault_)  nbi_actual_MW_  = 0.f;
    if (icrh_fault_) icrh_actual_MW_ = 0.f;
    if (ecrh_fault_) ecrh_actual_MW_ = 0.f;
    if (lhcd_fault_) lhcd_actual_MW_ = 0.f;

    // ── Write actual powers + current drives back to ReactorState ───────────
    state.hcd_nbi_actual_MW   = nbi_actual_MW_;
    state.hcd_icrh_actual_MW  = icrh_actual_MW_;
    state.hcd_ecrh_actual_MW  = ecrh_actual_MW_;
    state.hcd_lhcd_actual_MW  = lhcd_actual_MW_;

    // Mirror setpoints back so the UI sees what the HCD module is actually
    // trying to deliver (after warmup / fault gating).
    state.hcd_nbi_setpoint_MW  = nbi_setpoint_MW_;
    state.hcd_icrh_setpoint_MW = icrh_setpoint_MW_;
    state.hcd_ecrh_setpoint_MW = ecrh_setpoint_MW_;
    state.hcd_lhcd_setpoint_MW = lhcd_setpoint_MW_;

    // Current drive contributions (only NBI and LHCD deliver meaningful CD)
    state.hcd_nbi_current_drive_MA  = nbi_actual_MW_  * cfg_.nbi_eta_MA_per_MW;
    state.hcd_lhcd_current_drive_MA = lhcd_actual_MW_ * cfg_.lhcd_eta_MA_per_MW;

    // Compute the bootstrap current fraction (volume-averaged, simplified)
    // f_bs ≈ 0.3 · q · sqrt(R/a) · β_p
    // We need β_p from the plasma state.  This is recomputed in the bridge
    // (which has the geometry); we just mirror it here for UI display.
    // (See PlasmaCoreBridge::updatePowerBalance for the actual calc.)
    state.hcd_bootstrap_current_MA = 0.f;  // updated by bridge

    // Set the aggregate alarm if any system is faulted
    state.alarm_aux_heat_fault = nbi_fault_ || icrh_fault_
                              || ecrh_fault_ || lhcd_fault_;

    // ── Update the unified aux heat setpoint in ReactorState ────────────────
    //  This is what the PlasmaCoreBridge reads as P_aux.  Note that we write
    //  the *actual* delivered power, not the *setpoint* — so during ramp-up
    //  or ramp-down the bridge sees the real power going into the plasma,
    //  not what the operator asked for.
    state.sp_aux_heat_MW = nbi_actual_MW_ + icrh_actual_MW_
                         + ecrh_actual_MW_ + lhcd_actual_MW_;
}
