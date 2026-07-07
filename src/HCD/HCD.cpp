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
    nbi_warmup_started_  = false;
    lhcd_warmup_started_ = false;
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
    //
    //  Two command paths feed into the H&CD systems:
    //
    //  1. OPERATOR path: the UI sliders/presets write per-system setpoints
    //     directly to state.hcd_*_setpoint_MW.  These are the operator's
    //     "I want this much power from NBI/ICRH/ECRH/LHCD" commands.
    //
    //  2. CONTROLLER path: the ControlSystem's temperature PID writes a
    //     *total* aux-heat demand to state.sp_aux_heat_MW.  This is the
    //     "I need X MW of total heating to reach my T_e setpoint" command.
    //
    //  The two paths are reconciled here: each enabled system's setpoint is
    //  the MINIMUM of (operator's per-system setpoint, PID-allocated share).
    //  The PID-allocated share = sp_aux_heat_MW × (system's max / total max
    //  of enabled systems).  This means:
    //
    //    - If the PID demands less than the operator asked for, the PID wins
    //      (the plasma is hot enough — we don't need full operator power).
    //    - If the PID demands more than the operator asked for, the operator
    //      wins (the PID can't force a system past its setpoint cap).
    //    - If the PID is at zero (T_e above setpoint), all systems ramp down.
    //
    //  This closes the temperature control loop: T_e error → PID output →
    //  sp_aux_heat_MW → per-system setpoints → actual MW → P_aux in the
    //  power balance → T_e.  Previously sp_aux_heat_MW was written but never
    //  read, so the loop was open and T_e was uncontrolled.

    // Sync operator enable/disable commands from ReactorState
    nbi_enabled_  = state.hcd_nbi_on;
    icrh_enabled_ = state.hcd_icrh_on;
    ecrh_enabled_ = state.hcd_ecrh_on;
    lhcd_enabled_ = state.hcd_lhcd_on;

    // ── NBI/LHCD warmup: trigger on rising edge of enable ──────────────────
    //  When the operator enables NBI (or LHCD) and the warmup timer is at
    //  zero AND no actual power has been delivered yet, start the warmup
    //  countdown.  This makes the "5 s NBI neutraliser" and "30 s LHCD
    //  klystron filament" warmup times actually take effect.
    if (nbi_enabled_ && nbi_warmup_remaining_s_ <= 0.f && nbi_actual_MW_ < 0.01f
        && !nbi_warmup_started_) {
        nbi_warmup_remaining_s_ = cfg_.nbi_warmup_s;
        nbi_warmup_started_ = true;
    }
    if (lhcd_enabled_ && lhcd_warmup_remaining_s_ <= 0.f && lhcd_actual_MW_ < 0.01f
        && !lhcd_warmup_started_) {
        lhcd_warmup_remaining_s_ = cfg_.lhcd_warmup_s;
        lhcd_warmup_started_ = true;
    }
    // Reset warmup-started flag when system is disabled, so re-enabling
    // triggers a fresh warmup cycle.
    if (!nbi_enabled_)  nbi_warmup_started_  = false;
    if (!lhcd_enabled_) lhcd_warmup_started_ = false;

    // Read operator per-system setpoints (clamped to per-system max)
    float op_nbi  = std::clamp(state.hcd_nbi_setpoint_MW,  0.f, cfg_.nbi_max_MW);
    float op_icrh = std::clamp(state.hcd_icrh_setpoint_MW, 0.f, cfg_.icrh_max_MW);
    float op_ecrh = std::clamp(state.hcd_ecrh_setpoint_MW, 0.f, cfg_.ecrh_max_MW);
    float op_lhcd = std::clamp(state.hcd_lhcd_setpoint_MW, 0.f, cfg_.lhcd_max_MW);

    // ── PID allocation: distribute sp_aux_heat_MW proportionally ───────────
    //  Each enabled system gets a share proportional to its max capacity.
    //  E.g. if NBI (16.5) and ECRH (24) are enabled and PID demands 20 MW,
    //  NBI gets 20×16.5/40.5 = 8.1 MW, ECRH gets 20×24/40.5 = 11.9 MW.
    float total_max = 0.f;
    if (nbi_enabled_  && !nbi_fault_)  total_max += cfg_.nbi_max_MW;
    if (icrh_enabled_ && !icrh_fault_) total_max += cfg_.icrh_max_MW;
    if (ecrh_enabled_ && !ecrh_fault_) total_max += cfg_.ecrh_max_MW;
    if (lhcd_enabled_ && !lhcd_fault_) total_max += cfg_.lhcd_max_MW;

    float pid_demand = std::max(0.f, state.sp_aux_heat_MW);

    if (total_max > 0.1f) {
        if (nbi_enabled_  && !nbi_fault_)
            nbi_setpoint_MW_  = std::min(op_nbi,  pid_demand * cfg_.nbi_max_MW  / total_max);
        else                nbi_setpoint_MW_  = 0.f;
        if (icrh_enabled_ && !icrh_fault_)
            icrh_setpoint_MW_ = std::min(op_icrh, pid_demand * cfg_.icrh_max_MW / total_max);
        else                icrh_setpoint_MW_ = 0.f;
        if (ecrh_enabled_ && !ecrh_fault_)
            ecrh_setpoint_MW_ = std::min(op_ecrh, pid_demand * cfg_.ecrh_max_MW / total_max);
        else                ecrh_setpoint_MW_ = 0.f;
        if (lhcd_enabled_ && !lhcd_fault_)
            lhcd_setpoint_MW_ = std::min(op_lhcd, pid_demand * cfg_.lhcd_max_MW / total_max);
        else                lhcd_setpoint_MW_ = 0.f;
    } else {
        // No systems enabled — all setpoints zero
        nbi_setpoint_MW_  = 0.f;
        icrh_setpoint_MW_ = 0.f;
        ecrh_setpoint_MW_ = 0.f;
        lhcd_setpoint_MW_ = 0.f;
    }

    // Faulted systems get 0 (redundant with the if-checks above, but explicit)
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

    // ── NBI warmup countdown ────────────────────────────────────────────────
    //  The warmup timer is started by distributeDemand() on the rising edge
    //  of the enable flag.  Here we just count it down and gate the
    //  effective setpoint to 0 while warmup is in progress.
    if (nbi_warmup_remaining_s_ > 0.f) {
        nbi_warmup_remaining_s_ = std::max(0.f, nbi_warmup_remaining_s_ - dt);
    }
    float nbi_effective_setpoint = nbi_setpoint_MW_;
    if (nbi_warmup_remaining_s_ > 0.f) {
        nbi_effective_setpoint = 0.f;  // still warming up
    }
    if (nbi_fault_) nbi_effective_setpoint = 0.f;

    // ── NBI shine-through interlock ─────────────────────────────────────────
    //  A 1 MeV neutral beam is ionized by collisions with the plasma; at low
    //  line density the beam isn't stopped before it crosses the vessel and
    //  slams into the far wall at ~10 MW/m² (ITER's beams would melt their
    //  own beam dumps in seconds).  Real machines interlock the injectors on
    //  a minimum density — modelled here as a hard gate at 1.5×10¹⁹ m⁻³.
    //  Practical consequence for the operator: build density with gas puff /
    //  pellets (and heat with ECRH/ICRH/ohmic) BEFORE bringing on the beams.
    constexpr float NBI_MIN_DENSITY_M3 = 1.5e19f;
    bool shine_block = state.plasma_density_m3 < NBI_MIN_DENSITY_M3;
    state.nbi_shinethrough_block = shine_block;
    if (shine_block) nbi_effective_setpoint = 0.f;

    // ── LHCD warmup countdown ───────────────────────────────────────────────
    if (lhcd_warmup_remaining_s_ > 0.f) {
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

    // ── Do NOT mirror the reconciled setpoints back into ReactorState ───────
    //  BUG FIX: this used to write min(operator setpoint, PID share) back
    //  into state.hcd_*_setpoint_MW — the same fields the operator's sliders
    //  live in.  That made the reconciliation a RATCHET: any tick where the
    //  temperature PID demanded less than the operator (e.g. T_e briefly
    //  above its setpoint → PID output 0) permanently overwrote the
    //  operator's command with 0, and since the PID share is min'ed against
    //  the (now zero) operator value, the heating could never come back.
    //  The operator fields stay untouched; the reconciled values live in the
    //  module-internal *_setpoint_MW_ members and the delivered power is
    //  visible through hcd_*_actual_MW.

    // Current drive contributions (only NBI and LHCD deliver meaningful CD)
    state.hcd_nbi_current_drive_MA  = nbi_actual_MW_  * cfg_.nbi_eta_MA_per_MW;
    state.hcd_lhcd_current_drive_MA = lhcd_actual_MW_ * cfg_.lhcd_eta_MA_per_MW;

    // NOTE: state.hcd_bootstrap_current_MA is written by PlasmaCoreBridge
    // (it needs β_p and q_safety which only the bridge has).  Do NOT write
    // it here — doing so would overwrite the bridge's value with 0 every
    // tick, which was the previous bug.

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
