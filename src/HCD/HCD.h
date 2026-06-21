#pragma once
//
// src/HCD/HCD.h
// Heating and Current Drive (H&CD) actuators.
//
//  Four independent systems, mirroring the ITER H&CD suite:
//
//    NBI  — Neutral Beam Injection
//           1 MeV deuterium neutrals, 16.5 MW (3 injectors × 5.5 MW).
//           Heats ions (50/50 ion/electron split at 1 MeV / 20 keV plasma).
//           Drives co-current non-inductive current (~1 MA per injector).
//           Cannot fire until the neutraliser is ready (~5 s after start).
//
//    ICRH — Ion Cyclotron Resonance Heating
//           40–55 MHz, 20 MW (8 antennas × 2.5 MW).
//           Minority heating (3He in D-T) — on-axis ion heating.
//           No net current drive (use the minority scheme).
//           Fast ramp (≈0.1 s); can be pulsed for MHD control.
//
//    ECRH — Electron Cyclotron Resonance Heating
//           170 GHz, 24 MW (12 gyrotrons × 2 MW).
//           Heats electrons locally at the resonance layer.
//           Used for startup assist (pre-ionization), NTM control, and
//           profile tailoring.  Fastest ramp of the four (≈10 ms).
//
//    LHCD — Lower Hybrid Current Drive
//           5 GHz, 20 MW (2 launchers × 10 MW).
//           Pure current drive — drives off-axis non-inductive current
//           (~2.5 MA per launcher at full power) for steady-state scenarios.
//           Slow ramp (≈1 s) due to klystron filament warm-up.
//
//  Each system has:
//    - on/off flag (operator command)
//    - setpoint power (MW) — what the operator asked for
//    - actual power (MW) — what's actually delivered, ramp-limited
//    - current drive contribution (MA) — for NBI and LHCD only
//    - fault flag (overheated gyrotron, beamline valve closed, etc.)
//
//  The ControlSystem writes a TOTAL aux-heat demand to state.sp_aux_heat_MW.
//  This module distributes it across the enabled systems in proportion to
//  their setpoints (capped by per-system max).  The PlasmaCoreBridge reads
//  the per-system actual powers and sums them into P_aux for the power
//  balance — this replaces the old hardcoded `P_aux = 50 MW`.
//

#include "ReactorState.h"
#include "SimTime.h"
#include <string>
#include <algorithm>

struct HCDConfig {
    // Per-system maximum power
    float nbi_max_MW   = 16.5f;   // 3 injectors × 5.5 MW
    float icrh_max_MW  = 20.0f;   // 8 antennas × 2.5 MW
    float ecrh_max_MW  = 24.0f;   // 12 gyrotrons × 2 MW
    float lhcd_max_MW  = 20.0f;   // 2 launchers × 10 MW

    // Ramp rates (MW/s) — system-dependent
    float nbi_ramp_MW_s   = 5.0f;   // slow: neutraliser + accelerator column
    float icrh_ramp_MW_s  = 50.0f;  // moderate: antenna matching network
    float ecrh_ramp_MW_s  = 200.0f; // fast: solid-state modulator
    float lhcd_ramp_MW_s  = 10.0f;  // slow: klystron filament

    // Neutraliser / filament warm-up times (s)
    float nbi_warmup_s   = 5.0f;   // beamline + neutraliser gas fill
    float lhcd_warmup_s  = 30.0f;  // klystron filament warm-up

    // Current drive efficiencies (MA per MW)
    //  NBI:  ≈ 0.06 MA/MW at 1 MeV (ITER-class)
    //  LHCD: ≈ 0.13 MA/MW at 5 GHz  (off-axis, low collisionality)
    //  ICRH/ECRH: small/negligible for the schemes modelled here
    float nbi_eta_MA_per_MW = 0.06f;
    float lhcd_eta_MA_per_MW= 0.13f;

    // Fault thresholds
    float nbi_beamline_max_temp_K = 500.f;   // beamline component overtemp
    float icrh_antenna_max_Vswr   = 3.0f;    // antenna mismatch
    float ecrh_gyrotron_max_temp_K= 350.f;   // gyrotron collector temp
    float lhcd_klystron_max_V_kV  = 90.f;    // klystron HV limit
};

class HCDSystem {
public:
    explicit HCDSystem(const HCDConfig& cfg = {});

    // Updates per-system actual powers, current drives, and fault flags.
    // Reads state.sp_aux_heat_MW (from ControlSystem) and the per-system
    // on/off + setpoint fields.  Writes the per-system actual fields back
    // to ReactorState for the PlasmaCoreBridge and UI to consume.
    void update(ReactorState& state, const SimTime& t);

    // Cold-restart: zero all powers, clear all warmup timers and faults.
    // Per-system on/off and setpoints are also zeroed (an H&CD system
    // doesn't survive a cold restart — operator must re-arm each one).
    void reset();

    // Operator commands
    void enableNBI (bool on)  { nbi_enabled_  = on; if (!on) nbi_setpoint_MW_ = 0.f; }
    void enableICRH(bool on)  { icrh_enabled_ = on; if (!on) icrh_setpoint_MW_ = 0.f; }
    void enableECRH(bool on)  { ecrh_enabled_ = on; if (!on) ecrh_setpoint_MW_ = 0.f; }
    void enableLHCD(bool on)  { lhcd_enabled_ = on; if (!on) lhcd_setpoint_MW_ = 0.f; }

    void setNBI (float mw) { nbi_setpoint_MW_  = std::clamp(mw, 0.f, cfg_.nbi_max_MW); }
    void setICRH(float mw) { icrh_setpoint_MW_ = std::clamp(mw, 0.f, cfg_.icrh_max_MW); }
    void setECRH(float mw) { ecrh_setpoint_MW_ = std::clamp(mw, 0.f, cfg_.ecrh_max_MW); }
    void setLHCD(float mw) { lhcd_setpoint_MW_ = std::clamp(mw, 0.f, cfg_.lhcd_max_MW); }

    // Force a fault on a system (for UI "inject fault" buttons)
    void faultNBI (const std::string& reason);
    void faultICRH(const std::string& reason);
    void faultECRH(const std::string& reason);
    void faultLHCD(const std::string& reason);

    // Clear latched faults (operator "reset fault" button)
    void clearFaults();

    // Diagnostic accessors for the UI
    bool  nbiReady()  const { return nbi_warmup_remaining_s_ <= 0.f && !nbi_fault_; }
    bool  lhcdReady() const { return lhcd_warmup_remaining_s_ <= 0.f && !lhcd_fault_; }
    float nbiWarmupRemaining()  const { return nbi_warmup_remaining_s_; }
    float lhcdWarmupRemaining() const { return lhcd_warmup_remaining_s_; }
    bool  nbiFault()  const { return nbi_fault_; }
    bool  icrhFault() const { return icrh_fault_; }
    bool  ecrhFault() const { return ecrh_fault_; }
    bool  lhcdFault() const { return lhcd_fault_; }
    const std::string& nbiFaultReason()  const { return nbi_fault_reason_; }
    const std::string& icrhFaultReason() const { return icrh_fault_reason_; }
    const std::string& ecrhFaultReason() const { return ecrh_fault_reason_; }
    const std::string& lhcdFaultReason() const { return lhcd_fault_reason_; }

private:
    HCDConfig cfg_;

    // Per-system actual powers (ramp-limited)
    float nbi_actual_MW_   = 0.f;
    float icrh_actual_MW_  = 0.f;
    float ecrh_actual_MW_  = 0.f;
    float lhcd_actual_MW_  = 0.f;

    // Per-system operator settings (mirrored to ReactorState on update)
    bool  nbi_enabled_  = false;
    bool  icrh_enabled_ = false;
    bool  ecrh_enabled_ = false;
    bool  lhcd_enabled_ = false;
    float nbi_setpoint_MW_  = 0.f;
    float icrh_setpoint_MW_ = 0.f;
    float ecrh_setpoint_MW_ = 0.f;
    float lhcd_setpoint_MW_ = 0.f;

    // Warmup timers (count down to zero; system can't deliver power until 0)
    float nbi_warmup_remaining_s_  = 0.f;
    float lhcd_warmup_remaining_s_ = 0.f;
    // Warmup-started flags: set true when warmup is triggered, cleared when
    // the system is disabled.  Prevents re-triggering warmup every tick while
    // the system is enabled but hasn't delivered power yet.
    bool  nbi_warmup_started_  = false;
    bool  lhcd_warmup_started_ = false;

    // Latched fault flags + reason strings
    bool  nbi_fault_  = false;
    bool  icrh_fault_ = false;
    bool  ecrh_fault_ = false;
    bool  lhcd_fault_ = false;
    std::string nbi_fault_reason_;
    std::string icrh_fault_reason_;
    std::string ecrh_fault_reason_;
    std::string lhcd_fault_reason_;

    // Helper: ramp actual toward setpoint at given rate
    static float ramp(float actual, float setpoint, float rate_MW_s, float dt);

    // Helper: distribute the total aux heat demand across enabled systems
    void distributeDemand(ReactorState& state);
};
