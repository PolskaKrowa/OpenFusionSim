#pragma once
//
// include/ReactorState.h
// Shared reactor state — summary fields only.
// Detailed per-turbine, per-pump, per-circuit state lives in each module's own header.
//
//  REVISION HISTORY (Round 4)
//  ───────────────────────────
//  Added H&CD actuator state (NBI, ICRH, ECRH, LHCD) so the operator can
//  drive the plasma with the same four actuators a real tokamak plant has,
//  rather than a single hardcoded "P_aux = 50 MW" magic number inside the
//  plasma core.  The bridge now reads sp_aux_heat_MW (sum of the four H&CD
//  systems' current outputs) instead of 50 MW.
//
//  Added vacuum-vessel state so the operator must sequence the roughing →
//  turbo pumps and reach a base vacuum before plasma initiation is allowed
//  (this was the third user request: "feel of operating the whole reactor").
//
//  Added plasma-shape state (kappa, delta, R_out, Z_x) so the shape
//  controller has somewhere to publish its outputs.
//
//  Added disruption-mitigation state (MGI / SPI) so the operator can
//  manually trigger a mitigation injection when a disruption is detected,
//  rather than just watching it happen.
//
//  Added sp_aux_heat_MW — the unified H&CD setpoint the ControlSystem
//  writes and the PlasmaCoreBridge reads.
//

#include <cstdint>

enum class PlasmaStatus : uint8_t {
    Cold = 0, Initiating = 1, Burning = 2, Disrupting = 3, Quenched = 4,
};
enum class ReactorMode : uint8_t {
    Startup = 0, SteadyState = 1, Rampdown = 2, Emergency = 3,
};

// ─── Disruption cause codes (used by the alarm system to give the operator
//     a meaningful "cause" rather than a generic "disruption risk" string).
//     Each code maps to a specific physics trigger and a recommended action.
enum class DisruptionCause : uint8_t {
    None              = 0,
    LowQ              = 1,  // q_95 below safety limit (kink mode risk)
    Greenwald         = 2,  // density above Greenwald limit (radiative collapse)
    Troyon            = 3,  // beta_N above Troyon limit (ideal MHD kink)
    TearingMode       = 4,  // m/n = 2/1 or 3/2 island grew past W_crit
    VDE               = 5,  // vertical displacement event
    HaloCurrent       = 6,  // post-disruption halo current detected
    RunawayElectrons  = 7,  // Dreicer field exceeded (runaway seed)
    DensityLimit      = 8,  // Murphy density limit (radiative collapse)
    DisruptionImminent= 9,  // MHD module flagged disruption_ongoing
};

struct ReactorState {

    // ── Time ──────────────────────────────────────────────────────────────────
    double  time_s  = 0.0;
    float   dt_s    = 1e-3f;
    int     tick    = 0;
    ReactorMode mode = ReactorMode::Startup;

    // ── Plasma ────────────────────────────────────────────────────────────────
    PlasmaStatus plasma_status   = PlasmaStatus::Cold;
    float plasma_current_MA      = 0.f;
    float plasma_temp_keV        = 0.f;
    float electron_temp_keV      = 0.f;
    float plasma_density_m3      = 0.f;
    float beta                   = 0.f;
    float beta_N                 = 0.f;   // normalized beta (Troyon)
    float q_safety               = 0.f;
    float fusion_power_MW        = 0.f;
    float alpha_power_MW         = 0.f;
    float neutron_flux_m2s       = 0.f;
    float radiated_power_MW      = 0.f;
    float tau_E_s                = 0.f;   // energy confinement time (IPB98(y,2))
    bool  disruption_flag        = false;
    DisruptionCause disruption_cause = DisruptionCause::None;

    // ── Plasma shape (output of the shape controller, input to plasma physics)
    //  These give the operator the same "shape knobs" a real tokamak operator
    //  has — the plasma cross-section is parameterised by elongation κ,
    //  triangularity δ, outer major radius R_out and X-point height Z_x.
    //  ITER targets: κ ≈ 1.7, δ ≈ 0.33, R_out ≈ 6.2 m, Z_x ≈ -1.0 m.
    float kappa                   = 1.7f;  // elongation
    float delta                   = 0.33f; // triangularity
    float R_out_m                 = 6.2f;  // outer major radius
    float Z_x_m                   = -1.0f; // X-point vertical position [m]
    float plasma_volume_m3        = 840.f; // computed volume

    // ── Magnets ───────────────────────────────────────────────────────────────
    float B_toroidal_T           = 0.f;
    float B_poloidal_T           = 0.f;
    float coil_current_kA        = 0.f;
    float magnet_temp_K          = 4.5f;
    bool  quench_detected        = false;
    float stored_energy_GJ       = 0.f;

    // ── Central Solenoid (CS) and PF coil summary ─────────────────────────────
    //  ITER-class: 6 CS modules + 6 PF coils.  The CS drives the loop voltage
    //  via transformer action; the PF coils shape and position the plasma.
    //  Both are operated independently from the TF coils.
    float cs_current_kA          = 0.f;  // CS total current (sum of 6 modules)
    float cs_loop_voltage_V      = 0.f;  // induced loop voltage
    float pf_current_kA          = 0.f;  // PF total current (sum of 6 coils)
    float cs_flux_remaining_Wb   = 0.f;  // remaining swing in CS (Vs); ITER: ~90 Wb

    // ── Fuel ──────────────────────────────────────────────────────────────────
    float fuel_D_inventory_g     = 0.f;
    float fuel_T_inventory_g     = 0.f;
    float fuel_injection_rate    = 0.f;
    float pellet_frequency_Hz    = 0.f;
    float D_T_ratio              = 1.f;

    // ── Plasma exhaust (helium ash) ───────────────────────────────────────────
    float helium_fraction        = 0.f;
    float pump_throughput_Pa     = 0.f;
    float divertor_power_MW      = 0.f;
    float divertor_temp_K        = 300.f;
    bool  divertor_overtemp      = false;

    // ── In-vessel particle inventory (what's actually IN the plasma right now) ─
    float vessel_D_g             = 0.f;
    float vessel_T_g             = 0.f;
    float vessel_He_g            = 0.f;
    float vessel_impurity_g      = 0.f;
    float impurity_fraction      = 0.f;
    float exhaust_pumping_Ls     = 0.f;

    // ── Thermal hydraulics (blanket primary circuit) ──────────────────────────
    float coolant_inlet_temp_K   = 573.f;
    float coolant_outlet_temp_K  = 773.f;
    float coolant_flow_kg_s      = 0.f;
    float blanket_heat_MW        = 0.f;
    float tbr_current            = 0.f;
    float first_wall_temp_K      = 300.f;
    bool  thermal_runaway        = false;

    // ── Electrical summary (written by TurbineSystem + ElectricalGrid) ────────
    float gross_electric_MW      = 0.f;
    float net_electric_MW        = 0.f;
    float parasitic_load_MW      = 0.f;
    float Q_scientific           = 0.f;
    float grid_frequency_Hz      = 50.f;
    bool  grid_connected         = false;

    // Legacy single-turbine fields (kept for PlasmaCoreBridge 0D model)
    float turbine_rpm            = 0.f;
    float steam_pressure_MPa     = 0.f;
    float steam_temp_K           = 0.f;
    bool  turbine_trip           = false;

    // ── Molten salt summary ───────────────────────────────────────────────────
    float hot_tank_temp_K        = 850.f;
    float cold_tank_temp_K       = 600.f;
    float hot_tank_level_m       = 8.f;
    float cold_tank_level_m      = 8.f;
    float salt_flow_total_kg_s   = 0.f;

    // ── Helium system summary ─────────────────────────────────────────────────
    float cryostat_temp_K        = 4.5f;
    bool  cryo_ok                = true;
    float reactor_he_outlet_K    = 800.f;
    float magnet_he_temp_K       = 4.5f;
    bool  cryo_compressor_on     = false;

    // ── Control setpoints ─────────────────────────────────────────────────────
    float sp_plasma_current_MA   = 15.f;
    float sp_electron_temp_keV   = 20.f;
    float sp_density_m3          = 1e20f;
    float sp_fuel_rate           = 1.f;
    float sp_B_toroidal_T        = 5.3f;
    float sp_coolant_flow        = 1.f;
    float sp_aux_heat_MW         = 0.f;  // TOTAL auxiliary heating (sum of H&CD)
    float sp_kappa               = 1.7f; // shape setpoints
    float sp_delta               = 0.33f;
    bool  cmd_scram              = false;
    bool  cmd_disrupt_mitigation = false; // operator-triggered MGI/SPI

    // ── H&CD (Heating and Current Drive) actuators ────────────────────────────
    //  Four independent systems, each with its own on/off and power setpoint.
    //  The plasma power balance sums all four into P_aux.  This replaces the
    //  old hardcoded `P_aux_MW = 50.0f` in PlasmaCoreBridge.
    //
    //  NBI  — Neutral Beam Injection: 1 MeV, 16.5 MW (ITER-class).
    //         Heats ions; also drives co-current (non-inductive current).
    //  ICRH — Ion Cyclotron Resonance Heating: 40-55 MHz, 20 MW.
    //         Heats ions at harmonics of cyclotron freq; minority heating.
    //  ECRH — Electron Cyclotron Resonance Heating: 170 GHz, 24 MW.
    //         Heats electrons; used for startup assist + MHD control.
    //  LHCD — Lower Hybrid Current Drive: 5 GHz, 20 MW.
    //         Drives off-axis non-inductive current for advanced scenarios.
    bool  hcd_nbi_on             = false;
    float hcd_nbi_setpoint_MW    = 0.f;  // commanded power
    float hcd_nbi_actual_MW      = 0.f;  // ramped/limited actual power
    bool  hcd_icrh_on            = false;
    float hcd_icrh_setpoint_MW   = 0.f;
    float hcd_icrh_actual_MW     = 0.f;
    bool  hcd_ecrh_on            = false;
    float hcd_ecrh_setpoint_MW   = 0.f;
    float hcd_ecrh_actual_MW     = 0.f;
    bool  hcd_lhcd_on            = false;
    float hcd_lhcd_setpoint_MW   = 0.f;
    float hcd_lhcd_actual_MW     = 0.f;
    float hcd_nbi_current_drive_MA = 0.f;  // non-inductive current from NBI
    float hcd_lhcd_current_drive_MA= 0.f;  // non-inductive current from LHCD
    float hcd_bootstrap_current_MA= 0.f;  // bootstrap current fraction

    // ── Vacuum vessel pumping ─────────────────────────────────────────────────
    //  Vessel must be pumped down from atmospheric (~101325 Pa) to <1e-4 Pa
    //  before plasma initiation.  Sequence: roughing pump (101325 → 1 Pa),
    //  then turbo pump (1 → 1e-5 Pa).  A pre-plasma-initiation interlock
    //  blocks the INITIATE button until vessel_pressure_Pa < 1e-3 Pa.
    float vessel_pressure_Pa     = 101325.f;  // start at atmospheric
    bool  vessel_roughing_on     = false;
    bool  vessel_turbo_on        = false;
    bool  vessel_bakeout_on      = false;    // wall bakeout (423 K) for outgassing
    float vessel_wall_temp_K     = 300.f;    // vessel wall temperature
    bool  vessel_vacuum_ok       = false;    // true when pressure < 1e-3 Pa
    float vessel_total_outgas_Pa = 0.f;      // accumulated outgas load

    // ── Disruption mitigation ─────────────────────────────────────────────────
    //  Two independent systems the operator can fire when a disruption is
    //  detected or imminent:
    //    MGI — Massive Gas Injection: fast valve injects Ne/Ar, radiatively
    //          quenches the plasma on a ~1 ms timescale.  Reduces halo current.
    //    SPI — Shattered Pellet Injection: deuterium pellet shattered against
    //          a shatter ring, injects a spray of solid fragments.  More
    //          effective at high current (ITER-class).
    bool  mgi_armed              = false;   // operator-armed (ready to fire)
    bool  mgi_fired              = false;   // latched once fired
    float mgi_pressure_Pa        = 0.f;     // gas pressure in vessel after MGI
    bool  spi_armed              = false;
    bool  spi_fired              = false;
    float spi_pellet_mass_g      = 0.f;     // mass of last injected pellet
    float mitigation_force_N     = 0.f;     // peak halo force reduction

    // ── Wall conditioning ─────────────────────────────────────────────────────
    //  Periodic boronization coats the first wall with a thin boron film that
    //  getter oxygen and reduces impurity radiation.  After many discharges
    //  the film erodes; the operator can request a fresh boronization.
    float boronization_thickness_nm = 0.f;  // 0 = no film; ~200 nm fresh coat
    float last_boronization_time_s  = -1e9f; // when last applied (relative to sim time)
    bool  boronization_in_progress  = false;

    // ── Tritium plant ─────────────────────────────────────────────────────────
    //  Tritium extracted from the blanket isotope-separation column; periodic
    //  accountancy tracks how much T is in the plant vs. in the fuel store.
    float tritium_in_plant_g     = 0.f;     // inventory in ISS/SDU
    float tritium_recovery_rate_g_s = 0.f;  // from TES (Tritium Extraction System)
    bool  tritium_recovery_on    = false;
    float detritiation_flow_m3_s = 0.f;     // detritiation airflow rate

    // ── Alarms ────────────────────────────────────────────────────────────────
    bool  alarm_disruption       = false;
    bool  alarm_quench           = false;
    bool  alarm_overtemp         = false;
    bool  alarm_loss_of_coolant  = false;
    bool  alarm_low_tritium      = false;
    bool  alarm_vacuum_loss      = false;   // vessel vacuum degraded
    bool  alarm_runaway          = false;   // runaway electron risk
    bool  alarm_halo             = false;   // halo current > design
    bool  alarm_aux_heat_fault   = false;   // H&CD system tripped
    bool  alarm_vessel_breach    = false;   // vacuum vessel breach
};
