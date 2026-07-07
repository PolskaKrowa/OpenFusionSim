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

    // ── Two-temperature plasma + confinement regime (realism round 5) ─────────
    //  The plasma now carries separate electron and ion temperatures coupled
    //  by Coulomb collisional equilibration.  plasma_temp_keV keeps its old
    //  meaning of "ion temperature" (the OVERVIEW tab labels it Ti) and
    //  ion_temp_keV is an explicit alias for new code.
    float ion_temp_keV           = 0.f;   // T_i [keV] (== plasma_temp_keV)
    bool  h_mode                 = false; // true when confinement is H-mode
    float P_LH_MW                = 0.f;   // Martin-2008 L→H power threshold [MW]
    float H98_factor             = 1.f;   // confinement enhancement in use
    float greenwald_frac         = 0.f;   // n_e / n_Greenwald (1.0 = at limit)

    // ── Power-balance breakdown (all MW; written by PlasmaCoreBridge) ─────────
    //  Published so the UI can draw the full source/sink ledger instead of
    //  a single lumped "radiated power" number.
    float P_ohm_MW               = 0.f;   // ohmic heating (I_p² · R_plasma)
    float P_aux_actual_MW        = 0.f;   // H&CD power actually delivered
    float P_brem_MW              = 0.f;   // bremsstrahlung loss
    float P_sync_MW              = 0.f;   // synchrotron loss
    float P_line_MW              = 0.f;   // impurity line-radiation loss
    float P_cond_MW              = 0.f;   // conduction/convection loss (W/τ_E)
    float dW_dt_MW               = 0.f;   // net d(stored energy)/dt
    float loop_voltage_V         = 0.f;   // resistive loop voltage V = I_p·R_p

    // ── Sawtooth activity (q₀ < 1 internal-kink relaxation cycle) ─────────────
    float q0_estimate            = 3.f;   // on-axis safety factor estimate
    int   sawtooth_count         = 0;     // crashes this discharge
    float sawtooth_period_s      = 0.f;   // current crash period (0 = inactive)

    // ── ELM activity (Type-I edge-localized modes; H-mode only) ───────────────
    //  Every H-mode plasma builds an edge pedestal that periodically crashes,
    //  expelling a few percent of the stored energy onto the divertor in a
    //  sub-ms burst.  Unmitigated ITER-scale Type-I ELMs deliver tens of MJ
    //  per crash — enough to melt tungsten tiles.  Pellet pacing (firing the
    //  fuelling pellets faster than the natural ELM frequency) triggers more
    //  frequent, proportionally smaller ELMs: ΔW_ELM · f_ELM ≈ const.
    float elm_freq_Hz            = 0.f;   // current ELM frequency (0 = no ELMs)
    float elm_size_MJ            = 0.f;   // energy expelled by the last crash
    int   elm_count              = 0;     // ELM crashes this discharge
    bool  elm_paced              = false; // pellet pacing currently active
    //  Pending divertor energy pulse [MJ] — written by the plasma core on
    //  each ELM crash, consumed (zeroed) by the Helium divertor thermal model.
    float elm_div_pulse_MJ       = 0.f;

    // ── Runaway electron beam (disruption consequence) ────────────────────────
    //  During a fast, UNMITIGATED current quench the toroidal E-field spike
    //  avalanche-multiplies any seed electrons into a multi-MA relativistic
    //  beam.  When the quench completes, the beam loses position control and
    //  strikes the first wall, depositing its magnetic energy (hundreds of
    //  MJ at ITER scale) in a localized spot.  MGI/SPI mitigation raises the
    //  density and collisionally suppresses the avalanche — this is one of
    //  the two reasons the mitigation systems exist (the other is halo loads).
    float re_beam_MA             = 0.f;   // runaway beam current
    bool  re_wall_strike         = false; // latched: beam hit the wall
    //  Pending localized first-wall energy dump [MJ] — consumed by
    //  ThermalHydraulics::updateFirstWall.
    float re_wall_energy_MJ      = 0.f;
    int   re_strike_count        = 0;     // wall strikes this campaign

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
    //  CS volt-second budget.  The central solenoid is a transformer with a
    //  FINITE flux swing (ITER: ~90 Wb available for resistive flux
    //  consumption during the burn).  The plasma consumes flux at
    //  dΨ/dt = V_loop × (inductive fraction of I_p); when the swing is
    //  exhausted, the inductive drive is gone — I_p decays toward whatever
    //  the non-inductive systems (NBI/LHCD current drive + bootstrap) can
    //  carry.  This is THE pulse-length limit of a conventional tokamak and
    //  the entire reason steady-state scenarios chase non-inductive current.
    //  Owned by PlasmaCoreBridge (it knows V_loop and the CD currents).
    float cs_flux_remaining_Wb   = 90.f; // remaining swing in CS (Vs)
    bool  cs_flux_exhausted      = false;// no inductive drive left

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
    //  Neutron damage bookkeeping.  14 MeV neutrons displace atoms in the
    //  first-wall steel; damage is measured in dpa (displacements per atom).
    //  Rule of thumb for RAFM steel behind a W armour: ~10 dpa per MW·a/m²
    //  of neutron wall load.  ITER first wall is designed for ~1-3 dpa
    //  lifetime; DEMO needs materials good to 50+.  Accumulates over the
    //  whole campaign (not reset per discharge — only on cold restart).
    float neutron_wall_load_MW_m2 = 0.f;  // instantaneous neutron wall load
    float fw_fluence_MWa_m2       = 0.f;  // integrated wall load [MW·a/m²]
    float fw_dpa                  = 0.f;  // accumulated displacement damage

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
    //  NBI shine-through interlock: 1 MeV neutrals fired into a thin plasma
    //  are not ionized before they cross it — they slam into the far wall at
    //  full power density.  Real machines interlock the beams on a minimum
    //  line density; below n_e ≈ 1.5×10¹⁹ m⁻³ the HCD module gates NBI to 0.
    bool  nbi_shinethrough_block = false;

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
    //  Auto-MGI enable — lives in ReactorState so the UI checkbox persists
    //  (it used to be a local variable that reset to false every frame,
    //  making auto-MGI impossible to actually enable).
    bool  dm_auto_mgi            = false;
    //  Radiator inventory delivered by MGI/SPI (fractional impurity of Ne/Ar
    //  equivalents).  Consumed by the PlasmaCoreBridge line-radiation model:
    //  this is HOW mitigation kills the plasma — a radiative thermal quench
    //  through the actual radiation physics, not a scripted status flip.
    float dm_injected_impurity   = 0.f;
    //  Latched when a quench proceeds under mitigation control: the current
    //  quench is slower/gentler (τ_CQ ≈ 100 ms vs ~30 ms unmitigated),
    //  which is the entire point of firing MGI/SPI.
    bool  dm_mitigated           = false;
    // Pressure rise requested by DM module (MGI gas injection).  VacuumVessel
    // reads this and applies it to its internal pressure_Pa_ via
    // forcePressureRise().  Without this indirection, DM writing
    // state.vessel_pressure_Pa directly was overwritten by VacuumVessel::update
    // on the next tick.
    float dm_pressure_rise_Pa    = 0.f;

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
