#pragma once

//
// include/ReactorState.h
// Shared reactor state struct passed by reference between all modules each tick.
//
// Modules READ their inputs from ReactorState and WRITE their outputs back
// into ReactorState.  The Control module is the only one allowed to write
// to the setpoint fields.  All other modules write to their output fields only.
//
// All physical quantities are in SI units unless a comment says otherwise.
//

#include <array>
#include <cstdint>

// ─── Plasma status ─────────────────────────────────────────────────────────────
enum class PlasmaStatus : uint8_t {
    Cold        = 0,    // no plasma
    Initiating  = 1,    // ohmic heating / startup
    Burning     = 2,    // self-sustaining D-T burn (Q > 1)
    Disrupting  = 3,    // disruption in progress
    Quenched    = 4,    // disruption handled / plasma lost
};

// ─── Reactor operational mode ─────────────────────────────────────────────────
enum class ReactorMode : uint8_t {
    Startup     = 0,
    SteadyState = 1,
    Rampdown    = 2,
    Emergency   = 3,
};

// ─── Full reactor state ────────────────────────────────────────────────────────
struct ReactorState {

    // ── Simulation time ──────────────────────────────────────────────────────
    double  time_s       = 0.0;     // wall-clock sim time [s]
    float   dt_s         = 1e-3f;   // current timestep [s]
    int     tick         = 0;       // integer tick counter
    ReactorMode mode     = ReactorMode::Startup;

    // ── Plasma (written by PlasmaCore) ───────────────────────────────────────
    PlasmaStatus plasma_status = PlasmaStatus::Cold;
    float   plasma_current_MA  = 0.0f;   // Ip [MA]
    float   plasma_temp_keV    = 0.0f;   // volume-average ion temperature [keV]
    float   electron_temp_keV  = 0.0f;   // volume-average electron temperature [keV]
    float   plasma_density_m3  = 0.0f;   // volume-average n_e [m^-3]
    float   beta               = 0.0f;   // plasma beta (pressure / magnetic pressure)
    float   q_safety           = 0.0f;   // safety factor q at 95% flux surface
    float   fusion_power_MW    = 0.0f;   // total fusion power [MW]
    float   alpha_power_MW     = 0.0f;   // alpha heating power [MW]
    float   neutron_flux_m2s   = 0.0f;   // first-wall neutron flux [n/m^2/s]
    float   radiated_power_MW  = 0.0f;   // total radiated power (Bremsstrahlung + impurities) [MW]
    bool    disruption_flag    = false;  // true when disruption is imminent

    // ── Magnets (written by Magnets) ─────────────────────────────────────────
    float   B_toroidal_T      = 0.0f;   // toroidal field on axis [T]
    float   B_poloidal_T      = 0.0f;   // poloidal field [T]
    float   coil_current_kA   = 0.0f;   // TF coil current [kA]
    float   magnet_temp_K     = 4.5f;   // superconducting coil temperature [K]
    bool    quench_detected   = false;  // SC quench flag
    float   stored_energy_GJ  = 0.0f;  // magnetic energy stored [GJ]

    // ── Fuel (written by Fuel) ───────────────────────────────────────────────
    float   fuel_D_inventory_g   = 0.0f;  // deuterium gas inventory [g]
    float   fuel_T_inventory_g   = 0.0f;  // tritium inventory [g]
    float   fuel_injection_rate  = 0.0f;  // current injection rate [Pa·m^3/s]
    float   pellet_frequency_Hz  = 0.0f;  // pellet injector frequency [Hz]
    float   D_T_ratio            = 1.0f;  // D:T fuelling ratio (1.0 = 50:50)

    // ── Helium (written by Helium) ────────────────────────────────────────────
    float   helium_fraction      = 0.0f;  // He ash fraction of plasma (dilutes burn)
    float   pump_throughput_Pa   = 0.0f;  // divertor pump throughput [Pa·m^3/s]
    float   divertor_power_MW    = 0.0f;  // heat load on divertor [MW]
    float   divertor_temp_K      = 300.f; // divertor tile temperature [K]
    bool    divertor_overtemp    = false; // tile over-temperature alarm

    // ── ThermalHydraulics (written by ThermalHydraulics) ──────────────────────
    float   coolant_inlet_temp_K   = 573.0f;  // primary coolant inlet [K]
    float   coolant_outlet_temp_K  = 773.0f;  // primary coolant outlet [K]
    float   coolant_flow_kg_s      = 0.0f;    // primary coolant mass flow [kg/s]
    float   blanket_heat_MW        = 0.0f;    // thermal power extracted by blanket [MW]
    float   tbr_current            = 0.0f;    // tritium breeding ratio (dimensionless)
    float   first_wall_temp_K      = 300.0f;  // peak first-wall temperature [K]
    bool    thermal_runaway        = false;   // blanket overheating flag

    // ── SteamPower (written by SteamPower) ───────────────────────────────────
    float   gross_electric_MW    = 0.0f;    // gross electrical output [MW]
    float   parasitic_load_MW    = 0.0f;    // recirculating power (magnets, pumps) [MW]
    float   net_electric_MW      = 0.0f;    // net to grid [MW]
    float   Q_scientific         = 0.0f;    // fusion gain Q = P_fusion / P_heating
    float   turbine_rpm          = 0.0f;    // turbine shaft speed [RPM]
    float   steam_pressure_MPa   = 0.0f;    // live steam pressure [MPa]
    float   steam_temp_K         = 0.0f;    // live steam temperature [K]
    bool    turbine_trip         = false;   // turbine protection trip

    // ── Control setpoints (written only by Control) ──────────────────────────
    float   sp_plasma_current_MA  = 15.0f;  // Ip setpoint [MA]  (ITER-class: 15 MA)
    float   sp_electron_temp_keV  = 20.0f;  // T_e setpoint [keV]
    float   sp_density_m3         = 1e20f;  // n_e setpoint [m^-3]
    float   sp_fuel_rate          = 1.0f;   // fuelling rate setpoint [normalised 0–1]
    float   sp_B_toroidal_T       = 5.3f;   // B_T setpoint [T]   (ITER: 5.3 T)
    float   sp_coolant_flow       = 1.0f;   // coolant flow setpoint [normalised 0–1]
    bool    cmd_scram             = false;  // emergency shutdown command

    // ── Alarms ────────────────────────────────────────────────────────────────
    bool    alarm_disruption      = false;
    bool    alarm_quench          = false;
    bool    alarm_overtemp        = false;
    bool    alarm_loss_of_coolant = false;
    bool    alarm_low_tritium     = false;
};