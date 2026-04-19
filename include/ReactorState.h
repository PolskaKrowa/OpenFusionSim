#pragma once
//
// include/ReactorState.h
// Shared reactor state — summary fields only.
// Detailed per-turbine, per-pump, per-circuit state lives in each module's own header.
//

#include <cstdint>

enum class PlasmaStatus : uint8_t {
    Cold = 0, Initiating = 1, Burning = 2, Disrupting = 3, Quenched = 4,
};
enum class ReactorMode : uint8_t {
    Startup = 0, SteadyState = 1, Rampdown = 2, Emergency = 3,
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
    float q_safety               = 0.f;
    float fusion_power_MW        = 0.f;
    float alpha_power_MW         = 0.f;
    float neutron_flux_m2s       = 0.f;
    float radiated_power_MW      = 0.f;
    bool  disruption_flag        = false;

    // ── Magnets ───────────────────────────────────────────────────────────────
    float B_toroidal_T           = 0.f;
    float B_poloidal_T           = 0.f;
    float coil_current_kA        = 0.f;
    float magnet_temp_K          = 4.5f;
    bool  quench_detected        = false;
    float stored_energy_GJ       = 0.f;

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

    // ── Thermal hydraulics (blanket primary circuit) ──────────────────────────
    float coolant_inlet_temp_K   = 573.f;
    float coolant_outlet_temp_K  = 773.f;
    float coolant_flow_kg_s      = 0.f;
    float blanket_heat_MW        = 0.f;
    float tbr_current            = 0.f;
    float first_wall_temp_K      = 300.f;
    bool  thermal_runaway        = false;

    // ── Electrical summary (written by TurbineSystem + ElectricalGrid) ────────
    float gross_electric_MW      = 0.f;   // sum over all turbines
    float net_electric_MW        = 0.f;   // gross minus on-site loads
    float parasitic_load_MW      = 0.f;
    float Q_scientific           = 0.f;
    float grid_frequency_Hz      = 50.f;
    bool  grid_connected         = false; // at least one generator on bus

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
    float reactor_he_outlet_K    = 800.f; // reactor structure He cooling outlet
    float magnet_he_temp_K       = 4.5f;
    bool  cryo_compressor_on     = false;

    // ── Control setpoints ─────────────────────────────────────────────────────
    float sp_plasma_current_MA   = 15.f;
    float sp_electron_temp_keV   = 20.f;
    float sp_density_m3          = 1e20f;
    float sp_fuel_rate           = 1.f;
    float sp_B_toroidal_T        = 5.3f;
    float sp_coolant_flow        = 1.f;
    bool  cmd_scram              = false;

    // ── Alarms ────────────────────────────────────────────────────────────────
    bool  alarm_disruption       = false;
    bool  alarm_quench           = false;
    bool  alarm_overtemp         = false;
    bool  alarm_loss_of_coolant  = false;
    bool  alarm_low_tritium      = false;
};