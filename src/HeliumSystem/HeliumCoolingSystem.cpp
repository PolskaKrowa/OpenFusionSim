//
// src/HeliumSystem/HeliumCoolingSystem.cpp
//
#include "HeliumCoolingSystem.h"
#include <cmath>
#include <algorithm>

static constexpr float CP_HE_WARM   = 5193.f;   // J/(kg·K) warm He
static constexpr float CP_HE_CRYO   = 5260.f;   // J/(kg·K) supercritical 4.5K
static constexpr float SIGMA        = 5.67e-8f;  // Stefan-Boltzmann

HeliumCoolingSystem::HeliumCoolingSystem()
{
    // Warm He circuit defaults
    auto& rc = s_.reactor_circuit;
    for (auto& p : rc.pump) {
        p.inlet_temp_K  = 700.f;
        p.outlet_temp_K = 800.f;
        p.pressure_MPa  = 8.f;
    }

    // Cryo circuit defaults
    auto& mc = s_.magnet_circuit;
    mc.refrigerator_on = false;
    for (auto& p : mc.cryo_pump) {
        p.inlet_temp_K  = 4.3f;
        p.outlet_temp_K = 4.7f;
        p.pressure_MPa  = 0.3f;
    }

    // Cryostat defaults
    s_.cryostat.roughing_pump_on = false;
    s_.cryostat.turbo_pump_on    = false;
}

// ─── Reactor warm He circuit ──────────────────────────────────────────────────
void HeliumCoolingSystem::updateReactorCircuit(float dt,
                                                float fusion_power_MW,
                                                float first_wall_temp_K)
{
    auto& rc = s_.reactor_circuit;

    float total_flow = 0.f;
    for (int i = 0; i < 2; i++) {
        auto& p = rc.pump[i];
        if (p.running && !p.trip) {
            p.flow_kg_s = p.speed_frac * 200.f; // each pump up to 200 kg/s
            p.power_MW  = p.flow_kg_s * 0.05f;  // ~50 kJ/kg compression
            total_flow += p.flow_kg_s;
        } else {
            p.flow_kg_s = 0.f;
            p.power_MW  = 0.f;
        }
    }

    if (total_flow < 1.f) {
        rc.lo_flow_alarm = true;
        return;
    }
    rc.lo_flow_alarm = false;

    // Heat absorbed from first wall structure (~3% of fusion power)
    float Q_fw_MW  = fusion_power_MW * 0.03f;
    // ΔT = Q / (ṁ · cp)
    float dT       = Q_fw_MW * 1e6f / (total_flow * CP_HE_WARM);
    rc.outlet_temp_K = rc.inlet_temp_K + dT;
    rc.outlet_temp_K = std::clamp(rc.outlet_temp_K, rc.inlet_temp_K, 1200.f);

    rc.heat_removed_MW     = Q_fw_MW;
    rc.aux_HX_duty_MW      = Q_fw_MW; // all goes to air coolers
    rc.hi_outlet_temp_alarm= (rc.outlet_temp_K > 900.f);

    // HX cools He back to inlet temp (air coolers, simple model)
    rc.inlet_temp_K = std::max(rc.outlet_temp_K - dT * 0.98f, 680.f);

    (void)first_wall_temp_K;
    (void)dt;
}

// ─── Cryostat vacuum ─────────────────────────────────────────────────────────
void HeliumCoolingSystem::updateCryostat(float dt)
{
    auto& cv = s_.cryostat;

    // Vacuum degrades slowly (virtual outgassing at 5e-9 Pa/s)
    cv.vacuum_pressure_Pa += 5e-9f * dt;

    // Roughing pump drives pressure down from atmosphere to ~1 Pa quickly
    if (cv.roughing_pump_on)
        cv.vacuum_pressure_Pa = std::max(cv.vacuum_pressure_Pa - 0.1f * dt, 1.f);

    // Turbo-molecular pump drives from 1 Pa to high-vacuum regime
    if (cv.turbo_pump_on && cv.vacuum_pressure_Pa < 10.f)
        cv.vacuum_pressure_Pa = std::max(cv.vacuum_pressure_Pa - 0.001f * dt, 1e-5f);

    cv.vacuum_ok = (cv.vacuum_pressure_Pa < 0.1f);

    // Static heat in-leak to cold mass
    // Radiation: Q = σ · ε_eff · A · (T_hot⁴ - T_cold⁴)
    float eps_eff = 0.001f; // MLI ~30 layers, ε_eff very low
    float T4_diff = std::pow(cv.outer_wall_K, 4.f) - std::pow(4.5f, 4.f);
    cv.heat_leak_W = SIGMA * eps_eff * cv.surface_area_m2 * T4_diff;
    // Add conduction through supports (~300 W typical)
    cv.heat_leak_W += 300.f;

    // Thermal shield maintained by cryo system
    cv.thermal_shield_K = cv.vacuum_ok ? 80.f : 150.f;
}

// ─── Magnet cryo circuit (4.5 K) ─────────────────────────────────────────────
void HeliumCoolingSystem::updateMagnetCircuit(float dt,
                                               float magnet_temp_K,
                                               float stored_energy_GJ)
{
    auto& mc = s_.magnet_circuit;

    // Heat loads on cold mass
    float Q_static_W  = s_.cryostat.heat_leak_W;             // cryostat in-leak
    float Q_eddy_W    = stored_energy_GJ * 50.f;              // eddy currents in coil structure
    float Q_joints_W  = 200.f;                                 // inter-coil joint resistances
    mc.total_heat_load_W = Q_static_W + Q_eddy_W + Q_joints_W;

    // Cryo refrigerator: removes heat from cold mass
    // COP at 4.5 K: roughly 1/250 (need 250 W electrical per W removed)
    float Q_removed_W  = 0.f;
    if (mc.refrigerator_on) {
        mc.cryo_refrigerator_MW = mc.total_heat_load_W / 1000.f * 250.f / 1e3f; // MW
        mc.cryo_refrigerator_MW = std::clamp(mc.cryo_refrigerator_MW, 0.f, 80.f);
        Q_removed_W = mc.total_heat_load_W;
    } else {
        mc.cryo_refrigerator_MW = 0.f;
    }

    // Cryo pumps circulate supercritical He
    for (int i = 0; i < 2; i++) {
        auto& p = mc.cryo_pump[i];
        if (p.running && !p.trip) {
            p.flow_kg_s = p.speed_frac * 5.f;   // each ~5 kg/s supercritical He
            p.power_MW  = p.flow_kg_s * 2e-3f;  // very low power (circulators)
        } else {
            p.flow_kg_s = 0.f;
        }
    }

    // Cold mass temperature dynamics
    // dT/dt = (Q_in - Q_removed) / C_cold_mass
    // C for all coils: ~500 kJ/K (large mass of Nb3Sn + structural steel + He)
    float C_cold = 5e5f; // J/K
    float dT_dt  = (mc.total_heat_load_W - Q_removed_W) / C_cold;
    mc.cold_mass_temp_K += dT_dt * dt;
    mc.cold_mass_temp_K  = std::clamp(mc.cold_mass_temp_K, 2.f, 50.f);

    // Temperature tracking
    mc.return_temp_K  = mc.cold_mass_temp_K + 0.2f;
    mc.supply_temp_K  = mc.cold_mass_temp_K - 0.2f;

    mc.lo_temp_alarm  = (mc.cold_mass_temp_K < 4.0f);
    mc.hi_temp_alarm  = (mc.cold_mass_temp_K > 5.5f);
    mc.lo_pressure_alarm = (mc.pressure_MPa < 0.15f);

    (void)magnet_temp_K;
    (void)stored_energy_GJ;
    (void)dt;
}

// ─── Main update ──────────────────────────────────────────────────────────────
void HeliumCoolingSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;
    updateReactorCircuit(dt, state.fusion_power_MW, state.first_wall_temp_K);
    updateCryostat(dt);
    updateMagnetCircuit(dt, state.magnet_temp_K, state.stored_energy_GJ);

    // Write summaries
    state.reactor_he_outlet_K  = s_.reactor_circuit.outlet_temp_K;
    state.cryostat_temp_K      = s_.magnet_circuit.cold_mass_temp_K;
    state.magnet_he_temp_K     = s_.magnet_circuit.cold_mass_temp_K;
    state.cryo_ok              = (s_.magnet_circuit.cold_mass_temp_K < 6.f) &&
                                  s_.cryostat.vacuum_ok;
    state.cryo_compressor_on   = s_.magnet_circuit.refrigerator_on;
}

void HeliumCoolingSystem::startReactorCooling()
{
    for (auto& p : s_.reactor_circuit.pump) { p.running = true; p.speed_frac = 1.f; }
}
void HeliumCoolingSystem::stopReactorCooling()
{
    for (auto& p : s_.reactor_circuit.pump) p.running = false;
}
void HeliumCoolingSystem::startCryoplant()
{
    s_.magnet_circuit.refrigerator_on = true;
    for (auto& p : s_.magnet_circuit.cryo_pump) { p.running = true; p.speed_frac = 1.f; }
}
void HeliumCoolingSystem::stopCryoplant()
{
    s_.magnet_circuit.refrigerator_on = false;
    for (auto& p : s_.magnet_circuit.cryo_pump) p.running = false;
}
void HeliumCoolingSystem::startCryostatPumping()
{
    s_.cryostat.roughing_pump_on = true;
    // After 5 min (handled outside) switch to turbo
    s_.cryostat.turbo_pump_on = (s_.cryostat.vacuum_pressure_Pa < 5.f);
}