//
// src/TurbineSystem/TurbineUnit.cpp
//

#include "TurbineUnit.h"
#include <cmath>
#include <algorithm>

static constexpr float PI          = 3.14159265f;
static constexpr float RPM_NOMINAL = 3000.f;    // 50 Hz, 2-pole
static constexpr float FREQ_NOM    = 50.f;
static constexpr float P_RATED_MW  = 1200.f;
static constexpr float H_INERTIA   = 6.f;       // inertia constant [s]

void TurbineUnitController::initPreheaters()
{
    // 4 stages: extraction pressures roughly halved each stage
    float temps[4] = {380.f, 430.f, 470.f, 500.f}; // outlet temps [K]
    float extr[4]  = {0.5f,  1.5f,  4.0f,  8.0f};  // extraction pressures [MPa]
    for (int i = 0; i < 4; i++) {
        s.preheater[i].enabled              = true;
        s.preheater[i].fw_outlet_temp_K     = temps[i];
        s.preheater[i].extraction_frac      = 0.05f;
        s.preheater[i].drain_level_m        = 0.5f;
        s.preheater[i].heat_transferred_MW  = 0.f;
    }
    s.fw_temp_after_ph_K = 500.f;
}

// ─── State machine ────────────────────────────────────────────────────────────
void TurbineUnitController::runStateMachine(float dt, float grid_freq_Hz)
{
    switch (s.state) {
    case TurbineState::Offline:
        s.msiv_open = false;
        s.governor_demand = 0.f;
        s.generator.exciter_on = false;
        break;

    case TurbineState::RollingUp:
        // Open MSIV to ~10%, let turbine spin up
        s.msiv_open = true;
        s.msiv_setpoint = 0.12f;
        s.generator.exciter_on = true;
        // Transition when rpm close to nominal
        if (s.rpm > RPM_NOMINAL * 0.98f && s.rpm < RPM_NOMINAL * 1.02f) {
            s.state = TurbineState::Synchronizing;
            sync_timer_s_ = 0.f;
            phase_error_integrated_ = 0.f;
        }
        break;

    case TurbineState::Synchronizing: {
        s.msiv_open = true;
        sync_timer_s_ += dt;

        // Fine-tune RPM to match grid
        float freq_err = grid_freq_Hz - s.generator.frequency_Hz;
        // Adjust governor slightly to correct frequency
        s.governor_demand = std::clamp(0.12f + freq_err * 0.002f, 0.05f, 0.25f);
        s.msiv_setpoint   = s.governor_demand;

        // Check sync conditions
        float dfreq  = std::abs(s.generator.frequency_Hz - grid_freq_Hz);
        float dphase = std::abs(s.generator.phase_rad
                                - std::fmodf(2.f * PI * grid_freq_Hz * sync_timer_s_, 2.f * PI));
        dphase = std::min(dphase, 2.f * PI - dphase); // wrap

        if (breaker_close_requested_ && dfreq < 0.15f && dphase < 0.26f) { // 15° deg
            s.generator.breaker_closed = true;
            s.state = TurbineState::Online;
            breaker_close_requested_ = false;
        }
        // Timeout safety: if can't sync in 5 minutes, trip
        if (sync_timer_s_ > 300.f) {
            s.state = TurbineState::Tripping;
        }
        break;
    }

    case TurbineState::Online:
        // Governor tracks frequency
        {
            float freq_err = grid_freq_Hz - s.generator.frequency_Hz;
            // Droop control: 4% droop → 50 MW change per 0.1 Hz error
            float droop_demand = s.governor_demand + freq_err * 0.01f;
            s.governor_demand = std::clamp(droop_demand, 0.f, 1.f);
            s.msiv_setpoint   = s.governor_demand;
            s.msiv_open       = true;
        }
        break;

    case TurbineState::Runback:
        s.governor_demand = std::max(s.governor_demand - 0.05f * dt, 0.1f);
        s.msiv_setpoint   = s.governor_demand;
        if (s.governor_demand <= 0.1f) s.state = TurbineState::Online;
        break;

    case TurbineState::Tripping:
        s.msiv_open           = false;
        s.msiv_setpoint       = 0.f;
        s.governor_demand     = 0.f;
        s.bypass_valve_pos    = std::min(s.bypass_valve_pos + 0.5f * dt, 1.f);
        if (s.generator.breaker_closed && s.generator.power_MW < 10.f) {
            s.generator.breaker_closed = false;
        }
        if (s.rpm < 50.f) {
            s.state = TurbineState::Tripped;
            s.bypass_valve_pos = 0.f;
            s.msiv_trip_latch  = true;
        }
        break;

    case TurbineState::Tripped:
        s.msiv_open           = false;
        s.generator.exciter_on= false;
        s.generator.breaker_closed = false;
        break;
    }

    // Relief valve: open on overpressure
    s.relief_valve_open = (s.sg_pressure_MPa > s.relief_setpoint_MPa);

    // Overspeed trip
    if (s.rpm > RPM_NOMINAL * 1.10f) {
        s.overspeed_trip = true;
        s.alarm_turb_trip = true;
        s.state = TurbineState::Tripping;
    }
}

// ─── Steam Generator ──────────────────────────────────────────────────────────
void TurbineUnitController::updateSteamGenerator(float dt)
{
    // Feedwater pump operation
    float fw_total_flow = 0.f;
    for (int i = 0; i < 2; i++) {
        auto& p = s.fw_pump[i];
        if (p.running && !p.trip) {
            // Full-speed pump provides up to 1000 kg/s at design point
            p.flow_kg_s      = 1000.f * p.speed_frac
                             * (s.state == TurbineState::Online ? 1.f : 0.3f);
            p.discharge_MPa  = s.sg_pressure_MPa + 2.f;
            p.power_MW       = p.flow_kg_s * 2e6f / (3600.f * 0.82f) * 1e-6f; // 2 MJ/t lift
            fw_total_flow   += p.flow_kg_s;
        } else {
            p.flow_kg_s = 0.f; p.power_MW = 0.f; p.discharge_MPa = 0.f;
        }
    }

    // SG pressure dynamics
    // Heat in from molten salt; heat out as steam
    float h_steam_J_kg  = 2.8e6f; // ~2.8 MJ/kg at 18 MPa (approx enthalpy of steam)
    float h_fw_J_kg     = 1000.f * (s.fw_temp_after_ph_K - 273.f); // approx liquid enthalpy
    float Q_to_steam_MW = fw_total_flow * (h_steam_J_kg - h_fw_J_kg) * 1e-6f;

    // Net heat → changes pressure (volume ~300 m³ of steam/water)
    float dQ_net = s.sg_heat_input_MW - Q_to_steam_MW;
    // dP/dt ≈ dQ_net / (cp_water * volume_equiv)
    float dP_dt  = dQ_net * 1e6f / (1e7f); // ~10 GJ/MPa for large SG
    s.sg_pressure_MPa += dP_dt * dt;
    s.sg_pressure_MPa  = std::clamp(s.sg_pressure_MPa, 0.05f, 22.f);

    // Steam temperature from pressure (saturation, approximation)
    s.sg_steam_temp_K = 273.f + 100.f * std::pow(s.sg_pressure_MPa / 0.1013f, 0.25f);
    s.sg_steam_temp_K = std::min(s.sg_steam_temp_K, 650.f);

    // Steam production
    s.sg_steam_flow_kg_s = fw_total_flow; // conservation of mass (simplified)

    // SG level: rises with feedwater, falls with steam production
    float dLevel = (fw_total_flow - s.sg_steam_flow_kg_s) / 800.f * dt; // 800 kg/m
    s.sg_level_m = std::clamp(s.sg_level_m + dLevel, 0.f, 20.f);
    s.alarm_lo_sg_level = (s.sg_level_m < 3.f);
    s.alarm_hi_sg_pressure = (s.sg_pressure_MPa > 19.f);
}

// ─── MSIV and Steam Path ──────────────────────────────────────────────────────
void TurbineUnitController::updateMSIVAndSteamPath(float dt)
{
    // MSIV position ramps toward setpoint
    float target = (s.msiv_open && !s.msiv_trip_latch) ? s.msiv_setpoint : 0.f;
    float rate   = (target > s.msiv_position) ? 0.3f : 1.0f; // opens slow, closes fast
    s.msiv_position += std::clamp(target - s.msiv_position,
                                   -rate * dt, rate * dt);
    s.msiv_position  = std::clamp(s.msiv_position, 0.f, 1.f);

    // Steam flow through turbine (governor valve further throttles)
    float sg_flow   = s.sg_steam_flow_kg_s;
    float gov_flow  = sg_flow * s.msiv_position * s.governor_demand;

    // Bypass takes some if valve open
    float bypass_flow = sg_flow * s.bypass_valve_pos * 0.5f;

    // Relief valve venting
    float relief_flow = s.relief_valve_open
                      ? (s.sg_pressure_MPa - s.relief_setpoint_MPa) * 20.f : 0.f;
    relief_flow = std::max(relief_flow, 0.f);

    s.steam_flow_to_turbine = std::max(gov_flow - bypass_flow, 0.f);
}

// ─── Preheaters ───────────────────────────────────────────────────────────────
void TurbineUnitController::updatePreheaters(float dt)
{
    float fw_temp = 310.f; // cold feedwater from hotwell/condensate pump
    float total_extract = 0.f;

    for (int i = 0; i < 4; i++) {
        auto& ph = s.preheater[i];
        if (ph.enabled && s.state == TurbineState::Online) {
            // Extracts a fraction of steam from the turbine
            ph.heat_transferred_MW = s.steam_flow_to_turbine
                                   * ph.extraction_frac * 2.0e6f * 1e-6f;
            // Heats feedwater toward ph outlet temp
            fw_temp = fw_temp + (ph.fw_outlet_temp_K - fw_temp) * 0.85f;
            total_extract += ph.extraction_frac;
        } else {
            ph.heat_transferred_MW = 0.f;
        }
    }
    s.fw_temp_after_ph_K = fw_temp;
    (void)total_extract;
    (void)dt;
}

// ─── Turbine Shaft ────────────────────────────────────────────────────────────
void TurbineUnitController::updateTurbineShaft(float dt)
{
    if (s.steam_flow_to_turbine < 0.1f && s.state != TurbineState::Online) {
        // Coast-down: exponential decay
        s.rpm = std::max(s.rpm * (1.f - 0.05f * dt), 0.f);
        s.shaft_power_MW = 0.f;
        return;
    }

    // Steam enthalpy drop across turbine
    float h_in_J_kg  = 2.8e6f; // HPT inlet
    float h_out_J_kg = 2.1e6f; // LPT exhaust (varies with condenser pressure)
    float eta_turb   = 0.88f;
    float dh         = (h_in_J_kg - h_out_J_kg) * eta_turb;
    s.shaft_power_MW = s.steam_flow_to_turbine * dh * 1e-6f;

    // Swing equation: df/dt = (P_mech - P_elec) / (2H * P_base / f_nom)
    float P_elec = s.generator.breaker_closed ? s.generator.power_MW : 0.f;
    float df_dt  = (s.shaft_power_MW - P_elec) * FREQ_NOM / (2.f * H_INERTIA * P_RATED_MW);

    float freq = s.rpm / RPM_NOMINAL * FREQ_NOM;
    freq += df_dt * dt;
    freq  = std::max(freq, 0.f);
    s.rpm = freq / FREQ_NOM * RPM_NOMINAL;

    // Friction/windage loss
    s.rpm = std::max(s.rpm - 0.1f * dt, 0.f);
}

// ─── Condenser ────────────────────────────────────────────────────────────────
void TurbineUnitController::updateCondenser(float dt)
{
    auto& c = s.condenser;

    // Air ingress: slowly leaks in, CAR pumps and SJAE remove it
    c.air_fraction += 1e-6f * dt;     // slow ingress
    int car_count = 0;
    for (int i = 0; i < 4; i++) if (c.car_pump[i]) car_count++;
    c.air_fraction -= car_count * 5e-7f * dt;
    if (c.sjae_running) {
        c.air_fraction -= 2e-6f * dt;
        c.sjae_steam_kg_s = 0.8f;     // SJAE uses ~0.8 kg/s of auxiliary steam
    } else {
        c.sjae_steam_kg_s = 0.f;
    }
    c.air_fraction = std::clamp(c.air_fraction, 0.f, 0.10f);

    // Condenser pressure: set by cooling water temperature + air ingress
    // Saturation pressure at cooling water outlet (~CW_in + 8 K approach)
    float T_CW_out = c.cooling_water_temp_K + 8.f;
    // Antoine-like for water near 300 K: P_sat [kPa] ≈ exp(16.5 - 3800/(T-46)) * 0.1
    float T_C = T_CW_out - 273.f;
    float P_sat_kPa = std::exp(16.5f - 3800.f / (T_C + 227.f)) * 0.13f;
    // Air partial pressure
    float P_air_kPa = c.air_fraction * P_sat_kPa * 50.f;
    // Target total pressure
    float P_target  = P_sat_kPa + P_air_kPa;
    // Ramp actual pressure toward target
    float rate      = (P_target > c.pressure_kPa) ? 2.f : 0.5f;
    c.pressure_kPa += std::clamp(P_target - c.pressure_kPa, -rate*dt, rate*dt);
    c.pressure_kPa  = std::clamp(c.pressure_kPa, 2.f, 50.f);

    // Saturation temperature at condenser pressure
    c.temp_K = 273.f + 100.f * std::pow(c.pressure_kPa / 101.3f, 0.25f);
    c.temp_K = std::max(c.temp_K, c.cooling_water_temp_K);

    // Condensate pump flow
    c.condensate_flow_kg_s = c.condensate_pump_speed * 2200.f;

    s.alarm_lo_condenser_vac = (c.pressure_kPa > 15.f);
    (void)dt;
}

// ─── Hotwell ──────────────────────────────────────────────────────────────────
void TurbineUnitController::updateHotwell(float dt)
{
    auto& hw = s.hotwell;
    auto& c  = s.condenser;

    // Condensate in (from condenser) minus feedwater pumps out
    float fw_out = 0.f;
    for (int i = 0; i < 2; i++) fw_out += s.fw_pump[i].flow_kg_s;

    float in_flow  = c.condensate_flow_kg_s;
    float out_flow = fw_out;

    // Makeup and drain
    float makeup = hw.makeup_valve ? hw.makeup_flow_kg_s : 0.f;
    float drain  = hw.drain_valve  ? hw.drain_flow_kg_s  : 0.f;

    float dLevel = (in_flow + makeup - out_flow - drain) / (hw.area_m2 * 1000.f) * dt;
    hw.level_m   = std::clamp(hw.level_m + dLevel, 0.f, 3.5f);

    // Auto level control (deadband ±0.3 m around 1.5 m)
    hw.lo_level_alarm = (hw.level_m < 0.5f);
    hw.hi_level_alarm = (hw.level_m > 3.0f);
    hw.makeup_valve   = (hw.level_m < 1.0f);
    hw.drain_valve    = (hw.level_m > 2.5f);
}

// ─── Generator ────────────────────────────────────────────────────────────────
void TurbineUnitController::updateGenerator(float dt, float grid_freq_Hz)
{
    auto& g = s.generator;
    if (!g.exciter_on) { g.frequency_Hz = 0.f; g.power_MW = 0.f; return; }

    // Generator frequency tracks shaft (instantaneous)
    g.frequency_Hz = s.rpm / RPM_NOMINAL * FREQ_NOM;

    // Phase integrates from frequency
    g.phase_rad += 2.f * PI * g.frequency_Hz * dt;
    if (g.phase_rad > 2.f * PI) g.phase_rad -= 2.f * PI;

    // When connected: power = shaft - generator losses
    if (g.breaker_closed) {
        float eta_gen  = 0.986f;
        g.power_MW     = s.shaft_power_MW * eta_gen;
        g.power_MW     = std::clamp(g.power_MW, 0.f, P_RATED_MW * 1.05f);
        g.reactive_MVAR= g.power_MW * 0.15f; // ~0.15 lagging pf correction
    } else {
        g.power_MW      = 0.f;
        g.reactive_MVAR = 0.f;
    }

    g.overspeed_trip = (s.rpm > RPM_NOMINAL * 1.10f);
    (void)grid_freq_Hz;
}

// ─── Alarms ───────────────────────────────────────────────────────────────────
void TurbineUnitController::checkAlarms()
{
    s.alarm_turb_trip = s.msiv_trip_latch || s.overspeed_trip ||
                        s.generator.overspeed_trip || s.generator.diff_trip;
}

// ─── Main Update ─────────────────────────────────────────────────────────────
void TurbineUnitController::update(float dt, float grid_frequency_Hz)
{
    runStateMachine(dt, grid_frequency_Hz);
    updateSteamGenerator(dt);
    updateMSIVAndSteamPath(dt);
    updatePreheaters(dt);
    updateTurbineShaft(dt);
    updateCondenser(dt);
    updateHotwell(dt);
    updateGenerator(dt, grid_frequency_Hz);
    checkAlarms();
}

// ─── Operator Commands ────────────────────────────────────────────────────────
void TurbineUnitController::cmdStart()
{
    if (s.state == TurbineState::Offline || s.state == TurbineState::Tripped) {
        if (!s.msiv_trip_latch) {
            s.state = TurbineState::RollingUp;
            s.generator.exciter_on = true;
            // Start both FW pumps
            for (auto& p : s.fw_pump) { p.running = true; p.speed_frac = 1.f; }
            // Start condenser systems
            for (auto& car : s.condenser.car_pump) car = true;
            s.condenser.sjae_running = true;
            s.condenser.condensate_pump_speed = 1.f;
        }
    }
}

void TurbineUnitController::cmdStop()
{
    if (s.state == TurbineState::Online || s.state == TurbineState::Synchronizing)
        s.state = TurbineState::Runback;
}

void TurbineUnitController::cmdTrip()
{
    s.state = TurbineState::Tripping;
    s.msiv_trip_latch = true;
}

void TurbineUnitController::cmdReset()
{
    if (s.state == TurbineState::Tripped || s.state == TurbineState::Offline) {
        s.msiv_trip_latch  = false;
        s.overspeed_trip   = false;
        s.generator.overspeed_trip = false;
        s.generator.diff_trip = false;
        s.alarm_turb_trip  = false;
        s.state = TurbineState::Offline;
    }
}

void TurbineUnitController::cmdCloseBreakerRequest()
{
    if (s.state == TurbineState::Synchronizing)
        breaker_close_requested_ = true;
}

void TurbineUnitController::cmdOpenBreaker()
{
    s.generator.breaker_closed = false;
    if (s.state == TurbineState::Online)
        s.state = TurbineState::Runback;
}