//
// src/VacuumVessel/VacuumVessel.cpp
//  Vacuum vessel pumping and conditioning system.
//

#include "VacuumVessel.h"
#include <cmath>
#include <algorithm>

VacuumVesselSystem::VacuumVesselSystem(const VacuumConfig& cfg)
    : cfg_(cfg)
{}

void VacuumVesselSystem::reset()
{
    pressure_Pa_         = 101325.f;  // back to atmospheric
    wall_temp_K_         = 300.f;
    wall_outgas_factor_  = 1.0f;      // unbaked
    roughing_on_         = false;
    turbo_on_            = false;
    bakeout_on_          = false;
    bakeout_complete_    = false;
    bakeout_elapsed_s_   = 0.f;
    boronization_thickness_nm_ = 0.f;
    boronization_in_progress_  = false;
    boronization_elapsed_s_    = 0.f;
    breach_detected_     = false;
    breach_reason_.clear();
}

void VacuumVesselSystem::startTurbo()
{
    // Interlock: can't start turbo above 10 Pa (bearing damage)
    if (pressure_Pa_ > cfg_.turbo_max_Pa) return;
    turbo_on_ = true;
}

void VacuumVesselSystem::startBoronization()
{
    // Can only boronize if vacuum is good and plasma is cold (no plasma to
    // interfere with the diborane gas)
    if (pressure_Pa_ > cfg_.initiation_max_Pa) return;
    boronization_in_progress_ = true;
    boronization_elapsed_s_   = 0.f;
}

void VacuumVesselSystem::triggerVesselBreach()
{
    breach_detected_ = true;
    breach_reason_   = "Vacuum vessel breach detected (loss-of-vacuum accident)";
    // Open a "hole" — pressure rises rapidly toward atmospheric
    pressure_Pa_ = std::max(pressure_Pa_, 100.f);
}

void VacuumVesselSystem::clearBreach()
{
    breach_detected_ = false;
    breach_reason_.clear();
}

float VacuumVesselSystem::bakeoutProgress() const
{
    if (cfg_.bakeout_duration_s <= 0.f) return 1.f;
    return std::clamp(bakeout_elapsed_s_ / cfg_.bakeout_duration_s, 0.f, 1.f);
}

float VacuumVesselSystem::effectiveOutgasRate() const
{
    // Outgassing rate is high when wall is unbaked (factor=1), low when baked
    // (factor→0.001).  Decay is exponential with bakeout time.
    float baked_factor = wall_outgas_factor_;
    return cfg_.outgas_rate_unbaked * baked_factor
         + cfg_.outgas_rate_Pa_m3_s_m2 * (1.f - baked_factor);
}

void VacuumVesselSystem::updateBakeout(float dt)
{
    if (!bakeout_on_) {
        // Cool back to room temp when bakeout is off
        wall_temp_K_ += (300.f - wall_temp_K_) * 0.01f * dt;
        return;
    }

    // Heat wall toward bakeout temperature
    wall_temp_K_ += (cfg_.bakeout_temp_K - wall_temp_K_) * 0.005f * dt;

    // Bakeout progress: outgassing factor decays exponentially with time at temp
    if (wall_temp_K_ > cfg_.bakeout_temp_K - 10.f) {
        bakeout_elapsed_s_ += dt;
        // Exponential decay of outgassing load (1 → 0.001 over ~12 h)
        constexpr float tau = 12.f * 3600.f;  // 12 h
        wall_outgas_factor_ = 0.001f + 0.999f * expf(-bakeout_elapsed_s_ / tau);

        if (bakeout_elapsed_s_ >= cfg_.bakeout_duration_s) {
            bakeout_complete_ = true;
            // Auto-stop bakeout when complete
            bakeout_on_ = false;
        }
    }
}

void VacuumVesselSystem::updateBoronization(float dt)
{
    if (!boronization_in_progress_) return;

    // Boronization deposits ~50 nm/h of boron film using diborane gas
    // discharge.  During boronization, vacuum is poor (~1 Pa) due to the
    // diborane — plasma cannot run.
    constexpr float deposition_rate_nm_per_s = 50.f / 3600.f;
    boronization_thickness_nm_ += deposition_rate_nm_per_s * dt;
    boronization_elapsed_s_ += dt;

    // Hold pressure at ~1 Pa during boronization
    pressure_Pa_ += (1.f - pressure_Pa_) * 0.1f * dt;

    if (boronization_elapsed_s_ >= boronization_duration_s) {
        boronization_in_progress_ = false;
    }
}

void VacuumVesselSystem::updatePumping(float dt, float plasma_load_Pa_m3_s)
{
    if (breach_detected_) {
        // Pressure rises rapidly toward atmospheric through the breach
        constexpr float breach_rate = 1.0f;  // 1/s — fast leak
        pressure_Pa_ += (101325.f - pressure_Pa_) * breach_rate * dt;
        return;
    }

    if (boronization_in_progress_) return;  // pressure controlled by boronization

    // Gas load: outgassing from walls + plasma gas load (if plasma running)
    float outgas_load = effectiveOutgasRate() * cfg_.vessel_surface_m2;
    float total_load = outgas_load + plasma_load_Pa_m3_s;

    // Roughing pump: S = speed * (P > roughing_min)
    float roughing_throughput = 0.f;
    if (roughing_on_ && pressure_Pa_ > cfg_.roughing_min_Pa
                     && pressure_Pa_ < cfg_.roughing_max_Pa) {
        // Throughput = S * P (constant volumetric speed)
        roughing_throughput = cfg_.roughing_speed_m3s * pressure_Pa_;
    }

    // Turbo pump: S = speed * (P < turbo_max)
    float turbo_throughput = 0.f;
    if (turbo_on_ && pressure_Pa_ < cfg_.turbo_max_Pa) {
        // Turbo has constant volumetric speed down to its ultimate pressure
        float effective_P = std::max(pressure_Pa_, cfg_.turbo_min_Pa);
        turbo_throughput = cfg_.turbo_speed_m3s * effective_P;
    }

    // dP/dt = (load - throughput) / V
    float net_throughput = total_load - roughing_throughput - turbo_throughput;
    float dP = net_throughput / cfg_.vessel_volume_m3 * dt;
    pressure_Pa_ += dP;
    pressure_Pa_ = std::max(pressure_Pa_, 1e-8f);  // clamp above zero

    // Auto-stop turbo if pressure rises above its limit (bearing protection)
    if (turbo_on_ && pressure_Pa_ > cfg_.turbo_max_Pa) {
        turbo_on_ = false;
    }
}

void VacuumVesselSystem::checkBreach(float dt)
{
    if (breach_detected_) return;

    // Detect breach: pressure rising rapidly without pumps on, or rising
    // above 1000 Pa unexpectedly during operation
    static float prev_pressure = 101325.f;
    float dP_dt = (pressure_Pa_ - prev_pressure) / std::max(dt, 1e-6f);
    prev_pressure = pressure_Pa_;

    if (pressure_Pa_ > cfg_.breach_max_Pa
        && !roughing_on_ && !turbo_on_
        && pressure_Pa_ > 100.f) {
        // Already at high pressure — probably just haven't pumped down yet
        return;
    }

    // If pressure is rising > 1000 Pa/s while pumps are on, suspect breach
    if (dP_dt > 1000.f && (roughing_on_ || turbo_on_)
        && pressure_Pa_ < cfg_.breach_max_Pa) {
        triggerVesselBreach();
    }
}

void VacuumVesselSystem::update(ReactorState& state, const SimTime& t)
{
    float dt = t.dt_s;

    // Mirror operator commands from ReactorState → internal state
    // (UI sets state.vessel_roughing_on etc.; we honor those)
    if (state.vessel_roughing_on) roughing_on_ = true;
    else                          roughing_on_ = false;
    if (state.vessel_turbo_on)    turbo_on_    = true;
    else                          turbo_on_    = false;
    if (state.vessel_bakeout_on)  bakeout_on_  = true;
    else                          bakeout_on_  = false;

    // Update bakeout (affects wall outgassing)
    updateBakeout(dt);

    // Update boronization (if requested)
    if (state.boronization_in_progress) startBoronization();
    updateBoronization(dt);

    // Plasma gas load: when plasma is running, it adds a small load to the
    // vessel (gas puffing + recycling).  Roughly proportional to fuel rate.
    float plasma_load = 0.f;
    if (state.plasma_status != PlasmaStatus::Cold) {
        plasma_load = state.sp_fuel_rate * 0.1f;  // 0.1 Pa·m³/s at full fuel rate
    }

    // Update pumping
    updatePumping(dt, plasma_load);

    // Check for vessel breach
    checkBreach(dt);

    // ── Write outputs to ReactorState ───────────────────────────────────────
    state.vessel_pressure_Pa   = pressure_Pa_;
    state.vessel_roughing_on   = roughing_on_;
    state.vessel_turbo_on      = turbo_on_;
    state.vessel_bakeout_on    = bakeout_on_;
    state.vessel_wall_temp_K   = wall_temp_K_;
    state.vessel_vacuum_ok     = vacuumOK();
    state.boronization_thickness_nm = boronization_thickness_nm_;
    state.boronization_in_progress = boronization_in_progress_;
    state.last_boronization_time_s = (boronization_elapsed_s_ > 0.f)
                                     ? (float)t.total_s : state.last_boronization_time_s;
    state.alarm_vacuum_loss    = !vacuumOK() && (pressure_Pa_ > 1e-2f);
    state.alarm_vessel_breach  = breach_detected_;
}
