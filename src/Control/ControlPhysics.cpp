//
// src/Control/ControlPhysics.cpp
//

#include "ControlPhysics.h"
#include <cmath>
#include <algorithm>

namespace ControlPhysics {

// ─── Physical constants ────────────────────────────────────────────────────────
static constexpr float PI      = 3.14159265358979f;
static constexpr float MU0     = 1.25663706e-6f;   // [H/m]
static constexpr float EPS0    = 8.85418782e-12f;  // [F/m]
static constexpr float E_CHARGE= 1.60217663e-19f;  // [C]
static constexpr float M_ELEC  = 9.10938370e-31f;  // [kg]
static constexpr float K_B     = 1.38064852e-23f;  // [J/K]

// ─── plasmaCurrent ─────────────────────────────────────────────────────────────
float plasmaCurrent(float eta_ohm_m,
                    float V_loop_V,
                    const ToroidalGeometry& geo)
{
    // Plasma cross-sectional area: A ≈ π * a² * κ  (elongated ellipse)
    float A_plasma = PI * geo.a_minor_m * geo.a_minor_m * geo.kappa;

    // Plasma column length: L ≈ 2π * R_major (toroidal circumference)
    float L_column = 2.0f * PI * geo.R_major_m;

    // Resistance:  R = η * L / A   [Ω]
    float R_plasma = eta_ohm_m * L_column / A_plasma;

    if (R_plasma < 1e-20f) return 0.0f;

    // Ohm's law:  I_p = V_loop / R_plasma
    return V_loop_V / R_plasma;
}

// ─── bootstrapCurrent ─────────────────────────────────────────────────────────
float bootstrapCurrent(float p_Pa,
                       float B_T,
                       float epsilon,
                       float q_safety,
                       float I_p_A)
{
    if (B_T < 1e-3f || I_p_A < 1.0f) return 0.0f;

    // Poloidal beta: β_p = 2μ₀ <p> / B_p²
    // B_poloidal ≈ μ₀ * I_p / (2π * a)   —  but we use β_p via energy relation:
    // β_p = (2μ₀ <p>) / (B_p²)  where B_p ~ μ₀ * I_p / (2π * a)
    // Simplified: β_p ~ (2μ₀ <p> * (2π * R)² ) / (μ₀² * I_p² * ...)
    // For game purposes use the standard form β_p = (8π / (μ₀ I_p²)) * ∫p dV
    // ∫p dV ≈ p * V_plasma  where V ~ 2π²Ra²κ
    float a = epsilon;   // caller passes ε = a/R, we need to separate — accept as ratio
    // Treat epsilon as the geometric inverse aspect ratio directly.

    // Sauter simplified bootstrap fraction fit (valid for 0.1 < ε < 0.5):
    //   f_bs ≈ C1 * sqrt(ε) / q * (1 + 1.4*ε / (1-ε))^-1 * β_p
    // We compute β_p from pressure and current directly.
    float B_p_approx = MU0 * I_p_A / (2.0f * PI * 2.0f);  // B_pol ~ μ₀Ip/(2πa); a~2m
    if (B_p_approx < 1e-6f) return 0.0f;

    float beta_p = (2.0f * MU0 * p_Pa) / (B_p_approx * B_p_approx);

    float C_bs   = 0.3f;  // empirical, circular cross-section
    float f_bs   = C_bs * sqrtf(epsilon) * beta_p
                 * (1.0f + 1.4f * epsilon / (1.0f - epsilon + 1e-6f));

    // Clamp: bootstrap can't exceed 100 % of plasma current
    f_bs = std::clamp(f_bs, 0.0f, 1.0f);

    return f_bs * I_p_A;  // [A]
}

// ─── disruption_dI_dt ─────────────────────────────────────────────────────────
float disruption_dI_dt(float L_plasma_H,
                       float I_p_A,
                       float tau_CQ_s)
{
    if (tau_CQ_s < 1e-6f) tau_CQ_s = 1e-6f;

    // Exponential decay model:  I(t) = I_p * exp(-t/τ_CQ)
    // dI/dt at t=0:  = -I_p / τ_CQ
    // The inductance appears in the energy: E = 0.5 * L * I²
    // Force on vessel: F = I * (dΦ/dt) = I * L * (dI/dt)
    (void)L_plasma_H;  // used in halo force calculation

    return -I_p_A / tau_CQ_s;  // [A/s]
}

float disruptionHaloForce(float L_plasma_H,
                          float I_p_A,
                          float tau_CQ_s,
                          float halo_fraction)
{
    float dI_dt = disruption_dI_dt(L_plasma_H, I_p_A, tau_CQ_s);
    float I_halo = halo_fraction * I_p_A;

    // F = I_halo * L * |dI/dt|  — vertical force on the vessel [N]
    return I_halo * L_plasma_H * fabsf(dI_dt);
}

// ─── feedbackControlPID ────────────────────────────────────────────────────────
float feedbackControlPID(float setpoint,
                         float measured,
                         float Kp,
                         float Ki,
                         float Kd,
                         float dt,
                         PIDState& state,
                         float output_min,
                         float output_max,
                         float windup_limit)
{
    float error      = setpoint - measured;
    float derivative = (dt > 1e-10f) ? (error - state.prev_error) / dt : 0.0f;

    // Integrate with anti-windup clamping
    state.integral += error * dt;
    state.integral  = std::clamp(state.integral, -windup_limit, windup_limit);

    float u = Kp * error
            + Ki * state.integral
            + Kd * derivative;

    // Clamp output and back-calculate integral to prevent windup at saturation
    float u_clamped = std::clamp(u, output_min, output_max);
    if (u != u_clamped && Ki > 1e-12f) {
        // Back-calculate: remove the excess that caused saturation
        state.integral -= (u - u_clamped) / Ki;
        state.integral  = std::clamp(state.integral, -windup_limit, windup_limit);
    }

    state.prev_error = error;
    state.output     = u_clamped;
    return u_clamped;
}

// ─── runawayElectronThreshold ──────────────────────────────────────────────────
RunawayResult runawayElectronThreshold(float E_Vm,
                                       float n_e_m3,
                                       float Z_eff,
                                       float T_eV,
                                       float coulomb_log)
{
    RunawayResult res{};

    if (n_e_m3 < 1.0f || T_eV < 0.1f) {
        res.runaway_risk = false;
        return res;
    }

    // ── Dreicer field ──────────────────────────────────────────────────────────
    // E_D = n_e * e³ * ln(Λ) / (4π * ε₀² * m_eLEC * v_th²)
    // v_th = sqrt(k_B * T / m_eLEC)  [but T in eV → T_J = T_eV * e]
    float T_J    = T_eV * E_CHARGE;                       // [J]
    float v_th2  = T_J / M_ELEC;                             // thermal speed squared [m²/s²]
    float E_D    = n_e_m3 * E_CHARGE * E_CHARGE * E_CHARGE * coulomb_log
                 / (4.0f * PI * EPS0 * EPS0 * M_ELEC * v_th2);
    res.E_Dreicer_Vm = E_D;

    // ── Connor-Hastie critical field ──────────────────────────────────────────
    // E_CH = (4π * ε₀² * m_eLEC * c² * n_e * e * ln(Λ)) / ...
    // Simplified form:  E_CH ≈ n_e * e³ * ln(Λ) / (4π * ε₀² * m_eLEC * c²)
    constexpr float C = 2.99792458e8f;
    float E_CH = n_e_m3 * E_CHARGE * E_CHARGE * E_CHARGE * coulomb_log
               / (4.0f * PI * EPS0 * EPS0 * M_ELEC * C * C);
    res.E_Connor_Hastie_Vm = E_CH;

    // ── Runaway growth rate ────────────────────────────────────────────────────
    // Avalanche growth: γ ≈ (E/E_CH - 1) / (τ_coll * ln(Λ))
    // where τ_coll = 4π ε₀² m_eLEC^2 v_th^3 / (n_e e^4 ln(Λ))
    float tau_coll = 4.0f * PI * EPS0 * EPS0 * M_ELEC * M_ELEC * powf(v_th2, 1.5f)
                   / (n_e_m3 * powf(E_CHARGE, 4.0f) * coulomb_log);

    if (E_Vm > E_CH && tau_coll > 1e-30f) {
        float E_ratio = E_Vm / E_CH;
        // Avalanche: gamma ~ (E/E_CH - 1)^0.5 / (τ_coll * ln(Λ))  (Connor-Hastie)
        // Accounts for Z_eff through the effective collision rate
        float Z_factor = (Z_eff + 1.0f) / 2.0f;
        res.growth_rate_s = sqrtf(E_ratio - 1.0f) / (Z_factor * tau_coll * coulomb_log);
        res.runaway_risk  = true;
    } else {
        // Below Connor-Hastie: check Dreicer primary generation
        if (E_Vm > E_CH * 0.01f) {
            // Primary Dreicer generation (exponentially suppressed below E_D)
            float x = E_D / std::max(E_Vm, 1e-10f);  // E_D/E
            float exponent = -sqrtf(x) * (1.0f + (Z_eff + 1.0f) / (2.0f));
            res.growth_rate_s = (n_e_m3 / tau_coll) * expf(std::max(exponent, -80.0f));
            res.runaway_risk  = (res.growth_rate_s > 1.0f); // >1 runaway/s
        } else {
            res.growth_rate_s = 0.0f;
            res.runaway_risk  = false;
        }
    }

    return res;
}

} // namespace ControlPhysics