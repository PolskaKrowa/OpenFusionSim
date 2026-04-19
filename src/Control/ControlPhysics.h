#pragma once

//
// src/Control/ControlPhysics.h
// Physics functions for plasma control systems.
//
//  Covers:
//    - Resistive and bootstrap current drive
//    - Disruption force estimates during current quench
//    - Generic PID feedback controller
//    - Runaway electron threshold (Dreicer field)
//

#include <cstdint>

namespace ControlPhysics {

// ─── Geometry descriptor for current calculations ─────────────────────────────
struct ToroidalGeometry {
    float R_major_m;    // major radius [m]
    float a_minor_m;    // minor radius [m]
    float kappa;        // plasma elongation (b/a)
    float delta;        // triangularity
};

// ─── plasmaCurrent ─────────────────────────────────────────────────────────────
//
//  Ohm's law in the plasma: J = σ * E  →  I_p = (V_loop / R_plasma)
//  where R_plasma = η * 2πR / A_plasma  (η = plasma resistivity [Ω·m])
//
//  Uses Spitzer resistivity: η = 5.2e-5 * Z_eff * ln(Λ) / T_keV^1.5  [Ω·m]
//  Here we accept η directly so the caller can supply neoclassical corrections.
//
//  Returns plasma current [A].
//
float plasmaCurrent(float eta_ohm_m,            // plasma resistivity [Ω·m]
                    float V_loop_V,             // loop voltage [V]
                    const ToroidalGeometry& geo);

// ─── bootstrapCurrent ─────────────────────────────────────────────────────────
//
//  Neoclassical bootstrap current fraction:
//    f_bs ≈ C_bs * sqrt(ε) * β_p
//  where ε = a/R (inverse aspect ratio), β_p = poloidal beta,
//  C_bs ≈ 0.3 (empirical prefactor for circular cross-section).
//
//  Full Sauter formula used here (simplified to 3-term fit).
//  Returns bootstrap current [A].
//
float bootstrapCurrent(float p_Pa,              // volume-average plasma pressure [Pa]
                       float B_T,               // toroidal field [T]
                       float epsilon,            // inverse aspect ratio a/R
                       float q_safety,          // safety factor
                       float I_p_A);            // total plasma current [A] (for β_p)

// ─── disruption_dI_dt ─────────────────────────────────────────────────────────
//
//  During a current quench, the plasma current decays in τ_CQ ~ 1–10 ms.
//  The induced halo current and Lorentz force on the structure are:
//    F = I_p * dI/dt * L_plasma   (electromagnetic force [N])
//  and the energy transferred to the wall:
//    E_halo ≈ 0.5 * L_plasma * I_p^2
//
//  Returns dI/dt [A/s] (negative — this is a decay rate).
//
float disruption_dI_dt(float L_plasma_H,        // plasma inductance [H]  (~10 μH)
                       float I_p_A,             // pre-disruption current [A]
                       float tau_CQ_s);         // current quench time [s]

//  Returns the peak halo current force on the vacuum vessel [N]
float disruptionHaloForce(float L_plasma_H,
                          float I_p_A,
                          float tau_CQ_s,
                          float halo_fraction);  // typical 0.1–0.4

// ─── feedbackControlPID ────────────────────────────────────────────────────────
//
//  Discrete-time PID controller with integral anti-windup.
//
//  u(t) = Kp*e + Ki*∫e dt + Kd*(de/dt)
//
//  The integrator is clamped to [-windup_limit, +windup_limit] to prevent
//  integral windup during saturation (e.g. during plasma startup before
//  the controlled variable has reached its setpoint).
//
struct PIDState {
    float integral   = 0.0f;
    float prev_error = 0.0f;
    float output     = 0.0f;
};

float feedbackControlPID(float setpoint,
                         float measured,
                         float Kp,
                         float Ki,
                         float Kd,
                         float dt,
                         PIDState& state,
                         float output_min = -1e9f,
                         float output_max =  1e9f,
                         float windup_limit = 1e6f);

// ─── runawayElectronThreshold ──────────────────────────────────────────────────
//
//  Dreicer electric field: E_D = n_e * e^3 * ln(Λ) / (4π ε₀² * m_e * v_th²)
//  Runaway electrons are accelerated when E > E_D.
//
//  In SI:  E_D [V/m] = 1.44e-9 * n_e * ln(Λ) / T_eV   (with T in eV)
//
//  Also computes the Connor-Hastie threshold accounting for radiation losses
//  and partial screening by impurities.
//
//  Returns:
//    > 0  → runaway growth rate [s^-1]  (if E > E_D)
//    = 0  → below threshold (safe)
//
struct RunawayResult {
    float E_Dreicer_Vm;         // Dreicer field [V/m]
    float E_Connor_Hastie_Vm;   // Connor-Hastie critical field [V/m]
    float growth_rate_s;        // runaway growth rate [1/s]; 0 if sub-threshold
    bool  runaway_risk;
};

RunawayResult runawayElectronThreshold(float E_Vm,          // applied electric field [V/m]
                                       float n_e_m3,        // electron density [m^-3]
                                       float Z_eff,         // effective ion charge
                                       float T_eV,          // electron temperature [eV]
                                       float coulomb_log);  // ln(Λ)

} // namespace ControlPhysics