#pragma once

//
// src/Magnets/MagnetPhysics.h
// Low-level magnet physics calculations.
//
//  Covers:
//    - Toroidal field on axis (Biot-Savart for solenoid)
//    - Poloidal field coil Lorentz forces
//    - Quench propagation (normal zone growth in Nb₃Sn / NbTi)
//    - Cryostat static heat leak
//    - Inductive stored energy
//

#include <cstdint>

namespace MagnetPhysics {

// ─── toroidalFieldOnAxis ───────────────────────────────────────────────────────
//
//  For a multi-turn toroidal solenoid (TF coil set):
//    B = μ₀ * N * I / (2π * R)
//
//  This is the on-axis (at R_major) value.  The 1/R dependence means the
//  field is stronger on the inner (high-field) side of the torus.
//
//  Returns B [T].
//
float toroidalFieldOnAxis(int   N_turns,     // total number of turns
                          float I_A,         // current per turn [A]
                          float R_major_m);  // major radius at axis [m]

// ─── poloidalFieldCoilForces ──────────────────────────────────────────────────
//
//  Lorentz force between two parallel circular coils (Neumann formula):
//    F = μ₀ * I₁ * I₂ * M / L
//  where M is the mutual inductance (computed via elliptic integrals —
//  approximated here with the Neumann formula for coaxial loops).
//
//  Positive F = repulsive (same current direction → attract for coaxial).
//  Returns net axial force [N].
//
struct CoilForceResult {
    float F_axial_N;       // net axial (vertical) force [N]
    float F_radial_N;      // net radial (centering / decentering) force [N]
    float mutual_H;        // mutual inductance between the two coils [H]
};

CoilForceResult poloidalFieldCoilForces(float I1_A,          // current in coil 1 [A]
                                         float I2_A,          // current in coil 2 [A]
                                         float R1_m,          // radius of coil 1 [m]
                                         float R2_m,          // radius of coil 2 [m]
                                         float d_m);          // axial separation [m]

// ─── quenchPropagation ────────────────────────────────────────────────────────
//
//  Normal-zone propagation velocity in a superconducting strand:
//    v_nz = (J * ρ_n / (C_v * (T_cs - T_op)))^0.5
//  where:
//    ρ_n  = normal-state resistivity [Ω·m]  (Cu stabiliser dominates)
//    C_v  = volumetric heat capacity [J/(m³·K)]
//    T_cs = current-sharing temperature (where J > J_c(T,B))
//    T_op = operating temperature
//
//  Also estimates time to quench detection and energy deposited in the hot spot.
//
enum class SCMaterial { Nb3Sn, NbTi, REBCO };

struct QuenchResult {
    float v_normal_zone_ms;   // normal-zone front velocity [m/s]
    float T_hotspot_K;        // peak hot-spot temperature [K]
    float t_detect_s;         // estimated detection time [s]
    float energy_deposited_J; // energy in the normal zone before dump [J]
    bool  quench_confirmed;   // true if T > T_critical
};

QuenchResult quenchPropagation(float T_K,           // local temperature [K]
                               float B_T,           // local field [T]
                               float J_Am2,         // operating current density [A/m²]
                               SCMaterial material,
                               float strand_length_m,
                               float detect_threshold_V); // quench voltage threshold [V]

// ─── cryostatHeatLeak ─────────────────────────────────────────────────────────
//
//  Static heat in-leak through insulation to the cold mass:
//    Q_leak = k_eff * A * (T_warm - T_cold) / thickness
//
//  For multi-layer insulation (MLI) + support strut contributions.
//  k_eff is the effective conductivity of the insulation system.
//
//  Returns total heat load [W].
//
struct CryostatHeatLeak {
    float Q_radiation_W;   // radiation through MLI [W]
    float Q_conduction_W;  // conduction through support struts/pipes [W]
    float Q_total_W;       // sum [W]
};

CryostatHeatLeak cryostatHeatLeak(float A_m2,              // cold-mass surface area [m²]
                                   float k_insulation_Wm_K, // effective k of insulation [W/m/K]
                                   float thickness_m,        // insulation thickness [m]
                                   float T_warm_K,           // warm boundary temperature [K]
                                   float T_cold_K,           // cold mass temperature [K]
                                   int   n_MLI_layers);      // number of MLI layers (0 = none)

// ─── inductiveStoredEnergy ────────────────────────────────────────────────────
//
//  Energy stored in an inductor:  E = ½ L I²
//  For the TF coil set: L is the total inductance of all coils in series.
//  For the plasma: L_plasma ≈ μ₀ * R * (ln(8R/a) - 2 + li/2)  (Shafranov formula)
//
//  Returns stored energy [J].
//
float inductiveStoredEnergy(float L_H,    // inductance [H]
                             float I_A);  // current [A]

// Plasma self-inductance via Shafranov formula
float plasmaInductance(float R_major_m,    // major radius [m]
                       float a_minor_m,    // minor radius [m]
                       float li,           // internal inductance per unit length (~ 0.7)
                       float kappa);       // elongation

} // namespace MagnetPhysics