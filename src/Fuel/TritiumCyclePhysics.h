#pragma once

//
// src/Fuel/TritiumCyclePhysics.h
// Tritium fuel cycle: breeding, inventory, decay, permeation, isotope separation.
//
//  All functions from TritiumCycle.md spec:
//    tritiumBurnFraction     — fraction of injected T actually burned
//    tritiumInventory        — dN_T/dt bookkeeping
//    tritiumDecay            — radioactive decay (t½ = 12.32 years)
//    tritiumPermeation       — T₂ migration through steel via Sieverts' law
//    isotopeSeparationWork   — separative work units for D/T/He-3/He-4 streams
//

#include <cstdint>

namespace TritiumCyclePhysics {

// ─── tritiumBurnFraction ──────────────────────────────────────────────────────
//
//  Fraction of the tritium inventory actually fused per confinement time:
//    f_burn = n_T * <σv> * τ_conf / (n_T + n_D)
//  For D-T at peak reactivity (~65 keV, but reactors run at 20 keV):
//    f_burn ≈ 0.02–0.05  (2–5 %) — most T passes through unburned.
//
//  Returns burn fraction [dimensionless, 0–1].
//
float tritiumBurnFraction(float n_T_m3,             // tritium number density [m^-3]
                           float sigma_v_m3s,        // reaction rate <σv> [m³/s]
                           float confinement_time_s); // energy confinement time [s]

// ─── tritiumInventory ─────────────────────────────────────────────────────────
//
//  Tritium inventory balance (atoms or moles):
//    dN_T/dt = Ṡ_breed - Ṡ_burn - Ṡ_decay - Ṡ_permeation - Ṡ_losses
//
//  All rates in [atoms/s] or equivalently [mol/s] if normalised by N_A.
//  Returns dN_T/dt [atoms/s].
//
struct TritiumRates {
    float production_s;     // breeding rate from Li(n,t) [atoms/s]
    float consumption_s;    // fusion burn rate [atoms/s]
    float decay_s;          // radioactive decay rate [atoms/s]
    float permeation_s;     // loss via wall permeation [atoms/s]
    float misc_losses_s;    // handling losses, exhaust, etc. [atoms/s]
};

float tritiumInventory(const TritiumRates& rates); // returns dN_T/dt [atoms/s]

// Overload: pass totals directly for simple bookkeeping
float tritiumInventory(float production_s,    // T bred per second [atoms/s]
                        float consumption_s,  // T fused per second [atoms/s]
                        float decay_s,        // T decayed per second [atoms/s]
                        float losses_s);      // all other losses [atoms/s]

// ─── tritiumDecay ─────────────────────────────────────────────────────────────
//
//  Tritium undergoes beta decay:
//    ³H → ³He + e⁻ + ν̄ₑ   (t½ = 12.32 years)
//
//  N(t) = N₀ * exp(-λ * t)    where  λ = ln(2) / t½
//
//  Returns remaining tritium atoms N(t) after time t.
//  Also provides the He-3 produced (=N₀ - N(t)).
//
struct DecayResult {
    float N_T_remaining;    // tritium atoms remaining
    float N_He3_produced;   // He-3 atoms produced (potentially useful as fuel for D-He3)
    float activity_Bq;      // activity in Becquerels [decays/s]
    float decay_heat_W;     // decay heat from beta particles [W]  (~324 mW/g T)
};

DecayResult tritiumDecay(float N_T_initial,   // initial tritium atom count
                          double t_s);         // elapsed time [s]

// ─── tritiumPermeation ────────────────────────────────────────────────────────
//
//  Hydrogen isotope permeation through metallic walls (Sieverts' law):
//    J = P_s * (sqrt(p_T2_upstream) - sqrt(p_T2_downstream)) / thickness
//
//  where P_s is the permeability [mol/(m·s·Pa^0.5)], temperature-dependent:
//    P_s(T) = P_s0 * exp(-E_a / (R * T))
//
//  Permeability data (316L stainless steel):
//    P_s0 = 4.0e-7  mol/(m·s·Pa^0.5)
//    E_a  = 49.7 kJ/mol
//
enum class WallMaterial { SS316L, Tungsten, EUROFER97, Inconel625 };

struct PermeationResult {
    float flux_mol_m2s;         // tritium flux [mol/(m²·s)]
    float total_loss_mol_s;     // total loss through entire wall area [mol/s]
    float total_loss_Bq;        // activity of permeated T [Bq/s]
    float permeability_m3s;     // permeability at this T [mol/(m·s·Pa^0.5)]
};

PermeationResult tritiumPermeation(float p_upstream_Pa,    // T₂ partial pressure upstream [Pa]
                                    float p_downstream_Pa, // T₂ partial pressure downstream [Pa]
                                    float T_wall_K,        // wall temperature [K]
                                    float thickness_m,     // wall thickness [m]
                                    float area_m2,         // total permeation area [m²]
                                    WallMaterial material);

// ─── isotopeSeparationWork ────────────────────────────────────────────────────
//
//  Separative Work Units (SWU) for isotope separation (e.g. D/T/He-3/He-4).
//  Used for sizing the Isotope Separation System (ISS) in the tritium plant.
//
//  SWU = P * V(x_p) + W * V(x_w) - F * V(x_f)
//  where V(x) = (2x - 1) * ln(x / (1-x))  is the value function
//        P, W, F = product, waste, feed flow rates [mol/s or g/s]
//        x_p, x_w, x_f = product, waste, feed fractions
//        α = separation factor (stage-wise enrichment ratio)
//
//  Returns SWU [mol·SWU/s] and minimum theoretical stages.
//
struct SeparationResult {
    float SWU_mol_s;        // separative work [mol·SWU/s]
    float W_waste_mol_s;    // waste stream flow [mol/s]
    float P_product_mol_s;  // product stream flow [mol/s]
    float N_stages;         // minimum theoretical stages (from α)
    float energy_kJ_mol;    // theoretical energy [kJ per mol of product]
};

SeparationResult isotopeSeparationWork(float F_feed_mol_s,   // feed stream [mol/s]
                                        float x_feed,         // feed composition (target species fraction)
                                        float x_product,      // product enrichment fraction
                                        float x_waste,        // waste tail fraction
                                        float alpha);         // separation factor (α > 1)

} // namespace TritiumCyclePhysics