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

// ─── tritiumPermeationTwoRegime ───────────────────────────────────────────────
//
//  Two-regime permeation model, distinguishing:
//    (a) Diffusion-limited (Sieverts, low pressure, atomic H):
//          J = P_s(T) · (√p_up - √p_down) / L
//        Valid when surface recombination is fast compared to bulk diffusion.
//    (b) Recombination-limited (high pressure, oxidized surfaces):
//          J = k_r(T) · C_surf²   where  C_surf = K_S(T) · √p
//        Valid for surfaces with oxide layers (which act as T barriers).
//
//  The function selects the regime based on a transition parameter
//  (ratio of characteristic times: τ_diff = L²/D, τ_rec = 1/(k_r · C_surf)).
//  Below 500 °C with oxidized surfaces, regime (b) dominates and reduces
//  permeation by 10²-10⁴× (which is why Al₂O₃ / Er₂O₃ barriers work).
//
//  Also includes a barrier transmission factor (0-1) for coated walls.
//
//  References:
//    [1] Causey, J. Nucl. Mater. 300, 91 (2002) — T permeation in fusion.
//    [2] Forcey, J. Nucl. Mater. 160, 117 (1988) — SS316L permeability.
//    [3] Frauenfelder, J. Vac. Sci. Technol. 6, 388 (1969) — W permeability.
//
struct TwoRegimePermeationResult {
    float flux_mol_m2s;             // total flux through the wall [mol/(m²·s)]
    float J_diffusion_limited;      // what regime (a) would give
    float J_recombination_limited;  // what regime (b) would give
    bool  regime_is_diffusion;      // true=(a), false=(b)
    float barrier_factor;           // 0-1 (1 = no barrier)
    float total_loss_mol_s;
    float total_loss_Bq;
};

TwoRegimePermeationResult tritiumPermeationTwoRegime(
    float p_upstream_Pa,
    float p_downstream_Pa,
    float T_wall_K,
    float thickness_m,
    float area_m2,
    WallMaterial material,
    float barrier_factor = 1.0f,        // 1.0 = no barrier; 1e-3 = Al₂O₃ coated
    bool  surface_oxidized = false);     // oxidized → recombination-limited

// ─── TBR self-sufficiency target ──────────────────────────────────────────────
//
//  Post-2020 consensus (Abdou 2021, Fusion Eng. Des. 167, 112374):
//    TBR_design ≥ 1.10   (replaces the old 1.05 floor)
//  This accounts for:
//    - 12.3-year tritium decay (~5.5%/yr loss)
//    - Permeation losses through piping (~2-3%)
//    - Hold-up inventory in the breeder and ISS (~3-5%)
//    - Measurement uncertainty (~3%)
//
//  Net TBR = breeding_rate / burn_rate  must exceed this for self-sufficiency.
//  For P_fus = 3 GWth, burn rate ≈ 1.7×10²⁰ T/s.
//
struct TBRResult {
    float tbr_calculated;        // TBR from blanket physics
    float tbr_required;          // ≥ 1.10
    bool  self_sufficient;       // TBR_calc ≥ TBR_req
    float doubling_time_yr;      // time to double the start-up inventory
};

TBRResult checkTBRSelfSufficiency(float tbr_calculated,
                                   float burn_rate_atoms_s,
                                   float initial_inventory_g);

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