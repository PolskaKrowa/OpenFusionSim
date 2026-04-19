//
// src/Fuel/TritiumCyclePhysics.cpp
//

#include "TritiumCyclePhysics.h"
#include <cmath>
#include <algorithm>

namespace TritiumCyclePhysics {

// Physical constants
static constexpr double LN2        = 0.693147180559945;
static constexpr double T_HALF_S   = 12.32 * 365.25 * 24.0 * 3600.0; // 12.32 years [s]
static constexpr double LAMBDA_T   = LN2 / T_HALF_S;                  // decay constant [1/s]
static constexpr double N_A        = 6.02214076e23;                    // Avogadro [mol^-1]
static constexpr double R_GAS      = 8.31446f;                         // [J/(mol·K)]
static constexpr double M_T_kg_mol = 3.016049e-3;                      // molar mass T [kg/mol]
static constexpr double E_BETA_J   = 5.685e-15;                        // mean beta energy per decay [J] (5.7 keV)

// ─── tritiumBurnFraction ──────────────────────────────────────────────────────
float tritiumBurnFraction(float n_T_m3, float sigma_v_m3s, float confinement_time_s)
{
    // Lawson-type burn fraction:
    //   f_burn = 1 - exp(-n_T * <σv> * τ_conf)  ≈ n_T * <σv> * τ_conf  (for small f)
    //
    // More precisely, in a 50:50 D-T plasma where n_T = n_D = n/2:
    //   f_burn = n_T * <σv> * τ_conf / (1 + n_T * <σv> * τ_conf)
    //   (recycling equilibrium)
    if (confinement_time_s <= 0.0f || sigma_v_m3s <= 0.0f) return 0.0f;

    float reaction_probability = n_T_m3 * sigma_v_m3s * confinement_time_s;
    float f_burn = reaction_probability / (1.0f + reaction_probability);

    return std::clamp(f_burn, 0.0f, 1.0f);
}

// ─── tritiumInventory ─────────────────────────────────────────────────────────
float tritiumInventory(const TritiumRates& r)
{
    return r.production_s
         - r.consumption_s
         - r.decay_s
         - r.permeation_s
         - r.misc_losses_s;
}

float tritiumInventory(float production_s, float consumption_s,
                        float decay_s, float losses_s)
{
    return production_s - consumption_s - decay_s - losses_s;
}

// ─── tritiumDecay ─────────────────────────────────────────────────────────────
DecayResult tritiumDecay(float N_T_initial, double t_s)
{
    DecayResult res{};

    // N(t) = N₀ * exp(-λ * t)
    double N_T = static_cast<double>(N_T_initial) * exp(-LAMBDA_T * t_s);
    N_T        = std::max(N_T, 0.0);

    res.N_T_remaining  = static_cast<float>(N_T);
    res.N_He3_produced = static_cast<float>(N_T_initial - N_T);

    // Activity: A = λ * N [Bq = decays/s]
    res.activity_Bq    = static_cast<float>(LAMBDA_T * N_T);

    // Decay heat: Q = A * E_beta  [W]
    // For 1 gram T: activity = 3.57e14 Bq → heat ≈ 0.324 W/g
    res.decay_heat_W   = static_cast<float>(res.activity_Bq * E_BETA_J);

    return res;
}

// ─── tritiumPermeation ────────────────────────────────────────────────────────

// Permeability parameters: P_s(T) = P_s0 * exp(-E_a / (R*T))  [mol/(m·s·Pa^0.5)]
struct PermData { double P_s0; double E_a_J_mol; };

static PermData getPermData(WallMaterial mat)
{
    switch (mat) {
        case WallMaterial::SS316L:
            // Forcey (1988): P_s0 = 4.0e-7, E_a = 49.7 kJ/mol
            return { 4.0e-7, 49700.0 };
        case WallMaterial::Tungsten:
            // Frauenfelder (1969): very low permeability
            return { 1.6e-7, 96200.0 };
        case WallMaterial::EUROFER97:
            // Reduced activation ferritic-martensitic: similar to SS
            return { 3.0e-7, 48500.0 };
        case WallMaterial::Inconel625:
            // Louthan (1975)
            return { 5.5e-7, 52800.0 };
        default:
            return { 4.0e-7, 49700.0 };
    }
}

PermeationResult tritiumPermeation(float p_upstream_Pa,
                                    float p_downstream_Pa,
                                    float T_wall_K,
                                    float thickness_m,
                                    float area_m2,
                                    WallMaterial material)
{
    PermeationResult res{};
    if (thickness_m < 1e-6f || T_wall_K < 200.0f) return res;

    auto pd = getPermData(material);

    // Permeability: P_s(T) = P_s0 * exp(-E_a / (R*T))  [mol/(m·s·Pa^0.5)]
    double P_s = pd.P_s0 * exp(-pd.E_a_J_mol / (R_GAS * T_wall_K));
    res.permeability_m3s = static_cast<float>(P_s);

    // Sieverts' law flux [mol/(m²·s)]:
    //   J = P_s * (sqrt(p_up) - sqrt(p_down)) / L
    float dp_sqrt = sqrtf(std::max(p_upstream_Pa,   0.0f))
                  - sqrtf(std::max(p_downstream_Pa, 0.0f));
    dp_sqrt = std::max(dp_sqrt, 0.0f);  // no back-permeation here

    res.flux_mol_m2s    = static_cast<float>(P_s) * dp_sqrt / thickness_m;
    res.total_loss_mol_s = res.flux_mol_m2s * area_m2;

    // Activity of permeated T [Bq/s]:
    //   A = (N_A * flux_mol_s) * λ
    res.total_loss_Bq = static_cast<float>(
        res.total_loss_mol_s * N_A * LAMBDA_T);

    return res;
}

// ─── isotopeSeparationWork ────────────────────────────────────────────────────

// Value function: V(x) = (2x - 1) * ln(x / (1-x))
static float valueFunction(float x)
{
    x = std::clamp(x, 1e-6f, 1.0f - 1e-6f);
    return (2.0f * x - 1.0f) * logf(x / (1.0f - x));
}

SeparationResult isotopeSeparationWork(float F_feed_mol_s,
                                        float x_feed,
                                        float x_product,
                                        float x_waste,
                                        float alpha)
{
    SeparationResult res{};
    x_feed    = std::clamp(x_feed,    1e-6f, 1.0f - 1e-6f);
    x_product = std::clamp(x_product, x_feed + 1e-5f, 1.0f - 1e-6f);
    x_waste   = std::clamp(x_waste,   1e-6f, x_feed - 1e-5f);

    // Material balance:
    //   F = P + W
    //   F * x_f = P * x_p + W * x_w
    float denom = x_product - x_waste;
    if (fabsf(denom) < 1e-10f) return res;

    float P_frac = (x_feed - x_waste) / denom;       // P/F
    float W_frac = 1.0f - P_frac;                     // W/F

    res.P_product_mol_s = F_feed_mol_s * P_frac;
    res.W_waste_mol_s   = F_feed_mol_s * W_frac;

    // SWU calculation: SWU/s = P*V(xp) + W*V(xw) - F*V(xf)
    res.SWU_mol_s = res.P_product_mol_s * valueFunction(x_product)
                  + res.W_waste_mol_s   * valueFunction(x_waste)
                  - F_feed_mol_s        * valueFunction(x_feed);
    res.SWU_mol_s = std::max(res.SWU_mol_s, 0.0f);

    // Minimum stages from separation factor:
    // N_min = ln((xp/(1-xp)) / (xw/(1-xw))) / ln(α)
    if (alpha > 1.0f + 1e-6f) {
        float enrichment_ratio = (x_product / (1.0f - x_product))
                               / (x_waste   / (1.0f - x_waste));
        res.N_stages = logf(enrichment_ratio) / logf(alpha);
    }

    // Theoretical minimum energy: kT * ln(α) per separation event
    // Practical distillation / palladium membrane plants: ~1000 kJ/mol T
    res.energy_kJ_mol = (res.SWU_mol_s > 0.0f && res.P_product_mol_s > 0.0f)
                      ? 1000.0f * res.SWU_mol_s / res.P_product_mol_s
                      : 0.0f;

    return res;
}

} // namespace TritiumCyclePhysics