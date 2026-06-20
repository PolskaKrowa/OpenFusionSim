//
// src/ThermalHydraulics/MoltenSaltPhysics.cpp
//

#include "MoltenSaltPhysics.h"
#include <cmath>
#include <algorithm>

namespace MoltenSaltPhysics {

static constexpr float G = 9.81f;  // gravitational acceleration [m/s²]

// ─── heatGenerationInBlanket ──────────────────────────────────────────────────
BlanketHeatGenResult heatGenerationInBlanket(float neutron_flux_m2s,
                                              float Li6_density_m3,
                                              float blanket_volume_m3,
                                              SaltType salt)
{
    BlanketHeatGenResult res{};

    // Li-6(n,t)α  cross-section at 14 MeV: σ ≈ 0.05 barn = 5e-30 m²
    // Q-value of the reaction: Q = 4.78 MeV = 7.66e-13 J
    constexpr float sigma_Li6_m2 = 5.0e-30f;
    constexpr float Q_Li6_J      = 7.66e-13f;

    float rxn_rate = neutron_flux_m2s * Li6_density_m3 * sigma_Li6_m2; // [reactions/m³/s]
    res.q_Li6_reaction_W_m3 = rxn_rate * Q_Li6_J;

    // Gamma heating from neutron moderation: roughly 30 % of total nuclear heating
    // This depends strongly on the blanket composition and geometry.
    // For FLiBe with Be(n,2n) multiplier: effective energy multiplication ~ 1.2
    float mult = (salt == SaltType::FLiBe || salt == SaltType::FLiBe_Li7) ? 1.25f : 1.10f;
    res.q_gamma_W_m3   = res.q_Li6_reaction_W_m3 * (mult - 1.0f);
    res.q_dot_W_m3     = res.q_Li6_reaction_W_m3 * mult;
    res.total_power_MW = res.q_dot_W_m3 * blanket_volume_m3 * 1e-6f;

    return res;
}

// ─── dittusBoelterNusselt ─────────────────────────────────────────────────────
float dittusBoelterNusselt(float Re, float Pr, bool heating)
{
    // Clamp to validity range
    Re = std::max(Re, 10000.0f);
    Pr = std::clamp(Pr, 0.6f, 160.0f);

    float n = heating ? 0.4f : 0.3f;
    return 0.023f * powf(Re, 0.8f) * powf(Pr, n);
}

// ─── moltenSaltViscosity ──────────────────────────────────────────────────────
float moltenSaltViscosity(float T_K, SaltType salt)
{
    // Arrhenius model: μ = A * exp(B / T)
    float A, B;
    switch (salt) {
        case SaltType::FLiBe:
        case SaltType::FLiBe_Li7:
            // Cantor (1973): valid 733–1073 K
            A = 1.16e-4f;
            B = 3755.0f;
            break;
        case SaltType::FLiNaK:
            // Janz (1988): valid 730–1570 K
            A = 4.0e-5f;
            B = 4170.0f;
            break;
        default:
            A = 1.16e-4f;
            B = 3755.0f;
    }
    T_K = std::max(T_K, 700.0f); // below melting point → not physical
    return A * expf(B / T_K);
}

// ─── moltenSaltThermalConductivity ────────────────────────────────────────────
float moltenSaltThermalConductivity(float T_K, SaltType salt)
{
    switch (salt) {
        case SaltType::FLiBe:
        case SaltType::FLiBe_Li7:
            // Williams (2006) polynomial fit: valid 733–1073 K
            return 0.629697f + 5.0e-4f * T_K;
        case SaltType::FLiNaK:
            // Approximately constant (Janz 1988)
            return 0.60f;
        default:
            return 0.629697f + 5.0e-4f * T_K;
    }
}

// ─── moltenSaltHeatCapacity ───────────────────────────────────────────────────
float moltenSaltHeatCapacity(float T_K, SaltType salt)
{
    // FLiBe: Sohal (2010), ~constant 700–1000 K
    (void)T_K;
    switch (salt) {
        case SaltType::FLiBe:
        case SaltType::FLiBe_Li7:
            return 2415.78f;    // J/(kg·K)
        case SaltType::FLiNaK:
            return 1905.0f;     // J/(kg·K)
        default:
            return 2415.78f;
    }
}

// ─── bulkTemperatureRise ──────────────────────────────────────────────────────
float bulkTemperatureRise(float q_dot_W, float mdot_kg_s, float cp_J_kgK)
{
    if (mdot_kg_s < 1e-6f || cp_J_kgK < 1.0f) return 0.0f;
    return q_dot_W / (mdot_kg_s * cp_J_kgK);
}

// ─── pressureDrop ─────────────────────────────────────────────────────────────
PressureDropResult pressureDrop(float L_m, float D_m, float rho_kg_m3,
                                 float v_m_s, float f_darcy, float mu_Pa_s)
{
    PressureDropResult res{};
    if (D_m < 1e-6f) return res;

    res.v_m_s = v_m_s;
    res.Re    = rho_kg_m3 * v_m_s * D_m / std::max(mu_Pa_s, 1e-10f);

    // If f_darcy not supplied (== 0), compute via Colebrook-White (smooth pipe)
    if (f_darcy <= 0.0f) {
        if (res.Re < 2300.0f) {
            // Laminar: f = 64 / Re
            res.f_darcy = 64.0f / std::max(res.Re, 1.0f);
        } else {
            // Turbulent, smooth pipe: Petukhov equation
            float lnRe  = logf(std::max(res.Re, 1.0f));
            float denom = 0.790f * lnRe - 1.64f;
            res.f_darcy = 1.0f / (denom * denom);
        }
    } else {
        res.f_darcy = f_darcy;
    }

    // Darcy-Weisbach: ΔP = f * (L/D) * (ρ v² / 2)
    res.dP_Pa       = res.f_darcy * (L_m / D_m) * (rho_kg_m3 * v_m_s * v_m_s / 2.0f);

    // Pumping power: P_pump = Q_vol * ΔP = A * v * ΔP
    float A_pipe    = 3.14159f * D_m * D_m / 4.0f;
    res.pump_power_W = A_pipe * v_m_s * res.dP_Pa;

    return res;
}

// ─── HartmannNumber ───────────────────────────────────────────────────────────
float HartmannNumber(float B_T, float L_m, float sigma_S_m, float mu_Pa_s)
{
    if (mu_Pa_s < 1e-12f) return 0.0f;
    return B_T * L_m * sqrtf(sigma_S_m / mu_Pa_s);
}

// ─── magnetohydrodynamicDrag ──────────────────────────────────────────────────
MHDResult magnetohydrodynamicDrag(float sigma_S_m, float B_T, float v_m_s,
                                   float L_char_m, float mu_Pa_s, float rho_kg_m3)
{
    MHDResult res{};
    res.Ha = HartmannNumber(B_T, L_char_m, sigma_S_m, mu_Pa_s);

    // MHD pressure drop (Hartmann regime, Ha >> 1):
    //   ΔP_MHD ≈ Ha² * μ * v / L²   (for rectangular duct, parallel walls)
    //
    // Full expression for Hartmann flow in rectangular duct:
    //   ΔP/L = (σ * B² * v) / (1 + 1/Ha)
    //   For Ha >> 1: ΔP/L → σ * B² * v  (fully MHD dominated)
    if (res.Ha < 0.1f) {
        res.dP_MHD_Pa = 0.0f;
    } else {
        float dP_per_L = sigma_S_m * B_T * B_T * v_m_s / (1.0f + 1.0f / res.Ha);
        res.dP_MHD_Pa  = dP_per_L * L_char_m; // approximate path length
    }

    // Velocity profile: MHD flattens the profile
    // At Ha >> 1: flat core with thin Hartmann layers (thickness δ ~ L/Ha)
    res.velocity_profile = std::clamp(res.Ha / (res.Ha + 10.0f), 0.0f, 1.0f);

    // Turbulence suppression: k_eff reduces below molecular k at high Ha
    // (MHD quenches turbulent eddies perpendicular to B)
    float k_mol = 0.629697f + 5.0e-4f * 800.0f; // typical FLiBe k at ~800 K
    float suppression = 1.0f / (1.0f + 0.01f * res.Ha);
    res.effective_k_W_mK = k_mol * (0.1f + 0.9f * suppression); // retain at least molecular k

    (void)rho_kg_m3;
    return res;
}

// ─── naturalConvectionFlow ────────────────────────────────────────────────────
NatConvResult naturalConvectionFlow(float rho_kg_m3, float beta_K, float deltaT_K,
                                     float L_m, float k_W_mK, float cp_J_kgK,
                                     float mu_Pa_s, float decay_power_W)
{
    NatConvResult res{};

    // Thermal diffusivity: α = k / (ρ * cp)
    float alpha_m2s = k_W_mK / (rho_kg_m3 * cp_J_kgK + 1e-10f);
    float nu_m2s    = mu_Pa_s / (rho_kg_m3 + 1e-10f);   // kinematic viscosity

    // Rayleigh number: Ra = g * β * ΔT * L³ / (ν * α)
    res.Ra = G * beta_K * deltaT_K * L_m * L_m * L_m / (nu_m2s * alpha_m2s + 1e-30f);

    // Nusselt number correlation (vertical cavity, Ra > 10^4):
    if (res.Ra < 1.0f) {
        res.Nu = 1.0f;
    } else if (res.Ra < 1e9f) {
        // Transition: Nu = 0.59 * Ra^(1/4)
        res.Nu = 0.59f * powf(res.Ra, 0.25f);
    } else {
        // Turbulent natural convection: Nu = 0.15 * Ra^(1/3)
        res.Nu = 0.15f * powf(res.Ra, 1.0f / 3.0f);
    }

    // Heat removed: Q = Nu * k * A * ΔT / L
    // Assume loop cross-section area ~ L² (rough square geometry)
    float A_loop = L_m * L_m;
    res.Q_removed_W = res.Nu * k_W_mK * A_loop * deltaT_K / L_m;

    // Buoyancy velocity (order of magnitude):  v ~ sqrt(g * β * ΔT * L)
    res.flow_velocity_ms = sqrtf(G * beta_K * deltaT_K * L_m);

    res.adequate_cooling = (res.Q_removed_W >= decay_power_W);
    return res;
}

} // namespace MoltenSaltPhysics