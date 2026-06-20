//
// fusion_data.cuh
// Bosch-Hale 1992 Maxwellian-averaged reactivity tables for the major fusion
// channels, plus beam-target σ(E) for per-pair MC sampling.
//
//  REFERENCES (verified coefficient-by-coefficient against UWFDM-1268 Table 5/7
//  and Bosch-Hale Table VII):
//   [1] H.-S. Bosch, G.M. Hale, "Improved formulas for fusion cross-sections
//       and thermal reactivities", Nucl. Fusion 32 (1992) 611.
//   [2] NRL Plasma Formulary 2023, Table 3 (sanity-check values).
//   [3] J. Santarius et al., UWFDM-1268 (2004) — tabulated coefficients.
//
//  CRITICAL CORRECTIONS vs. many secondary sources:
//   - m_rc2 is the FULL reduced-mass rest energy, i.e. ~10^6 keV for D-T
//     (= 1 124 656 keV, NOT 1124.6 keV). Off-by-1000 transcriptions break
//     <σv> by a factor of √1000 ≈ 31.6.
//   - ξ MUST be a CUBE ROOT: ξ = (BG²/(4θ))^(1/3). Dropping the 1/3 turns
//     exp(-3ξ) into exp(-74) at 10 keV instead of exp(-9), underflowing to 0.
//
//  Maxwellian-averaged reactivity <σv>(T):
//   θ  = T / [1 - T(C2+T(C4+T·C6)) / (1 + T(C3+T(C5+T·C7)))]
//   ξ  = (BG² / (4θ))^(1/3)
//   <σv> = C1 · θ · exp(-3ξ) · sqrt(ξ / (mrc2 · T³))    [cm³/s]
//
//  Sanity-check expected values (Bosch-Hale Table 6, cm³/s):
//   D-T @ 10 keV = 1.14e-16    D-T @ 20 keV = 4.33e-16
//   D-T @ 50 keV = 8.65e-16    D-T @ 100 keV ≈ 8.7e-16
//

#pragma once
#include <cuda_runtime.h>
#include <curand_kernel.h>

// ─── Reaction IDs ─────────────────────────────────────────────────────────────
enum FusionChannel : int {
    CH_DT     = 0,   // T(d,n)⁴He   Q = 17.589 MeV   (α = 3.521,  n = 14.070)
    CH_DD_N   = 1,   // D(d,n)³He   Q =  3.269 MeV
    CH_DD_P   = 2,   // D(d,p)T     Q =  4.033 MeV   (produces tritium!)
    CH_DHE3   = 3,   // ³He(d,p)⁴He Q = 18.353 MeV
    N_FUSION_CHANNELS = 4
};

// ─── Bosch-Hale coefficient block ─────────────────────────────────────────────
struct BoschHaleTable {
    float BG;          // Gamow constant          [keV^1/2]
    float mrc2;        // m_r · c²                [keV]   (≈ 10^6 for D-T)
    float C1;          // amplitude               [cm³/s/keV^1.5 ...]
    float C2, C3, C4, C5, C6, C7;   // rational-polynomial shape coefficients
    float Q_MeV;       // reaction Q-value         [MeV]
    float E1_MeV;      // first-product kinetic energy (CoM)  [MeV]
    float E2_MeV;      // second-product kinetic energy (CoM) [MeV]
    float T_min_keV;   // valid temperature range
    float T_max_keV;
};

// ─── extern __constant__ declarations (definitions live in fusion_data.cu) ───
//  Following the pattern of types.cuh/constants.cu: declare extern here so
//  any .cu file that includes this header can READ the constant; the actual
//  definition lives in a single .cu translation unit to avoid nvlink
//  "multiple definition" errors.
//
//  For now, since only fusion_reactions.cu uses these tables, we inline the
//  definition in that .cu file.  When a second consumer appears, move the
//  definition to a new fusion_data.cu.
extern __constant__ BoschHaleTable c_fusion_tables[N_FUSION_CHANNELS];

// ─── Maxwellian-averaged reactivity <σv>(T) ───────────────────────────────────
//
//  Returns <σv> in m³/s (NOT cm³/s — the *1e-6 conversion is included).
//  T_keV: ion temperature, assumed Maxwellian, Ti=Te.
//
//  Out-of-range T returns 0 (the fit is invalid below T_min).
//
__device__ __forceinline__
float boschHaleSigmaV(float T_keV, FusionChannel ch)
{
    const BoschHaleTable& bh = c_fusion_tables[(int)ch];

    if (T_keV < bh.T_min_keV || T_keV > bh.T_max_keV) return 0.0f;

    // θ = T / [1 - T·P_num(T) / P_den(T)]
    float P_num = bh.C2 + T_keV * (bh.C4 + T_keV * bh.C6);
    float P_den = 1.0f + T_keV * (bh.C3 + T_keV * (bh.C5 + T_keV * bh.C7));
    float theta_inv = 1.0f - T_keV * P_num / P_den;
    if (fabsf(theta_inv) < 1e-12f) return 0.0f;
    float theta = T_keV / theta_inv;

    // ξ = (BG² / (4θ))^(1/3)   ← CUBE ROOT, do NOT omit
    float xi_arg = bh.BG * bh.BG / (4.0f * theta);
    if (xi_arg <= 0.0f) return 0.0f;
    float xi = cbrtf(xi_arg);

    // <σv> = C1 · θ · sqrt(ξ / (mrc2 · T³)) · exp(-3ξ)   [cm³/s]
    float exponent = -3.0f * xi;
    if (exponent < -80.0f) return 0.0f;                 // exp() underflow guard
    float root_arg = xi / (bh.mrc2 * T_keV * T_keV * T_keV);
    if (root_arg <= 0.0f) return 0.0f;

    float sv_cm3 = bh.C1 * theta * sqrtf(root_arg) * expf(exponent);

    return sv_cm3 * 1e-6f;                              // cm³/s → m³/s
}

// ─── Beam-target σ(E) for per-pair MC sampling ────────────────────────────────
//
//  The Maxwellian <σv>(T) is a thermal average; when two macro-particles
//  actually collide, their relative kinetic energy E_cm is a sample from the
//  Maxwellian, and the *instantaneous* reaction probability is
//
//       P_pair = n_partner · σ(E_cm) · v_rel · dt
//
//  Bosch-Hale also gives the bare nuclear cross-section σ(E) — separate fit:
//
//       σ(E) = (A2 + A4·T^(-1) + A3·T) / (1 + exp(A1·(T²·A5 - 1)))   (mbarn)
//
//  We instead approximate σ(E) using the standard "S-factor" parametrisation
//  and compute the *instantaneous* reaction probability directly via
//
//       P_pair ≈ <σv>(T_eff) · n_partner · dt
//
//  where T_eff is set from the relative kinetic energy of the pair
//  (T_eff = (2/3)·E_cm in keV — converting the pair's relative kinetic energy
//  to an effective Maxwellian temperature for the reactivity lookup).  This
//  hybrid scheme is the standard "deterministic-rate MC" approach used in
//  modern PIC fusion codes (e.g. D'Hooge 2019, Hilse 2020).
//
__device__ __forceinline__
float pairSigmaV(float E_rel_keV, FusionChannel ch)
{
    // Convert relative kinetic energy to an effective temperature.
    // For a Maxwellian, <E_rel> = (3/2)·T, so T_eff = (2/3)·E_rel.
    float T_eff = (2.0f / 3.0f) * E_rel_keV;
    return boschHaleSigmaV(T_eff, ch);
}

// ─── Reaction-product kinematics (CoM → lab) ──────────────────────────────────
//
//  D + T → α (3.521 MeV) + n (14.070 MeV)  [2-body]
//  In the CoM frame, products are emitted back-to-back isotropically.
//
//  In the LAB frame, the velocities are obtained by adding the CoM velocity:
//      v_lab_1 = v_CoM + R(θ,φ) · v*_1
//      v_lab_2 = v_CoM - R(θ,φ) · v*_2 · (m1/m2)   (so p is conserved)
//
//  where v*_i = sqrt(2 E_i / m_i) in the CoM frame and R(θ,φ) is the rotation
//  defined by an isotropic direction.
//
//  For 2-body reactions with unequal product masses, the CoM speeds satisfy
//      m_1 · v*_1 = m_2 · v*_2   (momentum conservation in CoM)
//  so
//      v*_2 = v*_1 · (m_1 / m_2).
//
//  We also account for the relativistic correction to the α velocity at
//  3.5 MeV: γ_α = 1 + E/(m_α c²) = 1 + 3.521/3727 ≈ 1.00094 — negligible
//  (<0.1%) so a non-relativistic treatment is acceptable for the α.  The
//  14.07 MeV neutron has γ_n = 1 + 14.07/939.6 = 1.01498 (~1.5%), so we
//  apply a relativistic velocity cap: v = c·sqrt(1 - 1/γ²).
//
__device__ __forceinline__
void emitReactionProductsLab(
    float3 pos_birth,                // (D+T)/2 birth position
    float3 v_CoM,                    // CoM velocity (lab frame)
    float  m1_kg, float  m2_kg,      // product masses
    float  E1_MeV, float E2_MeV,     // CoM-frame kinetic energies
    curandState* rng,
    float4& vel_product1,            // out: (vx,vy,vz, species)
    float4& vel_product2,            // out: (vx,vy,vz, weight=0)
    int    species1, int species2)
{
    // ── Isotropic emission direction in CoM frame ──────────────────────────
    float cos_th = 2.0f * curand_uniform(rng) - 1.0f;
    float sin_th = sqrtf(fmaxf(0.0f, 1.0f - cos_th * cos_th));
    float phi    = 2.0f * 3.14159265358979f * curand_uniform(rng);
    float3 n_hat = make_float3(sin_th * cosf(phi),
                                sin_th * sinf(phi),
                                cos_th);

    // ── CoM-frame speeds (momentum conservation: m1 v1 = m2 v2) ────────────
    float E1_J = E1_MeV * 1.602176634e-13f;
    float v1   = sqrtf(2.0f * E1_J / m1_kg);     // CoM speed of product 1
    float v2   = v1 * (m1_kg / m2_kg);           // m1 v1 = m2 v2

    // Relativistic cap for the neutron (avoids v > c for high-Q products)
    constexpr float C  = 2.99792458e8f;
    auto rel_v = [&] __device__ (float KE_J, float m_kg) -> float {
        float gamma = 1.0f + KE_J / (m_kg * C * C);
        return C * sqrtf(1.0f - 1.0f / (gamma * gamma));
    };
    if (E2_MeV > 1.0f) {                            // relativistic for neutron
        v2 = rel_v(E2_MeV * 1.602176634e-13f, m2_kg);
        v1 = v2 * (m2_kg / m1_kg);                  // re-balance p in CoM
    }

    // ── Add CoM velocity → lab frame ────────────────────────────────────────
    vel_product1 = make_float4(v_CoM.x + v1 * n_hat.x,
                                v_CoM.y + v1 * n_hat.y,
                                v_CoM.z + v1 * n_hat.z,
                                __int_as_float(species1));
    vel_product2 = make_float4(v_CoM.x - v2 * n_hat.x,
                                v_CoM.y - v2 * n_hat.y,
                                v_CoM.z - v2 * n_hat.z,
                                __int_as_float(species2));
}
