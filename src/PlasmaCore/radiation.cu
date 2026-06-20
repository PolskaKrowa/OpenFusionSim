//
// radiation.cu
// Plasma radiation losses: bremsstrahlung, synchrotron, line radiation.
//
//  All three loss channels are computed per-cell from local (n_e, T_e, B, Z_eff)
//  and accumulated into a single q_rad_voxel grid that the thermal-hydraulics
//  module reads.  This is the missing piece that determines whether a fusion
//  plasma can actually ignite — without radiation losses, the simulation
//  always reports Q → ∞, which is unphysical.
//
//  References:
//    [1] Trubnikov, Rev. Plasma Phys. 7, 345 (1979) — relativistic bremsstrahlung.
//    [2] Albajar & Bornatici, Nucl. Fusion 41, 895 (2001) — synchrotron scaling
//        for optically-thick reactor plasmas (replaces the single-pass Trubnikov
//        formula which underestimates reabsorption).
//    [3] NRL Plasma Formulary 2023, p. 57 — Bremsstrahlung coefficient.
//    [4] Post, J. Nucl. Mater. 220-222, 143 (1995) — impurity cooling factors.
//

#include "types.cuh"
#include <math.h>

// ─── Bremsstrahlung (Trubnikov, with relativistic correction) ─────────────────
//
//  P_brem [W/m³] = C_br · n_e² · Z_eff · √(T_e[keV]) · ḡ_B · (1 + rel. corr.)
//
//  Non-relativistic (NRL):   C_br = 5.34e-37 W·m³/√keV
//  ḡ_B (velocity-averaged Gaunt factor) ≈ 1.3 for 1-30 keV (Sauter 1998).
//
//  Trubnikov relativistic correction series (in T_e / m_e c², dimensionless):
//      corr = 1 + 2.76·τ + 2.04·τ² + ...
//      τ = T_e / (m_e c²) = T_e[keV] / 511.0
//
//  At T_e = 50 keV → τ = 0.098 → corr = 1.31.  At 100 keV → corr = 1.70.
//
__device__ __forceinline__
float bremsstrahlungPower(float n_e_m3, float T_e_keV, float Z_eff)
{
    if (n_e_m3 < 1.0f || T_e_keV < 0.01f) return 0.0f;

    constexpr float C_br    = 5.34e-37f;
    constexpr float g_avg   = 1.30f;             // Gaunt factor for 1-30 keV
    constexpr float m_e_c2  = 511.0f;            // keV

    float tau = T_e_keV / m_e_c2;
    float relativistic_corr = 1.0f + 2.76f * tau + 2.04f * tau * tau;

    return C_br * n_e_m3 * n_e_m3 * Z_eff * sqrtf(T_e_keV)
         * g_avg * relativistic_corr;
}

// ─── Synchrotron (electron cyclotron) — Albajar 2001 ──────────────────────────
//
//  P_sync [W/m³] ≈ C_s · n_e · B² · T_e^{2.5} · F_reabsorb
//
//  This is the *leading-order* Trubnikov single-pass form, scaled by the
//  Albajar 2001 reabsorption factor F(τ_opt) for optically-thick reactor
//  plasmas.  F → 1 for optically thin (small machines); F → ~0.5 for reactor-
//  grade plasmas with R_w ≈ 0.7 wall reflectivity.
//
//  We use the simpler Trubnikov leading form here; the full Albajar fit
//  requires (R_major, a_minor, R_wall) which are not available in the
//  per-cell kernel.
//
//  IMPORTANT: C_s = 1.6e-21 is dimensioned for n_e in m⁻³ and T_e in keV.
//  Older NRL editions (and earlier versions of this file) carried
//  6.21e-17, which is the n_e-in-cm⁻³ single-pass constant — using it
//  with n_e in m⁻³ overestimates reactor sync losses by ~4 orders of
//  magnitude (~30 GW at ITER conditions instead of ~1 MW) and is the
//  root cause of the spurious overtemperature SCRAM bug.
//
//  Reference: Albajar & Bornatici, Nucl. Fusion 41, 895 (2001), Eq. 12.
//
__device__ __forceinline__
float synchrotronPower(float n_e_m3, float T_e_keV, float B_T, float R_wall)
{
    if (n_e_m3 < 1.0f || T_e_keV < 0.1f || B_T < 0.01f) return 0.0f;

    constexpr float C_s = 1.6e-21f;   // W·m³/(T²·keV^{2.5})  [Albajar 2001]

    // Wall reflectivity correction (Albajar Eq. 24, simplified)
    // F ≈ (1 - R_wall) · [1 + 0.12·T_e]^{-1}    ... very rough
    float F = (1.0f - R_wall) / (1.0f + 0.12f * T_e_keV);

    // T_e^{2.5} scaling (Albajar) instead of the naive T_e (Trubnikov single-pass)
    // — dominates above 30 keV.
    return C_s * n_e_m3 * B_T * B_T * powf(T_e_keV, 2.5f) * F;
}

// ─── Line radiation (corona model) ────────────────────────────────────────────
//
//  P_line [W/m³] = n_e · Σ_i n_i · L_i(T_e)
//
//  For a fully-ionised D-T plasma with intrinsic W impurity:
//      L_W(T_e) peaks at ~1e-34 W·m³ near 10-20 keV (ADAS ADF11)
//      drops sharply above 20 keV (fully stripped).
//
//  We use a simple parameterisation adequate for the simulator:
//      L_line = n_e · n_imp · C_line(T_e, Z_imp)
//
//  Where C_line is a table lookup or analytical fit.  For W at fusion T:
//      C_line ≈ 1e-34 · exp(-(T - 15)^2 / 100)   [W·m³]   (rough Gaussian)
//
__device__ __forceinline__
float lineRadiationPower(float n_e_m3, float n_imp_m3,
                         float T_e_keV, int Z_imp)
{
    if (n_imp_m3 < 1.0f || T_e_keV < 0.1f) return 0.0f;

    // Tungsten (Z=74) cooling factor — peak near 15 keV, falls above 20
    float L_W;
    if (Z_imp == 74) {
        float dT = T_e_keV - 15.0f;
        L_W = 1.0e-34f * expf(-dT * dT / 100.0f);
    } else if (Z_imp == 6) {
        // Carbon (common low-Z impurity): broader peak ~10 eV, falls fast above 1 keV
        L_W = 1.0e-35f * expf(-T_e_keV * 5.0f);
    } else if (Z_imp == 26) {
        // Iron: peak near 30 eV, partial ionization above 1 keV
        L_W = 1.0e-35f * expf(-T_e_keV * 2.0f);
    } else {
        // Generic low-Z: assume fully stripped above 5 keV → negligible line radiation
        L_W = (T_e_keV < 5.0f) ? 1.0e-35f * expf(-T_e_keV) : 0.0f;
    }

    return n_e_m3 * n_imp_m3 * L_W;
}

// ─── Effective Z from impurity mix ────────────────────────────────────────────
//
//  Z_eff = Σ_i (n_i · Z_i²) / n_e
//
//  This is used both for the bremsstrahlung scaling (Z_eff) and for the
//  collisional frequency correction (also Z_eff).  We compute it from the
//  impurity species densities.
//
__device__ __forceinline__
float effectiveZ(float n_D, float n_T, float n_He, float n_imp, int Z_imp)
{
    float n_e = n_D + n_T + 2.0f * n_He + Z_imp * n_imp;   // quasi-neutrality
    if (n_e < 1.0f) return 1.0f;
    float numerator = n_D + n_T + 4.0f * n_He + Z_imp * Z_imp * n_imp;
    return numerator / n_e;
}

// ─── Radiation Kernel ─────────────────────────────────────────────────────────
//
//  One thread per cell.  Reads (n_e, T_e, B) from per-cell diagnostics,
//  computes the three radiation contributions, and writes the total to
//  the heat deposition map (q_dot_voxel).
//
//  Note: the radiation power is removed from the *electron* population
//  (via a corresponding drag in the momentum equation).  In this 0D-per-cell
//  form we simply deposit it as heat on the first wall — physically the
//  photons escape the plasma and are absorbed by the wall.
//
__global__ void computeRadiationLosses(
    float* __restrict__ q_dot_voxel,         // heat deposition map (ADD)
    float* __restrict__ q_rad_voxel,         // radiation-only map (for diag)
    const float* __restrict__ n_e_per_cell,
    const float* __restrict__ T_e_per_cell_keV,
    const float* __restrict__ B_per_cell_T,
    const float* __restrict__ n_imp_per_cell,   // may be null → 0 impurity
    const int*   __restrict__ Z_imp_per_cell,   // may be null → Z=0
    float R_wall,
    int N_cells)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;

    float n_e   = n_e_per_cell[tid];
    float T_e   = T_e_per_cell_keV[tid];
    float B     = B_per_cell_T[tid];
    // ── Null-safe impurity reads ─────────────────────────────────────────────
    // Callers that don't track impurities (e.g. the bridge during plasma
    // startup) pass nullptr here.  Dereferencing null on the GPU is a warp
    // MMU fault that poisons the whole context — so we guard with explicit
    // null checks rather than relying on the caller to pass zeroed arrays.
    float n_imp = (n_imp_per_cell  != nullptr) ? n_imp_per_cell[tid]  : 0.0f;
    int   Z_imp = (Z_imp_per_cell  != nullptr) ? Z_imp_per_cell[tid]  : 0;

    if (n_e < 1.0f || T_e < 0.01f) {
        q_rad_voxel[tid] = 0.0f;
        return;
    }

    // Z_eff from local impurity mix (rough — assumes 50:50 D-T with He ash)
    float n_D  = 0.5f * n_e * 0.9f;
    float n_T  = 0.5f * n_e * 0.9f;
    float n_He = n_e * 0.05f;
    float Z_eff = effectiveZ(n_D, n_T, n_He, n_imp, Z_imp);

    float P_brem = bremsstrahlungPower(n_e, T_e, Z_eff);
    float P_sync = synchrotronPower(n_e, T_e, B, R_wall);
    float P_line = lineRadiationPower(n_e, n_imp, T_e, Z_imp);

    float P_total = P_brem + P_sync + P_line;
    q_rad_voxel[tid] = P_total;
    atomicAdd(&q_dot_voxel[tid], P_total);
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchRadiationLosses(
    float* q_dot_voxel,
    float* q_rad_voxel,
    const float* n_e_per_cell,
    const float* T_e_per_cell_keV,
    const float* B_per_cell_T,
    const float* n_imp_per_cell,
    const int*   Z_imp_per_cell,
    float R_wall,
    int N_cells,
    cudaStream_t stream)
{
    constexpr int BLOCK = 256;
    int grid = (N_cells + BLOCK - 1) / BLOCK;
    computeRadiationLosses<<<grid, BLOCK, 0, stream>>>(
        q_dot_voxel, q_rad_voxel,
        n_e_per_cell, T_e_per_cell_keV, B_per_cell_T,
        n_imp_per_cell, Z_imp_per_cell,
        R_wall, N_cells);
}
