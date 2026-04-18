//
// fusion_reactions.cu
// Monte Carlo sampling of D-T fusion reactions.
//
//  For each D-T particle pair in the same cell:
//    1. Compute centre-of-mass energy E_cm from relative velocity.
//    2. Look up <σv> from the Bosch-Hale table in constant memory.
//    3. Sample reaction probability P = n_T * <σv> * dt.
//    4. If reaction occurs (curand_uniform < P), spawn an alpha (3.5 MeV)
//       and a neutron (14.1 MeV) at the collision site.
//    5. Alpha particles are injected back into the ion particle arrays.
//    6. Neutrons are handed off to the separate MC neutron transport kernel.
//
//  Bosch-Hale table: 9-coefficient rational polynomial fit to <σv> vs T[keV].
//  See H.-S. Bosch & G.M. Hale, Nucl. Fusion 32, 611 (1992).
//

#include "types.cuh"
#include <math.h>

// ─── Bosch-Hale D-T coefficient table ─────────────────────────────────────────
//  <σv>_DT in units of m^3/s.  Valid range: 0.2–100 keV.
struct BoschHaleCoeffs {
    float BG;    // Gamow constant
    float mrc2;  // reduced mass * c^2 [keV]
    float C[7];  // rational polynomial numerator coefficients
    float B[4];  // denominator coefficients
};

__constant__ BoschHaleCoeffs c_DT_BH = {
    34.3827f,                                        // BG [keV^1/2]
    1124656.0f,                                      // mrc2 [keV]
    {6.661e-7f, 6.371e-2f, -3.136e-2f,              // C[0..2]
     4.343e-3f, -3.769e-3f, 5.372e-4f, -2.926e-5f}, // C[3..6]
    {0.0f, -2.291e-3f, 0.0f, 0.0f}                  // B[0..3] (simplified)
};

// ─── Bosch-Hale <σv> evaluation ───────────────────────────────────────────────
//  T_keV: ion temperature in keV (approximated from relative KE)
//  returns <σv> in m^3/s
__device__ float boschHaleSigmaV(float T_keV)
{
    if (T_keV < 0.2f || T_keV > 100.0f) return 0.0f;

    const BoschHaleCoeffs& bh = c_DT_BH;
    float theta_inv = 1.0f - T_keV * (bh.C[1] + T_keV * (bh.C[3] + T_keV * bh.C[5]))
                           / (1.0f + T_keV * (bh.C[2] + T_keV * (bh.C[4] + T_keV * bh.C[6])));
    // Guard against bad denominator
    if (fabsf(theta_inv) < 1e-10f) return 0.0f;
    float theta = T_keV / theta_inv;

    float xi = cbrtf(bh.BG * bh.BG / (4.0f * theta));

    // <σv> = C[0] * theta * sqrt(xi / (mrc2 * T_keV^3)) * exp(-3 * xi)  [cm^3/s]
    float exponent = -3.0f * xi;
    if (exponent < -80.0f) return 0.0f; // underflow guard

    float sigma_v_cm3 = bh.C[0] * theta
                      * sqrtf(xi / (bh.mrc2 * T_keV * T_keV * T_keV))
                      * expf(exponent);

    return sigma_v_cm3 * 1e-6f; // cm^3/s → m^3/s
}

// ─── D-T Reaction Energy partition ────────────────────────────────────────────
//  Q = 17.59 MeV total: alpha = 3.52 MeV, neutron = 14.07 MeV
//  Products emitted isotropically in CoM frame.
__device__ void dtReactionProducts(
    float3 posD, float3 posT,
    curandState* rng,
    ReactionProduct& out)
{
    // Birth position: midpoint of D-T pair
    out.pos = make_float4(0.5f * (posD.x + posT.x),
                          0.5f * (posD.y + posT.y),
                          0.5f * (posD.z + posT.z), 1.0f);

    // Isotropic emission direction in CoM frame
    float cos_th = 2.0f * curand_uniform(rng) - 1.0f;
    float sin_th = sqrtf(fmaxf(0.0f, 1.0f - cos_th * cos_th));
    float phi    = 2.0f * 3.14159265f * curand_uniform(rng);
    float3 dir = make_float3(sin_th * cosf(phi), sin_th * sinf(phi), cos_th);

    // Alpha: 3.52 MeV kinetic energy
    //  KE = 0.5 * m * v^2 → v = sqrt(2 KE / m)
    float KE_alpha_J   = 3.52e6f * 1.60217663e-19f;
    float m_alpha      = 4.0f * 1.67262192e-27f;
    float v_alpha      = sqrtf(2.0f * KE_alpha_J / m_alpha);
    out.vel_alpha      = make_float4(v_alpha * dir.x,
                                     v_alpha * dir.y,
                                     v_alpha * dir.z,
                                     __int_as_float(3)); // species=3 → alpha

    // Neutron: 14.07 MeV, emitted opposite the alpha in CoM
    float KE_n_J    = 14.07e6f * 1.60217663e-19f;
    float m_neutron = 1.67492750e-27f;
    float v_neutron = sqrtf(2.0f * KE_n_J / m_neutron);
    out.vel_neutron = make_float4(-v_neutron * dir.x,
                                  -v_neutron * dir.y,
                                  -v_neutron * dir.z, 0.0f);
    out.active = 1;
}

// ─── Fusion Reaction Sampling Kernel ─────────────────────────────────────────
//
//  Called once per timestep after particle sort.
//  Deuterium and tritium must be separated into their own sub-arrays
//  (pos_D/vel_D and pos_T/vel_T) OR species filtering applied inside kernel.
//  Here we take the filtered-array approach for clarity.
//
//  products output array must be pre-allocated to at least N_D (worst case
//  every D reacts).  products[i].active flags which slots are filled.
//
__global__ void sampleFusionReactions(
    const float4* __restrict__ pos_D,
    const float4* __restrict__ vel_D,
    const float4* __restrict__ pos_T,
    const float4* __restrict__ vel_T,
    const int* __restrict__ cell_D_start,   // CSR offsets for D particles
    const int* __restrict__ cell_T_start,   // CSR offsets for T particles
    ReactionProduct* __restrict__ products,  // output alpha+neutron births
    curandState* __restrict__ rng,
    float dt,
    GridParams grid,
    int N_cells)
{
    int cell_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell_id >= N_cells) return;

    int dStart = cell_D_start[cell_id],  dEnd = cell_D_start[cell_id + 1];
    int tStart = cell_T_start[cell_id],  tEnd = cell_T_start[cell_id + 1];
    int nD = dEnd - dStart, nT = tEnd - tStart;
    if (nD == 0 || nT == 0) return;

    curandState local_rng = rng[cell_id];

    float V_cell = grid.dx * grid.dy * grid.dz;

    // Loop over D particles in this cell; pair each with a random T
    for (int id = dStart; id < dEnd; id++) {
        float4 vD4 = vel_D[id];
        float4 pD4 = pos_D[id];

        // Pick a random T partner
        int it = tStart + (int)(curand_uniform(&local_rng) * nT);
        it = min(it, tEnd - 1);

        float4 vT4 = vel_T[it];
        float4 pT4 = pos_T[it];

        // Relative velocity → CoM kinetic energy
        float dvx = vD4.x - vT4.x;
        float dvy = vD4.y - vT4.y;
        float dvz = vD4.z - vT4.z;
        float vrel2 = dvx*dvx + dvy*dvy + dvz*dvz;

        float mu_DT = PC_MD * PC_MT / (PC_MD + PC_MT);
        float KE_J  = 0.5f * mu_DT * vrel2;
        float KE_keV = KE_J / (1.60217663e-19f * 1e3f);

        // <σv> from Bosch-Hale at this energy
        float sigma_v = boschHaleSigmaV(KE_keV);

        // Local T number density
        float n_T = (float)nT * pos_T[tStart].w / V_cell; // weighted

        // Reaction probability: P = n_T * <σv> * dt
        float P = n_T * sigma_v * dt;
        P = fminf(P, 1.0f); // cap at 1 (check timestep if frequently capped)

        if (curand_uniform(&local_rng) < P) {
            // Record the reaction
            float3 p3D = make_float3(pD4.x, pD4.y, pD4.z);
            float3 p3T = make_float3(pT4.x, pT4.y, pT4.z);
            dtReactionProducts(p3D, p3T, &local_rng, products[id]);
        } else {
            products[id].active = 0;
        }
    }

    rng[cell_id] = local_rng;
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchFusionReactions(
    const float4* pos_D, const float4* vel_D,
    const float4* pos_T, const float4* vel_T,
    const int* cell_D_start, const int* cell_T_start,
    ReactionProduct* products,
    curandState* rng,
    float dt,
    GridParams grid,
    int N_cells,
    cudaStream_t stream)
{
    constexpr int BLOCK = 128;
    int gridDim = (N_cells + BLOCK - 1) / BLOCK;
    sampleFusionReactions<<<gridDim, BLOCK, 0, stream>>>(
        pos_D, vel_D, pos_T, vel_T,
        cell_D_start, cell_T_start,
        products, rng, dt, grid, N_cells);
}
