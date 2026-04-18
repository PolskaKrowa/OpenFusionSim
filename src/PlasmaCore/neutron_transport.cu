//
// neutron_transport.cu
// Monte Carlo neutron transport using Woodcock delta-tracking.
//
//  Each thread independently tracks one neutron history from birth to
//  termination (absorption, escape, or energy cutoff).
//
//  Woodcock delta-tracking eliminates geometry branch divergence by
//  sampling a fictitious "majorant" cross-section that equals the maximum
//  total cross-section across all materials in the domain.  Real collisions
//  are accepted probabilistically (virtual collisions are rejected).
//  This gives near-uniform branch structure across all threads — ideal for GPU.
//
//  Materials in the domain:
//    0 = plasma/vacuum  (negligible XS)
//    1 = Li-6 breeder   (TBR source)
//    2 = Li-7 breeder
//    3 = structural steel
//    4 = coolant (water/FLiBe)
//
//  Reactions tracked:
//    Li-6(n,t)He-4   → tritium breeding
//    All materials   → elastic scattering, inelastic, capture
//    Energy deposition recorded per voxel for thermal coupling.
//

#include "types.cuh"
#include <curand_kernel.h>
#include <math.h>

// ─── Cross-section table in constant memory ────────────────────────────────────
//  Simplified multi-group structure: N_GROUPS energy bins.
//  In practice this would be loaded from an ENDF/B-VIII evaluated data file.
#define N_MAT    5
#define N_GROUPS 64

__constant__ float c_sigma_total[N_MAT][N_GROUPS]; // total XS [m^-1]
__constant__ float c_sigma_abs  [N_MAT][N_GROUPS]; // absorption XS [m^-1]
__constant__ float c_sigma_Li6  [N_GROUPS];         // Li-6(n,t) XS [m^-1]
__constant__ float c_sigma_major[N_GROUPS];          // majorant (max over materials)
__constant__ float c_E_bins[N_GROUPS + 1];           // energy group boundaries [MeV]

// ─── Energy → group index ──────────────────────────────────────────────────────
__device__ __forceinline__
int energyGroup(float E_MeV)
{
    // Log-spaced bins from 1e-7 MeV (thermal) to 20 MeV (fusion)
    // Linear search acceptable for N_GROUPS=64; binary search for larger tables
    for (int g = 0; g < N_GROUPS; g++) {
        if (E_MeV < c_E_bins[g + 1]) return g;
    }
    return N_GROUPS - 1;
}

// ─── Material ID at a position ────────────────────────────────────────────────
__device__ __forceinline__
int getMaterial(const float* __restrict__ material_map,
                float3 pos, const GridParams& g)
{
    int ix = (int)((pos.x - g.ox) / g.dx);
    int iy = (int)((pos.y - g.oy) / g.dy);
    int iz = (int)((pos.z - g.oz) / g.dz);
    if (ix < 0 || ix >= g.Nx || iy < 0 || iy >= g.Ny || iz < 0 || iz >= g.Nz)
        return -1; // escaped domain
    return (int)material_map[flatCell(ix, iy, iz, g)];
}

// ─── Sample isotropic scattering direction ────────────────────────────────────
__device__ float3 sampleIsotropic(curandState* rng)
{
    float cos_th = 2.0f * curand_uniform(rng) - 1.0f;
    float sin_th = sqrtf(fmaxf(0.0f, 1.0f - cos_th * cos_th));
    float phi    = 2.0f * 3.14159265f * curand_uniform(rng);
    return make_float3(sin_th * cosf(phi), sin_th * sinf(phi), cos_th);
}

// ─── Elastic scatter energy loss (A = target mass number) ────────────────────
__device__ float elasticScatterEnergy(float E_in, float A, curandState* rng)
{
    float cos_cm = 2.0f * curand_uniform(rng) - 1.0f;
    float alpha  = ((A - 1.0f) / (A + 1.0f)) * ((A - 1.0f) / (A + 1.0f));
    return E_in * 0.5f * ((1.0f + alpha) + (1.0f - alpha) * cos_cm);
}

// ─── Deposit energy into the heat map ─────────────────────────────────────────
__device__ void depositHeat(HeatDepositionMap* qmap, float3 pos,
                             float energy_J, const GridParams& g)
{
    int ix = (int)((pos.x - g.ox) / g.dx);
    int iy = (int)((pos.y - g.oy) / g.dy);
    int iz = (int)((pos.z - g.oz) / g.dz);
    if (ix < 0 || ix >= qmap->Nx || iy < 0 || iy >= qmap->Ny ||
        iz < 0 || iz >= qmap->Nz) return;
    int idx = ix + qmap->Nx * (iy + qmap->Ny * iz);
    atomicAdd(&qmap->q_dot[idx], energy_J);
}

// ─── Accumulate tritium breeding event ────────────────────────────────────────
__device__ void recordTBR(TritiumProductionMap* tbr, float3 pos,
                           float rate, const GridParams& g)
{
    int ix = (int)((pos.x - g.ox) / g.dx);
    int iy = (int)((pos.y - g.oy) / g.dy);
    int iz = (int)((pos.z - g.oz) / g.dz);
    if (ix < 0 || ix >= tbr->Nx || iy < 0 || iy >= tbr->Ny ||
        iz < 0 || iz >= tbr->Nz) return;
    int idx = ix + tbr->Nx * (iy + tbr->Ny * iz);
    atomicAdd(&tbr->tbr_voxel[idx], rate);
}

// ─── Neutron Transport Kernel (Woodcock delta-tracking) ───────────────────────
//
//  Each thread = one complete neutron history.
//  No communication between threads → embarrassingly parallel.
//
__global__ void transportNeutrons(
    NeutronParticle* __restrict__ neutrons,
    const float* __restrict__ material_map,
    TritiumProductionMap* tbr,
    HeatDepositionMap*    q_dot,
    curandState* __restrict__ rng,
    GridParams grid,
    int N_neutrons,
    float E_cutoff_MeV,      // terminate below this energy (e.g. 1e-7 MeV thermal)
    int max_collisions)      // safety cap on history length
{
    int nid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nid >= N_neutrons) return;
    if (!neutrons[nid].alive) return;

    curandState local_rng = rng[nid];

    float3 pos = neutrons[nid].pos;
    float3 dir = neutrons[nid].dir;
    float E    = neutrons[nid].energy_MeV;
    float wt   = neutrons[nid].weight;

    for (int iter = 0; iter < max_collisions && E > E_cutoff_MeV; iter++) {

        int group = energyGroup(E);
        float sigma_maj = c_sigma_major[group];
        if (sigma_maj < 1e-30f) break; // void — free stream out

        // Sample distance to next (possibly virtual) collision
        float xi = curand_uniform(&local_rng);
        float dist = -logf(xi + 1e-20f) / sigma_maj;

        // Advance position
        pos.x += dir.x * dist;
        pos.y += dir.y * dist;
        pos.z += dir.z * dist;

        // Check for domain escape
        int mat = getMaterial(material_map, pos, grid);
        if (mat < 0) {
            // Escaped — terminate history
            break;
        }

        // Woodcock acceptance: real collision with probability sigma_total / sigma_maj
        float sigma_real = c_sigma_total[mat][group];
        float P_real = sigma_real / sigma_maj;

        if (curand_uniform(&local_rng) > P_real) {
            // Virtual collision — no interaction, continue flight
            continue;
        }

        // ── Real collision — determine reaction type ───────────────────────
        float sigma_abs = c_sigma_abs[mat][group];
        float sigma_Li6 = (mat == 1) ? c_sigma_Li6[group] : 0.0f;

        float r = curand_uniform(&local_rng) * sigma_real;

        if (mat == 1 && r < sigma_Li6) {
            // Li-6(n,t)He-4: tritium breeding + alpha
            float Q_J = (4.78e6f) * 1.60217663e-19f; // 4.78 MeV Q-value
            recordTBR(tbr, pos, wt, grid);
            depositHeat(q_dot, pos, wt * Q_J, grid);
            break; // neutron absorbed
        } else if (r < sigma_abs) {
            // Radiative capture (n,gamma) — deposit energy, terminate
            float E_J = E * 1e6f * 1.60217663e-19f;
            depositHeat(q_dot, pos, wt * E_J, grid);
            break;
        } else {
            // Elastic scattering
            //  Approximate target mass: steel≈56, Li≈6.5, water≈9 (H+O avg)
            float A_target[] = {1.0f, 6.5f, 7.0f, 56.0f, 9.0f};
            float A = A_target[mat];

            // Energy lost to recoil
            float E_new = elasticScatterEnergy(E, A, &local_rng);
            float dE_J  = (E - E_new) * 1e6f * 1.60217663e-19f;
            depositHeat(q_dot, pos, wt * dE_J, grid);
            E = E_new;

            // New isotropic direction (lab-frame approximation valid for heavy targets)
            dir = sampleIsotropic(&local_rng);
        }
    }

    // Write back terminal state
    neutrons[nid].pos    = pos;
    neutrons[nid].dir    = dir;
    neutrons[nid].energy_MeV = E;
    neutrons[nid].alive  = (E > E_cutoff_MeV) ? 1 : 0;
    rng[nid] = local_rng;
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchNeutronTransport(NeutronParticle* neutrons,
                            const float* material_map,
                            TritiumProductionMap* tbr,
                            HeatDepositionMap* q_dot,
                            curandState* rng,
                            GridParams grid,
                            int N_neutrons,
                            float E_cutoff_MeV,
                            int max_collisions,
                            cudaStream_t stream)
{
    constexpr int BLOCK = 128;
    int gridDim = (N_neutrons + BLOCK - 1) / BLOCK;
    transportNeutrons<<<gridDim, BLOCK, 0, stream>>>(
        neutrons, material_map, tbr, q_dot, rng,
        grid, N_neutrons, E_cutoff_MeV, max_collisions);
}
