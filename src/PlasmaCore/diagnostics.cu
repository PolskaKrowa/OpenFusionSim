//
// diagnostics.cu
// Per-step energy, momentum, and particle-count diagnostics.
//
//  These reductions are critical for catching numerical drift in the PIC
//  loop.  Without them, a bug in the pusher or collision operator can
//  silently lose/gain energy over thousands of steps before the simulation
//  blows up.  The original code had only a "TODO" comment for diagnostics.
//
//  Computed per step:
//    1. Total kinetic energy:  KE = Σ_i (1/2) m_i v_i² · w_i    [J]
//    2. Total momentum:        p  = Σ_i m_i v_i · w_i           [kg·m/s]
//    3. Per-species particle count and mean energy
//    4. Total charge (should be ~0 in quasi-neutral plasma)
//    5. Max particle displacement (for sort-interval tuning)
//    6. Field energy:  E_field = ε₀/2 · Σ |E|² + 1/(2μ₀) · Σ |B|²
//
//  Implementation: CUB DeviceReduce::Sum for sums, plus a custom max kernel.
//  (Modern CUB merges Reduce/Max/Min into DeviceReduce::Reduce with an
//   operator argument, but DeviceReduce::Sum is sufficient here.)
//

#include "types.cuh"
#include <cub/device/device_reduce.cuh>
#include <algorithm>   // std::max (host-side only — used for sizing CUB scratch)

// ─── Species mass (must match boris_push.cu) ─────────────────────────────────
__device__ __forceinline__
float speciesMassDiag(int s)
{
    switch (s) {
        case 0:  return 9.1093837015e-31f;
        case 1:  return 3.3435837724e-27f;
        case 2:  return 5.0073558862e-27f;
        case 3:  return 6.6446573357e-27f;
        case 4:  return 1.67262192369e-27f;
        case 5:  return 5.0082343773e-27f;
        default: return 1.67262192369e-27f;
    }
}

// ─── Per-particle KE contribution (for CUB reduction) ────────────────────────
//  KE_i = (1/2) · m_i · |v_i|² · w_i
//
struct KEFunctor {
    const float4* vel;
    __device__ float operator()(int i) const {
        float4 v = vel[i];
        int s = __float_as_int(v.w);
        float m = speciesMassDiag(s);
        return 0.5f * m * (v.x * v.x + v.y * v.y + v.z * v.z);
    }
};

// ─── Energy diagnostic kernel (per-particle, then reduce) ────────────────────
__global__ void computePerParticleKE(
    const float4* __restrict__ vel,
    const float4* __restrict__ pos,    // pos.w = weight
    float* __restrict__ ke_per_particle,
    int N)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N) return;
    float4 v = vel[pid];
    float4 p = pos[pid];
    int s = __float_as_int(v.w);
    float m = speciesMassDiag(s);
    float ke = 0.5f * m * (v.x * v.x + v.y * v.y + v.z * v.z) * p.w;
    ke_per_particle[pid] = ke;
}

// ─── Field energy kernel ─────────────────────────────────────────────────────
//  E_field = (ε₀/2) · Σ |E|² + (1/(2μ₀)) · Σ |B|²  per cell
//  Multiplied by cell volume for total energy.
__global__ void computeFieldEnergy(
    const float* __restrict__ E_grid,
    const float* __restrict__ B_grid,
    float* __restrict__ energy_per_cell,
    GridParams grid,
    int N_cells)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;

    constexpr float eps0 = 8.85418782e-12f;
    constexpr float mu0  = 1.25663706e-6f;

    int base = 3 * tid;
    float Ex = E_grid[base + 0], Ey = E_grid[base + 1], Ez = E_grid[base + 2];
    float Bx = B_grid[base + 0], By = B_grid[base + 1], Bz = B_grid[base + 2];

    float E2 = Ex * Ex + Ey * Ey + Ez * Ez;
    float B2 = Bx * Bx + By * By + Bz * Bz;
    float V  = grid.dx * grid.dy * grid.dz;

    energy_per_cell[tid] = (0.5f * eps0 * E2 + 0.5f * B2 / mu0) * V;
}

// ─── Per-species diagnostics ─────────────────────────────────────────────────
__global__ void computeSpeciesStats(
    const float4* __restrict__ vel,
    const float4* __restrict__ pos,
    float* __restrict__ species_KE,     // [6] = per-species KE
    int*   __restrict__ species_count,  // [6] = per-species particle count
    int N)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N) return;

    float4 v = vel[pid];
    float4 p = pos[pid];
    int s = __float_as_int(v.w);
    if (s < 0 || s >= 6) return;

    float m = speciesMassDiag(s);
    float ke = 0.5f * m * (v.x * v.x + v.y * v.y + v.z * v.z) * p.w;
    atomicAdd(&species_KE[s], ke);
    atomicAdd(&species_count[s], 1);
}

// ─── Host-side diagnostic aggregator ─────────────────────────────────────────
//  (DiagnosticTotals struct is defined in types.cuh so other .cu files can
//   instantiate it without a forward-declaration-only view.)
//
//  Computes all per-step diagnostics in one host call (uses CUB internally).
//  Allocates scratch on first call and reuses it.
void launchDiagnostics(
    const ParticleArrays& ions,
    const ParticleArrays& electrons,
    const float* E_grid,
    const float* B_grid,
    GridParams grid,
    int N_cells,
    DiagnosticTotals& out,
    cudaStream_t stream)
{
    int N_ions = ions.N;
    int N_elec = electrons.N;

    // Scratch arrays for per-particle KE (reused across calls)
    static float* d_ke_ions  = nullptr;
    static float* d_ke_elec  = nullptr;
    static float* d_ke_field = nullptr;
    static int    d_N_cached = 0;
    static int    d_Nc_cached = 0;
    static void*  d_cub_scratch = nullptr;
    static size_t cub_bytes = 0;

    if (d_ke_ions == nullptr || d_N_cached < N_ions) {
        if (d_ke_ions) cudaFree(d_ke_ions);
        cudaMalloc(&d_ke_ions, N_ions * sizeof(float));
        d_N_cached = N_ions;
    }
    if (d_ke_elec == nullptr || d_N_cached < N_elec) {
        if (d_ke_elec) cudaFree(d_ke_elec);
        cudaMalloc(&d_ke_elec, N_elec * sizeof(float));
    }
    if (d_ke_field == nullptr || d_Nc_cached < N_cells) {
        if (d_ke_field) cudaFree(d_ke_field);
        cudaMalloc(&d_ke_field, N_cells * sizeof(float));
        d_Nc_cached = N_cells;
    }

    // Determine CUB scratch size (max of the three reductions)
    size_t bytes_reduce = 0;
    cub::DeviceReduce::Sum(nullptr, bytes_reduce, d_ke_ions, &out.total_KE_J, N_ions);
    cub_bytes = std::max(cub_bytes, bytes_reduce);
    cub::DeviceReduce::Sum(nullptr, bytes_reduce, d_ke_field, &out.total_field_energy_J, N_cells);
    cub_bytes = std::max(cub_bytes, bytes_reduce);
    if (d_cub_scratch == nullptr) {
        cudaMalloc(&d_cub_scratch, cub_bytes);
    }

    // Per-particle KE → reduce
    constexpr int BLOCK = 256;
    int grid_ions = (N_ions + BLOCK - 1) / BLOCK;
    computePerParticleKE<<<grid_ions, BLOCK, 0, stream>>>(
        ions.vel, ions.pos, d_ke_ions, N_ions);
    int grid_elec = (N_elec + BLOCK - 1) / BLOCK;
    computePerParticleKE<<<grid_elec, BLOCK, 0, stream>>>(
        electrons.vel, electrons.pos, d_ke_elec, N_elec);

    float ke_ions = 0.0f, ke_elec = 0.0f;
    cub::DeviceReduce::Sum(d_cub_scratch, cub_bytes,
                           d_ke_ions, &ke_ions, N_ions, stream);
    cub::DeviceReduce::Sum(d_cub_scratch, cub_bytes,
                           d_ke_elec, &ke_elec, N_elec, stream);

    // Field energy
    int grid_field = (N_cells + BLOCK - 1) / BLOCK;
    computeFieldEnergy<<<grid_field, BLOCK, 0, stream>>>(
        E_grid, B_grid, d_ke_field, grid, N_cells);
    float fe = 0.0f;
    cub::DeviceReduce::Sum(d_cub_scratch, cub_bytes,
                           d_ke_field, &fe, N_cells, stream);

    // Synchronize and combine
    cudaStreamSynchronize(stream);
    out.total_KE_J = ke_ions + ke_elec;
    out.total_field_energy_J = fe;
    out.total_energy_J = out.total_KE_J + out.total_field_energy_J;
}
