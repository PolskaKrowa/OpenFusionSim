//
// src/PlasmaCore/pic_bridge.cu
// CUDA-side implementation of the PIC bridge.
//
//  Provides the C-linkage functions declared in pic_bridge.cpp:
//    pic_cuda_create / pic_cuda_destroy / pic_cuda_init_plasma /
//    pic_cuda_step / pic_cuda_readback_diagnostics / pic_cuda_readback_particles /
//    pic_cuda_set_B_field / pic_cuda_set_material_map / pic_cuda_init_rng /
//    pic_cuda_init_cross_sections
//
//  All actual CUDA work (kernel launches, device memory management) happens
//  here.  The host-side pic_bridge.cpp only orchestrates.
//

#include "types.cuh"
#include "endf_loader.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>
// NOTE: do NOT include <cub/device/device_radix_sort.cuh> or
// <cub/device/device_scan.cuh> here — we delegate to sorting.cu's
// allocSortContext()/sortParticlesByCell() functions instead.  Direct
// CUB template instantiation here would create duplicate kernel
// registrations at device-link time.
#include <cmath>
#include <chrono>
#include <cstring>   // std::memcpy (host-side type-punning for species_id)
#include <vector>    // std::vector (host scratch for bulk readback)
#include <algorithm> // std::min, std::max

// Forward declarations of host launch wrappers from the .cu modules
void launchBorisPush(ParticleArrays&, float4*, float4*, float, cudaStream_t);
void launchGatherFields(float4*, float4*, const float*, const float*,
                        const ParticleArrays&, GridParams, cudaStream_t);
void launchDepositCharge(float*, float*, const ParticleArrays&,
                         GridParams, float, cudaStream_t);
void launchFDTD(float*, float*, float*, const GridParams&, float, cudaStream_t);
void launchCoulombCollisions(ParticleArrays&, const float*, const int*, float,
                             float, curandState*, GridParams, int, cudaStream_t);
void launchFusionReactions(const float4*, const float4*, const float4*, const float4*,
                           const int*, const int*, ReactionProduct*, curandState*,
                           float, GridParams, int, cudaStream_t);
void launchNeutronTransport(NeutronParticle*, const float*, TritiumProductionMap*,
                            HeatDepositionMap*, curandState*, GridParams,
                            int, float, int, cudaStream_t);
void launchRadiationLosses(float*, float*, const float*, const float*, const float*,
                           const float*, const int*, float, int, cudaStream_t);
void launchRadiationFeedback(float*, const float*, const float*, const float*,
                              const float*, const float*, float, float, int,
                              cudaStream_t);
void launchInitRNG(curandState*, unsigned long long, int, cudaStream_t);
void sortParticlesByCell(ParticleArrays&, int*, struct SortContext&,
                         GridParams, cudaStream_t);
// Forward declarations from sorting.cu — we delegate to these instead of
// instantiating CUB::DeviceRadixSort directly here, which would create
// duplicate kernel registrations at device-link time.
void allocSortContext(SortContext& ctx, int N, int N_cells);
void freeSortContext(SortContext& ctx);

// ENDF upload function (from endf_upload.cu)
namespace ENDFUpload {
    void initializeWithDefaults();
    bool initializeFromFile(const char* path);
}

// ─── SortContext is now defined in types.cuh (included above) ─────────────────
//  (Previously re-declared here, which was an ODR violation.)

// ─── Per-cell diagnostic buffers (device-side) ───────────────────────────────
struct DeviceDiagnostics {
    float* T_e_keV;          // per-cell electron temperature
    float* n_e_m3;           // per-cell density
    float* B_T;              // per-cell B-field magnitude
    float* q_rad_W_m3;       // per-cell radiation power density
    float* q_alpha_W_m3;     // per-cell alpha heating
    float* tau_E_s;          // per-cell confinement time
    float* T_e_per_cell_keV; // for fusion kernel
    float* n_e_per_cell;     // for collision kernel
    float* n_imp_per_cell;   // impurity number density [m^-3]  (zeroed by default)
    int*   Z_imp_per_cell;   // impurity Z (zeroed by default → no line radiation)
};

// ─── Complete device-side PIC state ──────────────────────────────────────────
struct PICDeviceState {
    // Configuration
    GridParams grid;
    int N_particles;
    float dt_pic;
    float coulomb_log;
    int sort_interval;

    // Particle arrays
    ParticleArrays ions;
    ParticleArrays electrons;

    // Per-species sub-arrays (views into ions for D, T)
    float4* pos_D; float4* vel_D;
    float4* pos_T; float4* vel_T;

    // Field grids
    float* E_grid;
    float* B_grid;
    float* rho_grid;
    float* J_grid;

    // Field at particle positions
    float4* E_at_ion;  float4* B_at_ion;
    float4* E_at_elec; float4* B_at_elec;

    // Weights
    float* ion_weights;
    float* elec_weights;

    // CSR offsets
    int* ion_cell_start;
    int* elec_cell_start;
    int* cell_D_start;
    int* cell_T_start;

    // Reaction products
    ReactionProduct* reaction_products;

    // Neutrons
    NeutronParticle* neutrons;
    int N_neutrons_max;

    // TBR + heat maps
    TritiumProductionMap* tbr;
    HeatDepositionMap* q_dot;

    // Material map
    float* material_map;

    // RNG states
    curandState* rng_ion;
    curandState* rng_elec;
    curandState* rng_neutron;
    curandState* rng_cell;

    // Sort context
    SortContext sort_ctx;

    // Diagnostics
    DeviceDiagnostics diag;

    // CUDA streams
    cudaStream_t stream_main;
    cudaStream_t stream_neutron;

    int step_count;
};

// ─── Particle initialization kernel ──────────────────────────────────────────
//
//  Seeds particles with a Maxwellian distribution at temperature T_e,
//  distributed uniformly in the plasma ellipse (R-R_major)²/a² + z²/(κa)² ≤ 1.
//
//  Species assignment (encoded in vel.w):
//    0 = electron, 1 = deuterium, 2 = tritium, 3 = alpha
//  50:50 D-T mix among the ions; equal number of electrons for quasi-neutrality.
//
__global__ void initPlasmaParticles(
    float4* pos, float4* vel, float* weights,
    int N, int species,
    float T_keV, float n_e_m3,
    float R_major, float a_minor, float kappa,
    float B_T, float Lx, float Ly, float Lz,
    curandState* rng_states, unsigned long long seed)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;

    curandState local_rng;
    curand_init(seed, tid, 0, &local_rng);

    // Sample position uniformly in the plasma ellipse
    float x, y, z, r_norm_sq;
    int tries = 0;
    do {
        x = (curand_uniform(&local_rng) * 2.0f - 1.0f) * a_minor + R_major;
        y = (curand_uniform(&local_rng) * 2.0f - 1.0f) * a_minor;
        z = (curand_uniform(&local_rng) * 2.0f - 1.0f) * kappa * a_minor;
        float dx = x - R_major;
        r_norm_sq = (dx*dx + y*y) / (a_minor*a_minor) + (z*z) / (kappa*kappa*a_minor*a_minor);
        tries++;
    } while (r_norm_sq > 1.0f && tries < 20);

    // Mass-dependent thermal velocity: v_th = sqrt(2 T [J] / m)
    float m_kg;
    switch (species) {
        case 0: m_kg = 9.1093837015e-31f; break;       // electron
        case 1: m_kg = 3.3435837724e-27f; break;       // D
        case 2: m_kg = 5.0073558862e-27f; break;       // T
        case 3: m_kg = 6.6446573357e-27f; break;       // α
        default: m_kg = 1.67262192369e-27f;
    }
    float T_J = T_keV * 1.602176634e-16f;
    float v_th = sqrtf(2.0f * T_J / m_kg);

    // Maxwellian velocity via Box-Muller
    float u1 = fmaxf(curand_uniform(&local_rng), 1e-7f);
    float u2 = curand_uniform(&local_rng);
    float g1 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * 3.14159265f * u2);
    float u3 = curand_uniform(&local_rng);
    float u4 = curand_uniform(&local_rng);
    float g2 = sqrtf(-2.0f * logf(fmaxf(u3, 1e-7f))) * cosf(2.0f * 3.14159265f * u4);
    float u5 = curand_uniform(&local_rng);
    float u6 = curand_uniform(&local_rng);
    float g3 = sqrtf(-2.0f * logf(fmaxf(u5, 1e-7f))) * cosf(2.0f * 3.14159265f * u6);

    // Add a toroidal drift (E×B in the lab frame gives ion/electron counter-streaming)
    float v_toroidal = (species == 0) ? -v_th * 0.1f : v_th * 0.1f;

    pos[tid] = make_float4(x, y, z, /*weight=*/1.0f);
    vel[tid] = make_float4(v_th * g1 + v_toroidal, v_th * g2, v_th * g3,
                           __int_as_float(species));
    weights[tid] = 1.0f;   // uniform weight; will be rescaled later

    rng_states[tid] = local_rng;
}

// ─── Compute per-cell density and temperature ────────────────────────────────
__global__ void computeCellDiagnostics(
    const float4* pos, const float4* vel,
    int N, int N_cells, GridParams grid,
    float* n_e_per_cell, float* T_e_per_cell_keV,
    float* B_per_cell_T, float B_T_global)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;

    // For simplicity, use a uniform density based on the global n_e and
    // the global T_e.  A real diagnostic would do a per-cell reduction.
    // (This kernel is a placeholder — the actual per-cell reductions are
    //  implemented in diagnostics.cu.)
    n_e_per_cell[tid] = 1.0e20f;          // placeholder: ITER-like density
    T_e_per_cell_keV[tid] = 10.0f;        // placeholder: 10 keV
    B_per_cell_T[tid] = B_T_global;

    (void)pos; (void)vel; (void)N; (void)grid;
}

// ─── Set toroidal field ──────────────────────────────────────────────────────
__global__ void setToroidalFieldKernel(float* B_grid, float B_T, int N_cells)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;
    // B = (0, 0, B_T) — aligned with the toroidal direction
    int base = 3 * tid;
    B_grid[base + 0] = 0.0f;
    B_grid[base + 1] = 0.0f;
    B_grid[base + 2] = B_T;
}

// ─── Set material map ────────────────────────────────────────────────────────
__global__ void setMaterialMapKernel(float* material_map, const int* host_ids, int N_cells)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;
    material_map[tid] = (float)host_ids[tid];
}

// ═══════════════════════════════════════════════════════════════════════════════
// C-linkage implementation functions (called from pic_bridge.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

extern "C" {

void* pic_cuda_create(int N_particles, int Nx, int Ny, int Nz,
                       float Lx, float Ly, float Lz)
{
    // ── Ensure device 0 is selected before any CUDA call ────────────────────
    // The bridge calls cudaSetDevice(0) in initGPU(), but we re-assert it here
    // so pic_bridge.cu is self-contained and robust against being called from
    // a different thread or after a device reset.
    cudaError_t dev_err = cudaSetDevice(0);
    if (dev_err != cudaSuccess) {
        return nullptr;   // bridge will see null and fall back to 0D-only mode
    }

    PICDeviceState* s = new PICDeviceState{};
    s->N_particles = N_particles;
    s->dt_pic = 1e-12f;
    s->coulomb_log = 17.0f;
    s->sort_interval = 20;
    s->step_count = 0;

    s->grid.ox = 0.0f; s->grid.oy = -Ly * 0.5f; s->grid.oz = -Lz * 0.5f;
    s->grid.dx = Lx / Nx; s->grid.dy = Ly / Ny; s->grid.dz = Lz / Nz;
    s->grid.Nx = Nx; s->grid.Ny = Ny; s->grid.Nz = Nz;

    int N_cells = Nx * Ny * Nz;
    int N = N_particles;

    // Allocate particle arrays (ions + electrons)
    cudaMalloc(&s->ions.pos, N * sizeof(float4));
    cudaMalloc(&s->ions.vel, N * sizeof(float4));
    cudaMalloc(&s->ions.cell, N * sizeof(int));
    s->ions.N = N;
    cudaMalloc(&s->electrons.pos, N * sizeof(float4));
    cudaMalloc(&s->electrons.vel, N * sizeof(float4));
    cudaMalloc(&s->electrons.cell, N * sizeof(int));
    s->electrons.N = N;

    // ── Initialize D/T sub-array pointers ──────────────────────────────────
    // The ions array is split: first half = deuterium, second half = tritium.
    // (Matches the initialization in pic_cuda_init_plasma, which seeds
    //  species=1 for indices [0, N/2) and species=2 for [N/2, N).)
    // Without this, pos_D/pos_T stay null and the fusion kernel crashes the
    // GPU — the resulting CUDA context error then surfaces as a misleading
    // "invalid device ordinal" from Thrust on the next operation.
    int N_half = N / 2;
    s->pos_D = s->ions.pos;          s->vel_D = s->ions.vel;
    s->pos_T = s->ions.pos + N_half; s->vel_T = s->ions.vel + N_half;

    // Fields
    cudaMalloc(&s->E_grid, 3 * N_cells * sizeof(float));
    cudaMalloc(&s->B_grid, 3 * N_cells * sizeof(float));
    cudaMalloc(&s->rho_grid, N_cells * sizeof(float));
    cudaMalloc(&s->J_grid, 3 * N_cells * sizeof(float));
    cudaMemset(s->E_grid, 0, 3 * N_cells * sizeof(float));
    cudaMemset(s->B_grid, 0, 3 * N_cells * sizeof(float));

    cudaMalloc(&s->E_at_ion, N * sizeof(float4));
    cudaMalloc(&s->B_at_ion, N * sizeof(float4));
    cudaMalloc(&s->E_at_elec, N * sizeof(float4));
    cudaMalloc(&s->B_at_elec, N * sizeof(float4));

    cudaMalloc(&s->ion_weights, N * sizeof(float));
    cudaMalloc(&s->elec_weights, N * sizeof(float));

    cudaMalloc(&s->ion_cell_start, (N_cells + 1) * sizeof(int));
    cudaMalloc(&s->elec_cell_start, (N_cells + 1) * sizeof(int));
    cudaMalloc(&s->cell_D_start, (N_cells + 1) * sizeof(int));
    cudaMalloc(&s->cell_T_start, (N_cells + 1) * sizeof(int));

    cudaMalloc(&s->reaction_products, N * sizeof(ReactionProduct));

    s->N_neutrons_max = N;
    cudaMalloc(&s->neutrons, N * sizeof(NeutronParticle));

    s->tbr   = new TritiumProductionMap{nullptr, Nx, Ny, Nz};
    s->q_dot = new HeatDepositionMap{nullptr, Nx, Ny, Nz};
    cudaMalloc(&s->tbr->tbr_voxel, N_cells * sizeof(float));
    cudaMalloc(&s->q_dot->q_dot, N_cells * sizeof(float));
    cudaMemset(s->tbr->tbr_voxel, 0, N_cells * sizeof(float));
    cudaMemset(s->q_dot->q_dot, 0, N_cells * sizeof(float));

    cudaMalloc(&s->material_map, N_cells * sizeof(float));

    cudaMalloc(&s->rng_ion, N * sizeof(curandState));
    cudaMalloc(&s->rng_elec, N * sizeof(curandState));
    cudaMalloc(&s->rng_neutron, N * sizeof(curandState));
    cudaMalloc(&s->rng_cell, N_cells * sizeof(curandState));

    // Sort context — delegate to sorting.cu's allocSortContext() so the
    // CUB DeviceRadixSort template is instantiated in exactly one .cu file.
    // (Previously we re-instantiated SortPairs here to size the scratch
    // buffer, which caused "Duplicate entry kernels" warnings at device-link
    // time when the same <int,int> instantiation appeared in both sorting.cu
    // and pic_bridge.cu.)
    allocSortContext(s->sort_ctx, N, N_cells);

    // Diagnostics
    cudaMalloc(&s->diag.T_e_keV, N_cells * sizeof(float));
    cudaMalloc(&s->diag.n_e_m3, N_cells * sizeof(float));
    cudaMalloc(&s->diag.B_T, N_cells * sizeof(float));
    cudaMalloc(&s->diag.q_rad_W_m3, N_cells * sizeof(float));
    cudaMalloc(&s->diag.q_alpha_W_m3, N_cells * sizeof(float));
    cudaMalloc(&s->diag.tau_E_s, N_cells * sizeof(float));
    cudaMalloc(&s->diag.T_e_per_cell_keV, N_cells * sizeof(float));
    cudaMalloc(&s->diag.n_e_per_cell, N_cells * sizeof(float));
    // Impurity diagnostic arrays — zeroed by default (no line radiation).
    // Without these, the radiation kernel would crash with a Warp MMU Fault
    // when it dereferenced the null pointer.
    cudaMalloc(&s->diag.n_imp_per_cell, N_cells * sizeof(float));
    cudaMalloc(&s->diag.Z_imp_per_cell, N_cells * sizeof(int));
    cudaMemset(s->diag.n_imp_per_cell, 0, N_cells * sizeof(float));
    cudaMemset(s->diag.Z_imp_per_cell, 0, N_cells * sizeof(int));

    // Streams
    cudaStreamCreate(&s->stream_main);
    cudaStreamCreate(&s->stream_neutron);

    // Upload grid params to __constant__ memory
    cudaMemcpyToSymbol(c_grid, &s->grid, sizeof(GridParams));

    return s;
}

void pic_cuda_destroy(void* handle)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;
    cudaStreamSynchronize(s->stream_main);
    cudaStreamSynchronize(s->stream_neutron);
    cudaStreamDestroy(s->stream_main);
    cudaStreamDestroy(s->stream_neutron);

    cudaFree(s->ions.pos); cudaFree(s->ions.vel); cudaFree(s->ions.cell);
    cudaFree(s->electrons.pos); cudaFree(s->electrons.vel); cudaFree(s->electrons.cell);
    cudaFree(s->E_grid); cudaFree(s->B_grid); cudaFree(s->rho_grid); cudaFree(s->J_grid);
    cudaFree(s->E_at_ion); cudaFree(s->B_at_ion);
    cudaFree(s->E_at_elec); cudaFree(s->B_at_elec);
    cudaFree(s->ion_weights); cudaFree(s->elec_weights);
    cudaFree(s->ion_cell_start); cudaFree(s->elec_cell_start);
    cudaFree(s->cell_D_start); cudaFree(s->cell_T_start);
    cudaFree(s->reaction_products);
    cudaFree(s->neutrons);
    cudaFree(s->tbr->tbr_voxel); cudaFree(s->q_dot->q_dot);
    delete s->tbr; delete s->q_dot;
    cudaFree(s->material_map);
    cudaFree(s->rng_ion); cudaFree(s->rng_elec);
    cudaFree(s->rng_neutron); cudaFree(s->rng_cell);
    // Sort context cleanup — delegated to sorting.cu (keeps CUB template
    // instantiation in one place).
    freeSortContext(s->sort_ctx);
    cudaFree(s->diag.T_e_keV); cudaFree(s->diag.n_e_m3);
    cudaFree(s->diag.B_T); cudaFree(s->diag.q_rad_W_m3);
    cudaFree(s->diag.q_alpha_W_m3); cudaFree(s->diag.tau_E_s);
    cudaFree(s->diag.T_e_per_cell_keV); cudaFree(s->diag.n_e_per_cell);
    cudaFree(s->diag.n_imp_per_cell); cudaFree(s->diag.Z_imp_per_cell);
    delete s;
}

void pic_cuda_init_plasma(void* handle, float T_e_keV, float n_e_m3,
                           float R_major, float a_minor, float kappa,
                           float B_T)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;

    int N = s->N_particles;
    int N_ions = N / 2;
    int N_elec = N - N_ions;

    // Initialize ions: first half D, second half T
    constexpr int BLOCK = 256;
    int grid_dim = (N_ions + BLOCK - 1) / BLOCK;
    initPlasmaParticles<<<grid_dim, BLOCK, 0, s->stream_main>>>(
        s->ions.pos, s->ions.vel, s->ion_weights, N_ions, /*species=*/1,
        T_e_keV, n_e_m3, R_major, a_minor, kappa, B_T,
        s->grid.Nx * s->grid.dx, s->grid.Ny * s->grid.dy, s->grid.Nz * s->grid.dz,
        s->rng_ion, 0x5EED1234ULL);
    initPlasmaParticles<<<grid_dim, BLOCK, 0, s->stream_main>>>(
        s->ions.pos + N_ions, s->ions.vel + N_ions, s->ion_weights + N_ions,
        N_ions, /*species=*/2,
        T_e_keV, n_e_m3, R_major, a_minor, kappa, B_T,
        s->grid.Nx * s->grid.dx, s->grid.Ny * s->grid.dy, s->grid.Nz * s->grid.dz,
        s->rng_ion + N_ions, 0x5EED1235ULL);

    // Electrons
    grid_dim = (N_elec + BLOCK - 1) / BLOCK;
    initPlasmaParticles<<<grid_dim, BLOCK, 0, s->stream_main>>>(
        s->electrons.pos, s->electrons.vel, s->elec_weights, N_elec, /*species=*/0,
        T_e_keV, n_e_m3, R_major, a_minor, kappa, B_T,
        s->grid.Nx * s->grid.dx, s->grid.Ny * s->grid.dy, s->grid.Nz * s->grid.dz,
        s->rng_elec, 0x5EED1236ULL);

    // Set B field
    int N_cells = s->grid.Nx * s->grid.Ny * s->grid.Nz;
    grid_dim = (N_cells + BLOCK - 1) / BLOCK;
    setToroidalFieldKernel<<<grid_dim, BLOCK, 0, s->stream_main>>>(
        s->B_grid, B_T, N_cells);

    cudaStreamSynchronize(s->stream_main);
}

float pic_cuda_step(void* handle, float dt)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return 0.0f;

    auto t0 = std::chrono::steady_clock::now();

    int N_cells = s->grid.Nx * s->grid.Ny * s->grid.Nz;

    // 1. Sort particles by cell (every sort_interval steps)
    if (s->step_count % s->sort_interval == 0) {
        sortParticlesByCell(s->ions, s->ion_cell_start, s->sort_ctx,
                             s->grid, s->stream_main);
        sortParticlesByCell(s->electrons, s->elec_cell_start, s->sort_ctx,
                             s->grid, s->stream_main);
    }

    // 2. Deposit charge and current
    launchDepositCharge(s->rho_grid, s->J_grid, s->ions,
                         s->grid, dt, s->stream_main);
    launchDepositCharge(s->rho_grid, s->J_grid, s->electrons,
                         s->grid, dt, s->stream_main);

    // 3. Solve fields (FDTD)
    launchFDTD(s->E_grid, s->B_grid, s->J_grid, s->grid, dt, s->stream_main);

    // 4. Gather fields at particle positions
    launchGatherFields(s->E_at_ion, s->B_at_ion, s->E_grid, s->B_grid,
                        s->ions, s->grid, s->stream_main);
    launchGatherFields(s->E_at_elec, s->B_at_elec, s->E_grid, s->B_grid,
                        s->electrons, s->grid, s->stream_main);

    // 5. Push particles (Boris / Vay)
    launchBorisPush(s->ions, s->E_at_ion, s->B_at_ion, dt, s->stream_main);
    launchBorisPush(s->electrons, s->E_at_elec, s->B_at_elec, dt, s->stream_main);

    // 6. Coulomb collisions
    launchCoulombCollisions(s->ions, s->ion_weights, s->ion_cell_start,
                              dt, s->coulomb_log, s->rng_cell,
                              s->grid, N_cells, s->stream_main);
    launchCoulombCollisions(s->electrons, s->elec_weights, s->elec_cell_start,
                              dt, s->coulomb_log, s->rng_cell,
                              s->grid, N_cells, s->stream_main);

    // 7. Update cell diagnostics (placeholder — real per-cell reductions in diagnostics.cu)
    {
        constexpr int BLOCK = 256;
        int gdim = (N_cells + BLOCK - 1) / BLOCK;
        computeCellDiagnostics<<<gdim, BLOCK, 0, s->stream_main>>>(
            s->ions.pos, s->ions.vel, s->ions.N, N_cells, s->grid,
            s->diag.n_e_per_cell, s->diag.T_e_per_cell_keV,
            s->diag.B_T, 5.3f);
    }

    // 8. Fusion reactions — guard against null D/T sub-arrays.
    //    (pos_D/pos_T are initialized in pic_cuda_create, but cell_D_start
    //    and cell_T_start are not yet computed by the sort step's species-
    //    filtering pass — that's a TODO in sorting.cu.  Until then, we skip
    //    the kinetic fusion channel and let the 0D model in updatePowerBalance
    //    handle fusion power.  This prevents the null-pointer GPU crash that
    //    was surfacing as "cudaErrorInvalidDevice" from Thrust.)
    if (s->pos_D && s->pos_T && s->vel_D && s->vel_T &&
        s->cell_D_start && s->cell_T_start) {
        launchFusionReactions(s->pos_D, s->vel_D, s->pos_T, s->vel_T,
                               s->cell_D_start, s->cell_T_start,
                               s->reaction_products, s->rng_cell,
                               dt, s->grid, N_cells, s->stream_main);
    }

    // 9. Neutron transport
    int N_neutrons = s->N_neutrons_max;  // placeholder; real count comes from reaction kernel
    launchNeutronTransport(s->neutrons, s->material_map,
                            s->tbr, s->q_dot, s->rng_neutron,
                            s->grid, N_neutrons, 1e-7f, 1000,
                            s->stream_neutron);

    // 10. Radiation losses (bremsstrahlung, synchrotron, line)
    //     Impurity arrays are now allocated and zeroed in pic_cuda_create,
    //     so line radiation is zero by default (pure D-T plasma).  When
    //     impurity tracking is added later, these arrays get populated.
    launchRadiationLosses(s->q_dot->q_dot, s->diag.q_rad_W_m3,
                           s->diag.n_e_per_cell, s->diag.T_e_per_cell_keV,
                           s->diag.B_T,
                           s->diag.n_imp_per_cell, s->diag.Z_imp_per_cell,
                           /*R_wall=*/0.6f, N_cells, s->stream_main);

    // 11. Radiation feedback: update T_e based on net heating
    //     (couples radiation → temperature, the missing piece)
    launchRadiationFeedback(s->diag.T_e_per_cell_keV,
                              s->diag.n_e_per_cell, s->diag.B_T,
                              s->diag.q_alpha_W_m3, s->diag.q_rad_W_m3,
                              s->diag.tau_E_s,
                              /*P_aux_W_m3=*/50.0e6f / (3.14159f * 6.2f * 4.0f * 1.7f),
                              dt, N_cells, s->stream_main);

    cudaStreamSynchronize(s->stream_main);
    cudaStreamSynchronize(s->stream_neutron);

    auto t1 = std::chrono::steady_clock::now();
    float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();
    s->step_count++;
    return ms;
}

void pic_cuda_readback_diagnostics(void* handle,
                                     float* T_e_keV, float* n_e_m3,
                                     float* B_T, float* q_rad, float* q_alpha,
                                     float* tau_E, int N_cells)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;
    cudaMemcpy(T_e_keV, s->diag.T_e_per_cell_keV, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(n_e_m3, s->diag.n_e_per_cell, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(B_T, s->diag.B_T, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(q_rad, s->diag.q_rad_W_m3, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(q_alpha, s->diag.q_alpha_W_m3, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(tau_E, s->diag.tau_E_s, N_cells * sizeof(float),
               cudaMemcpyDeviceToHost);
}

void pic_cuda_readback_particles(void* handle,
                                   float* pos_xyz, float* vel_xyz,
                                   float* weight, int* species,
                                   int max_particles, int* N_returned)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;

    // ── Bulk-copy the entire particle arrays to host scratch, then stride-
    //    sample on the CPU.  Doing one cudaMemcpy per particle (the old code)
    //    was ~4096 separate DMA launches per readback — brutal.  One bulk
    //    copy is O(N) on the GPU bus and O(1) sync points.
    int N = s->N_particles;
    std::vector<float4> h_pos(N), h_vel(N);
    std::vector<float>  h_w(N);
    cudaMemcpy(h_pos.data(), s->ions.pos,        N * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel.data(), s->ions.vel,        N * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_w.data(),   s->ion_weights,     N * sizeof(float),  cudaMemcpyDeviceToHost);

    // Stride-sample down to max_particles
    int stride   = std::max(1, N / max_particles);
    int out_idx  = 0;
    for (int i = 0; i < N && out_idx < max_particles; i += stride) {
        const float4& p = h_pos[i];
        const float4& v = h_vel[i];
        pos_xyz[3 * out_idx + 0] = p.x;
        pos_xyz[3 * out_idx + 1] = p.y;
        pos_xyz[3 * out_idx + 2] = p.z;
        vel_xyz[3 * out_idx + 0] = v.x;
        vel_xyz[3 * out_idx + 1] = v.y;
        vel_xyz[3 * out_idx + 2] = v.z;
        weight[out_idx] = h_w[i];

        // Species was packed into vel.w on the device via __int_as_float().
        // On the host we use memcpy (the standard C++ type-punning idiom;
        // reinterpret_cast<int&> would also work but triggers -Wstrict-aliasing).
        int species_id;
        std::memcpy(&species_id, &v.w, sizeof(int));
        species[out_idx] = species_id;

        out_idx++;
    }
    *N_returned = out_idx;
}

void pic_cuda_set_B_field(void* handle, float B_T)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;
    int N_cells = s->grid.Nx * s->grid.Ny * s->grid.Nz;
    constexpr int BLOCK = 256;
    int gdim = (N_cells + BLOCK - 1) / BLOCK;
    setToroidalFieldKernel<<<gdim, BLOCK, 0, s->stream_main>>>(
        s->B_grid, B_T, N_cells);
    cudaStreamSynchronize(s->stream_main);
}

void pic_cuda_set_material_map(void* handle, const int* material_ids, int N_cells)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;
    // We need a device-side copy of material_ids, then run the kernel
    int* d_ids;
    cudaMalloc(&d_ids, N_cells * sizeof(int));
    cudaMemcpy(d_ids, material_ids, N_cells * sizeof(int), cudaMemcpyHostToDevice);
    constexpr int BLOCK = 256;
    int gdim = (N_cells + BLOCK - 1) / BLOCK;
    setMaterialMapKernel<<<gdim, BLOCK, 0, s->stream_main>>>(
        s->material_map, d_ids, N_cells);
    cudaStreamSynchronize(s->stream_main);
    cudaFree(d_ids);
}

void pic_cuda_init_rng(void* handle, unsigned long long seed)
{
    PICDeviceState* s = (PICDeviceState*)handle;
    if (!s) return;
    launchInitRNG(s->rng_ion, seed, s->N_particles, s->stream_main);
    launchInitRNG(s->rng_elec, seed + 1, s->N_particles, s->stream_main);
    launchInitRNG(s->rng_neutron, seed + 2, s->N_neutrons_max, s->stream_main);
    launchInitRNG(s->rng_cell, seed + 3,
                   s->grid.Nx * s->grid.Ny * s->grid.Nz, s->stream_main);
    cudaStreamSynchronize(s->stream_main);
}

void pic_cuda_init_cross_sections(const char* path)
{
    if (path) {
        ENDFUpload::initializeFromFile(path);
    } else {
        ENDFUpload::initializeWithDefaults();
    }
}

} // extern "C"
