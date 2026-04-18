//
// sorting.cu
// Particle sort by cell index and CSR offset computation.
//
//  Sorting by cell is essential for:
//    - Coulomb collision kernel (pairs within same cell)
//    - Fusion reaction kernel  (D-T pairs within same cell)
//    - Field deposition cache performance
//
//  Uses Thrust (sort_by_key) and CUB (ExclusiveSum) — both ship with
//  the CUDA toolkit and impose no external dependencies.
//
//  Sorting strategy:
//    - Compute cell index for every particle from its position.
//    - Sort all particle arrays by that key.
//    - Compute CSR cell_start[] offsets via exclusive prefix sum on counts.
//    - Re-sort every N steps (not every step) — particles don't move far.
//

#include "types.cuh"
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/sequence.h>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_radix_sort.cuh>

// ─── Compute cell index for all particles ─────────────────────────────────────
__global__ void computeCellIndices(
    int* __restrict__ cell_ids,
    const float4* __restrict__ pos,
    GridParams grid,
    int N)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N) return;

    float4 p = pos[pid];
    int ix = (int)((p.x - grid.ox) / grid.dx);
    int iy = (int)((p.y - grid.oy) / grid.dy);
    int iz = (int)((p.z - grid.oz) / grid.dz);

    // Clamp to valid range (particles that exited domain go to cell 0 / boundary)
    ix = max(0, min(ix, grid.Nx - 1));
    iy = max(0, min(iy, grid.Ny - 1));
    iz = max(0, min(iz, grid.Nz - 1));

    cell_ids[pid] = flatCell(ix, iy, iz, grid);
}

// ─── SoA permutation: reorder all particle arrays by sorted index ─────────────
__global__ void permuteParticles(
    float4* __restrict__ pos_out,
    float4* __restrict__ vel_out,
    const float4* __restrict__ pos_in,
    const float4* __restrict__ vel_in,
    const int* __restrict__ perm,   // permutation: new[i] = old[perm[i]]
    int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int src = perm[i];
    pos_out[i] = pos_in[src];
    vel_out[i] = vel_in[src];
}

// ─── Count particles per cell (for CSR offset computation) ────────────────────
__global__ void countParticlesPerCell(
    int* __restrict__ cell_counts,  // output: [N_cells], zeroed before call
    const int* __restrict__ cell_ids,
    int N_particles)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N_particles) return;
    atomicAdd(&cell_counts[cell_ids[pid]], 1);
}

// ─── SortContext: scratch memory for CUB ──────────────────────────────────────
struct SortContext {
    void*  cub_temp   = nullptr;
    size_t cub_bytes  = 0;
    int*   cell_ids   = nullptr;    // [N] particle cell keys (primary)
    int*   cell_ids_alt = nullptr;  // [N] double-buffer
    int*   perm       = nullptr;    // [N] sort permutation
    int*   perm_alt   = nullptr;    // [N] double-buffer
    float4* pos_tmp   = nullptr;    // [N] temporary for permute
    float4* vel_tmp   = nullptr;    // [N] temporary for permute
    int*   counts     = nullptr;    // [N_cells] particles per cell
    int    N_cells    = 0;
    int    N_max      = 0;
};

// ─── Allocate sort scratch memory ────────────────────────────────────────────
void allocSortContext(SortContext& ctx, int N, int N_cells)
{
    ctx.N_max   = N;
    ctx.N_cells = N_cells;

    cudaMalloc(&ctx.cell_ids,     N * sizeof(int));
    cudaMalloc(&ctx.cell_ids_alt, N * sizeof(int));
    cudaMalloc(&ctx.perm,         N * sizeof(int));
    cudaMalloc(&ctx.perm_alt,     N * sizeof(int));
    cudaMalloc(&ctx.pos_tmp,      N * sizeof(float4));
    cudaMalloc(&ctx.vel_tmp,      N * sizeof(float4));
    cudaMalloc(&ctx.counts,       (N_cells + 1) * sizeof(int)); // +1 for scan

    // CUB DeviceRadixSort scratch
    cub::DeviceRadixSort::SortPairs(nullptr, ctx.cub_bytes,
                                    ctx.cell_ids, ctx.cell_ids_alt,
                                    ctx.perm, ctx.perm_alt, N);
    cudaMalloc(&ctx.cub_temp, ctx.cub_bytes);
}

void freeSortContext(SortContext& ctx)
{
    cudaFree(ctx.cell_ids);    cudaFree(ctx.cell_ids_alt);
    cudaFree(ctx.perm);        cudaFree(ctx.perm_alt);
    cudaFree(ctx.pos_tmp);     cudaFree(ctx.vel_tmp);
    cudaFree(ctx.counts);      cudaFree(ctx.cub_temp);
}

// ─── Main Sort + CSR Offset Routine ──────────────────────────────────────────
//
//  After this call:
//    parts.pos / parts.vel are sorted in-place by ascending cell index.
//    parts.cell[i] holds the cell id of particle i.
//    cell_start[c]     = index of first particle in cell c.
//    cell_start[N_cells] = N_particles (sentinel).
//
void sortParticlesByCell(ParticleArrays& parts,
                         int* cell_start,        // [N_cells+1], device
                         SortContext& ctx,
                         GridParams grid,
                         cudaStream_t stream)
{
    int N = parts.N;
    constexpr int BLOCK = 256;
    int gridK = (N + BLOCK - 1) / BLOCK;

    // 1. Compute cell index for each particle
    computeCellIndices<<<gridK, BLOCK, 0, stream>>>(
        ctx.cell_ids, parts.pos, grid, N);

    // 2. Build initial permutation [0, 1, 2, ..., N-1]
    thrust::sequence(thrust::device_ptr<int>(ctx.perm),
                     thrust::device_ptr<int>(ctx.perm) + N, 0);

    // 3. Sort (cell_id, index) pairs by cell_id — CUB radix sort is fastest
    cub::DeviceRadixSort::SortPairs(ctx.cub_temp, ctx.cub_bytes,
                                    ctx.cell_ids, ctx.cell_ids_alt,
                                    ctx.perm,     ctx.perm_alt,
                                    N, 0, 32, stream);

    // After sort: ctx.cell_ids_alt has sorted keys, ctx.perm_alt has new->old map
    // Copy sorted cell ids back to parts.cell
    cudaMemcpyAsync(parts.cell, ctx.cell_ids_alt, N * sizeof(int),
                    cudaMemcpyDeviceToDevice, stream);

    // 4. Permute pos and vel arrays
    permuteParticles<<<gridK, BLOCK, 0, stream>>>(
        ctx.pos_tmp, ctx.vel_tmp,
        parts.pos, parts.vel,
        ctx.perm_alt, N);

    cudaMemcpyAsync(parts.pos, ctx.pos_tmp, N * sizeof(float4),
                    cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(parts.vel, ctx.vel_tmp, N * sizeof(float4),
                    cudaMemcpyDeviceToDevice, stream);

    // 5. Count particles per cell
    cudaMemsetAsync(ctx.counts, 0, (ctx.N_cells + 1) * sizeof(int), stream);
    countParticlesPerCell<<<gridK, BLOCK, 0, stream>>>(
        ctx.counts, ctx.cell_ids_alt, N);

    // 6. Exclusive prefix sum → CSR offsets
    //    cell_start[c] = number of particles in cells [0, c)
    size_t scan_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes,
                                  ctx.counts, cell_start, ctx.N_cells + 1);

    // Reuse cub_temp if large enough, else fall back
    void* scan_buf = ctx.cub_temp;
    if (scan_bytes > ctx.cub_bytes) {
        cudaMalloc(&scan_buf, scan_bytes); // rare: only if scan > sort scratch
    }
    cub::DeviceScan::ExclusiveSum(scan_buf, scan_bytes,
                                  ctx.counts, cell_start, ctx.N_cells + 1,
                                  stream);
    if (scan_buf != ctx.cub_temp) cudaFree(scan_buf);
}
