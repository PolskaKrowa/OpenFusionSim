//
// field_interpolation.cu
// Gathers E and B fields from the grid to particle positions.
// Uses trilinear interpolation over the 8 surrounding cell corners.
// Particles are scattered in space → cache-hostile; shared memory
// staging is used to amortise global memory latency per tile of particles.
//

#include "types.cuh"

// ─── Constants ─────────────────────────────────────────────────────────────────
// E_grid and B_grid layout: 3 components interleaved per cell node.
//   index for component c at node (ix,iy,iz):
//     3 * flatCell(ix,iy,iz) + c
//   c: 0=x, 1=y, 2=z

// ─── Device helper: trilinear fetch of a 3-component vector at a grid node ────
__device__ __forceinline__
float3 fetchGridVec(const float* __restrict__ grid, int ix, int iy, int iz,
                    const GridParams& g)
{
    // Clamp to valid range (absorbing boundary)
    ix = max(0, min(ix, g.Nx - 1));
    iy = max(0, min(iy, g.Ny - 1));
    iz = max(0, min(iz, g.Nz - 1));
    int base = 3 * flatCell(ix, iy, iz, g);
    return make_float3(grid[base], grid[base + 1], grid[base + 2]);
}

// ─── Device helper: trilinear interpolation ───────────────────────────────────
__device__ __forceinline__
float3 trilinear(const float* __restrict__ grid,
                 int ix, int iy, int iz,
                 float fx, float fy, float fz,
                 const GridParams& g)
{
    // Weights for 8 corners: (1-f)(1-f)(1-f) pattern
    float3 c000 = fetchGridVec(grid, ix,   iy,   iz,   g);
    float3 c100 = fetchGridVec(grid, ix+1, iy,   iz,   g);
    float3 c010 = fetchGridVec(grid, ix,   iy+1, iz,   g);
    float3 c110 = fetchGridVec(grid, ix+1, iy+1, iz,   g);
    float3 c001 = fetchGridVec(grid, ix,   iy,   iz+1, g);
    float3 c101 = fetchGridVec(grid, ix+1, iy,   iz+1, g);
    float3 c011 = fetchGridVec(grid, ix,   iy+1, iz+1, g);
    float3 c111 = fetchGridVec(grid, ix+1, iy+1, iz+1, g);

    float omfx = 1.0f - fx;
    float omfy = 1.0f - fy;
    float omfz = 1.0f - fz;

    // Blend along x
    float3 c00 = make_float3(omfx * c000.x + fx * c100.x,
                             omfx * c000.y + fx * c100.y,
                             omfx * c000.z + fx * c100.z);
    float3 c10 = make_float3(omfx * c010.x + fx * c110.x,
                             omfx * c010.y + fx * c110.y,
                             omfx * c010.z + fx * c110.z);
    float3 c01 = make_float3(omfx * c001.x + fx * c101.x,
                             omfx * c001.y + fx * c101.y,
                             omfx * c001.z + fx * c101.z);
    float3 c11 = make_float3(omfx * c011.x + fx * c111.x,
                             omfx * c011.y + fx * c111.y,
                             omfx * c011.z + fx * c111.z);

    // Blend along y
    float3 c0  = make_float3(omfy * c00.x + fy * c10.x,
                             omfy * c00.y + fy * c10.y,
                             omfy * c00.z + fy * c10.z);
    float3 c1  = make_float3(omfy * c01.x + fy * c11.x,
                             omfy * c01.y + fy * c11.y,
                             omfy * c01.z + fy * c11.z);

    // Blend along z
    return make_float3(omfz * c0.x + fz * c1.x,
                       omfz * c0.y + fz * c1.y,
                       omfz * c0.z + fz * c1.z);
}

// ─── Gather Kernel ────────────────────────────────────────────────────────────
//
//  For each particle: compute its cell and fractional offsets, then
//  trilinearly interpolate E and B from the surrounding 8 nodes.
//
//  Shared memory staging: each block loads a tile of particle positions into
//  shared memory first to reduce register pressure and allow some prefetching
//  before the scattered grid reads.
//
#define GATHER_BLOCK 128

__global__ void gatherFields(
    float4* __restrict__ E_at_particle,         // output: E at particle positions
    float4* __restrict__ B_at_particle,         // output: B at particle positions
    const float* __restrict__ E_grid,           // 3D E grid, 3 floats per node
    const float* __restrict__ B_grid,           // 3D B grid, 3 floats per node
    const float4* __restrict__ pos,             // particle positions [x,y,z,w]
    GridParams grid,
    int N_particles)
{
    __shared__ float4 s_pos[GATHER_BLOCK];

    int pid_global = blockIdx.x * GATHER_BLOCK + threadIdx.x;

    // Stage positions into shared memory
    if (pid_global < N_particles)
        s_pos[threadIdx.x] = pos[pid_global];
    __syncthreads();

    if (pid_global >= N_particles) return;

    float4 p = s_pos[threadIdx.x];

    // Compute cell indices and fractional offsets
    int ix, iy, iz;
    float fx, fy, fz;
    worldToCell(p.x, p.y, p.z, grid, ix, iy, iz, fx, fy, fz);

    // Gather fields
    float3 E = trilinear(E_grid, ix, iy, iz, fx, fy, fz, grid);
    float3 B = trilinear(B_grid, ix, iy, iz, fx, fy, fz, grid);

    // Write results (w component unused, set 0)
    E_at_particle[pid_global] = make_float4(E.x, E.y, E.z, 0.0f);
    B_at_particle[pid_global] = make_float4(B.x, B.y, B.z, 0.0f);
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchGatherFields(float4* E_at_particle,
                        float4* B_at_particle,
                        const float* E_grid,
                        const float* B_grid,
                        const ParticleArrays& parts,
                        GridParams grid,
                        cudaStream_t stream)
{
    int gridDim = (parts.N + GATHER_BLOCK - 1) / GATHER_BLOCK;
    gatherFields<<<gridDim, GATHER_BLOCK, 0, stream>>>(
        E_at_particle, B_at_particle,
        E_grid, B_grid,
        parts.pos, grid, parts.N);
}
