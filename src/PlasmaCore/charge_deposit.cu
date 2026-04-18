//
// charge_deposit.cu
// Deposits charge density (rho) and current density (J) from particles
// to the grid using trilinear weighting (first-order cloud-in-cell).
//
// Race conditions: multiple particles deposit to the same cell corners.
// Strategy:
//   - Primary:  sort particles by cell, then use warp-level reductions
//               to accumulate within a warp before a single atomicAdd.
//   - Fallback: direct atomicAdd per corner (native float on Ampere+).
//
// This file implements the direct atomicAdd approach, which is simpler
// and performant on Ampere and newer architectures where native float
// atomicAdd is hardware-accelerated.
//

#include "types.cuh"

// ─── Deposit helper: distribute a scalar quantity to 8 grid corners ───────────
__device__ __forceinline__
void depositToCorners(float* __restrict__ grid_scalar, // target scalar grid
                      int ix, int iy, int iz,
                      float fx, float fy, float fz,
                      float value,
                      const GridParams& g)
{
    float omfx = 1.0f - fx;
    float omfy = 1.0f - fy;
    float omfz = 1.0f - fz;

    // 8 corner weights
    float w[2][2][2];
    w[0][0][0] = omfx * omfy * omfz;
    w[1][0][0] =   fx * omfy * omfz;
    w[0][1][0] = omfx *   fy * omfz;
    w[1][1][0] =   fx *   fy * omfz;
    w[0][0][1] = omfx * omfy *   fz;
    w[1][0][1] =   fx * omfy *   fz;
    w[0][1][1] = omfx *   fy *   fz;
    w[1][1][1] =   fx *   fy *   fz;

    // Deposit with atomicAdd to avoid race conditions
    for (int dz = 0; dz <= 1; dz++) {
        int cz = iz + dz;
        if (cz < 0 || cz >= g.Nz) continue;
        for (int dy = 0; dy <= 1; dy++) {
            int cy = iy + dy;
            if (cy < 0 || cy >= g.Ny) continue;
            for (int dx = 0; dx <= 1; dx++) {
                int cx = ix + dx;
                if (cx < 0 || cx >= g.Nx) continue;
                int idx = flatCell(cx, cy, cz, g);
                atomicAdd(&grid_scalar[idx], value * w[dx][dy][dz]);
            }
        }
    }
}

// ─── Deposit helper: distribute a 3-vector to 8 grid corners ──────────────────
__device__ __forceinline__
void depositToCorners3(float* __restrict__ grid_vec3, // target: 3 floats per cell
                       int ix, int iy, int iz,
                       float fx, float fy, float fz,
                       float jx, float jy, float jz,
                       const GridParams& g)
{
    float omfx = 1.0f - fx;
    float omfy = 1.0f - fy;
    float omfz = 1.0f - fz;

    for (int dz = 0; dz <= 1; dz++) {
        int cz = iz + dz;
        float wz = (dz == 0) ? omfz : fz;
        if (cz < 0 || cz >= g.Nz) continue;

        for (int dy = 0; dy <= 1; dy++) {
            int cy = iy + dy;
            float wy = ((dy == 0) ? omfy : fy) * wz;
            if (cy < 0 || cy >= g.Ny) continue;

            for (int dx = 0; dx <= 1; dx++) {
                int cx = ix + dx;
                float w = ((dx == 0) ? omfx : fx) * wy;
                if (cx < 0 || cx >= g.Nx) continue;

                int base = 3 * flatCell(cx, cy, cz, g);
                atomicAdd(&grid_vec3[base + 0], jx * w);
                atomicAdd(&grid_vec3[base + 1], jy * w);
                atomicAdd(&grid_vec3[base + 2], jz * w);
            }
        }
    }
}

// ─── Charge & Current Deposition Kernel ───────────────────────────────────────
//
//  rho_grid: flat cell array [1 float per cell], units: C/m^3
//  J_grid:   flat cell array [3 floats per cell, interleaved], units: A/m^2
//
__global__ void depositCharge(
    float* __restrict__ rho_grid,
    float* __restrict__ J_grid,
    const float4* __restrict__ pos,     // [x,y,z, weight]
    const float4* __restrict__ vel,     // [vx,vy,vz, species]
    GridParams grid,
    float dt,
    int N_particles)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N_particles) return;

    float4 p = pos[pid];
    float4 v = vel[pid];

    float weight = p.w;
    int species  = __float_as_int(v.w);

    // Charge per real particle for this species
    float q_real;
    switch (species) {
        case 0:  q_real = -PC_E; break;         // electron
        case 1:  q_real =  PC_E; break;         // deuterium (Z=1)
        case 2:  q_real =  PC_E; break;         // tritium   (Z=1)
        case 3:  q_real =  2.0f * PC_E; break;  // alpha     (Z=2)
        default: q_real =  PC_E; break;
    }

    // Volume of a cell [m^3]
    float cell_vol = grid.dx * grid.dy * grid.dz;

    // Charge density contribution: rho = (N_real * q) / V_cell
    float q_macro = weight * q_real / cell_vol;

    // Current density: J = rho * v (pre-volume-normalised)
    float jx = q_macro * v.x;
    float jy = q_macro * v.y;
    float jz = q_macro * v.z;

    // Cell index + fractional offsets
    int ix, iy, iz;
    float fx, fy, fz;
    worldToCell(p.x, p.y, p.z, grid, ix, iy, iz, fx, fy, fz);

    // Deposit rho and J to 8 surrounding corners
    depositToCorners (rho_grid, ix, iy, iz, fx, fy, fz, q_macro, grid);
    depositToCorners3(J_grid,   ix, iy, iz, fx, fy, fz, jx, jy, jz, grid);
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchDepositCharge(float* rho_grid,
                         float* J_grid,
                         const ParticleArrays& parts,
                         GridParams grid,
                         float dt,
                         cudaStream_t stream)
{
    // Zero the grids before deposition each step
    int n_cells = grid.Nx * grid.Ny * grid.Nz;
    cudaMemsetAsync(rho_grid, 0, n_cells * sizeof(float), stream);
    cudaMemsetAsync(J_grid,   0, n_cells * 3 * sizeof(float), stream);

    constexpr int BLOCK = 256;
    int gridDim = (parts.N + BLOCK - 1) / BLOCK;
    depositCharge<<<gridDim, BLOCK, 0, stream>>>(
        rho_grid, J_grid,
        parts.pos, parts.vel,
        grid, dt, parts.N);
}
