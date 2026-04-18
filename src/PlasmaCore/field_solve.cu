//
// field_solve.cu
// Two field solvers as described in the spec:
//
//   Option A — Spectral (FFT-based) Poisson solver.
//              Most accurate; suitable for electrostatic / quasi-neutral cases.
//              Uses cuFFT for 3D real-to-complex and complex-to-real transforms.
//
//   Option B — Finite-Difference Time-Domain (FDTD) Maxwell solver.
//              Full electromagnetic; uses the Yee-grid staggering where E and B
//              are offset by half a cell and half a timestep.
//
// Both solvers operate on field grids stored in global memory.
//

#include "types.cuh"
#include <cufft.h>
#include <cmath>

// ═══════════════════════════════════════════════════════════════════════════════
// OPTION A: SPECTRAL POISSON SOLVER
// ═══════════════════════════════════════════════════════════════════════════════

// ─── k-space Poisson solve kernel ─────────────────────────────────────────────
//
//  In k-space:  phi_k = -rho_k / (eps0 * k^2)
//  k^2 = kx^2 + ky^2 + kz^2 (precomputed in k2_grid for speed).
//
//  k2_grid has the same layout as the half-complex cuFFT output:
//    dimensions [Nx, Ny, Nz/2+1], Fortran-order (cuFFT convention).
//
__global__ void poissonSolveKspace(
    cufftComplex* __restrict__ phi_k,
    const cufftComplex* __restrict__ rho_k,
    const float* __restrict__ k2_grid,
    float eps0,
    int Nx, int Ny, int NzHalf)  // NzHalf = Nz/2+1
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = Nx * Ny * NzHalf;
    if (tid >= total) return;

    float k2 = k2_grid[tid];

    if (k2 < 1e-30f) {
        // DC mode: zero potential (no net charge assumed)
        phi_k[tid] = make_cuFloatComplex(0.0f, 0.0f);
    } else {
        float scale = -1.0f / (eps0 * k2);
        phi_k[tid] = make_cuFloatComplex(rho_k[tid].x * scale,
                                         rho_k[tid].y * scale);
    }
}

// ─── Gradient of phi → E field ────────────────────────────────────────────────
//
//  After inverse FFT: E = -grad(phi), computed via finite differences.
//  Central differences, periodic boundary conditions.
//
__global__ void gradPhiToE(
    float* __restrict__ E_grid,     // output: 3 components per node
    const float* __restrict__ phi,  // scalar potential
    float inv2dx, float inv2dy, float inv2dz,
    int Nx, int Ny, int Nz)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= Nx || iy >= Ny || iz >= Nz) return;

    // Periodic neighbours
    int ixp = (ix + 1) % Nx,  ixm = (ix - 1 + Nx) % Nx;
    int iyp = (iy + 1) % Ny,  iym = (iy - 1 + Ny) % Ny;
    int izp = (iz + 1) % Nz,  izm = (iz - 1 + Nz) % Nz;

    auto idx = [&](int x, int y, int z){ return x + Nx * (y + Ny * z); };

    float Ex = -(phi[idx(ixp,iy,iz)] - phi[idx(ixm,iy,iz)]) * inv2dx;
    float Ey = -(phi[idx(ix,iyp,iz)] - phi[idx(ix,iym,iz)]) * inv2dy;
    float Ez = -(phi[idx(ix,iy,izp)] - phi[idx(ix,iy,izm)]) * inv2dz;

    int base = 3 * idx(ix, iy, iz);
    E_grid[base + 0] = Ex;
    E_grid[base + 1] = Ey;
    E_grid[base + 2] = Ez;
}

// ─── k^2 precomputation ───────────────────────────────────────────────────────
__global__ void buildK2Grid(
    float* __restrict__ k2_grid,    // output [Nx * Ny * (Nz/2+1)]
    float Lx, float Ly, float Lz,  // physical domain lengths [m]
    int Nx, int Ny, int Nz)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    int NzH = Nz / 2 + 1;
    if (ix >= Nx || iy >= Ny || iz >= NzH) return;

    // Wave-numbers (negative frequencies for upper half)
    float kx = (ix <= Nx / 2) ? (2.0f * M_PI * ix / Lx)
                               : (2.0f * M_PI * (ix - Nx) / Lx);
    float ky = (iy <= Ny / 2) ? (2.0f * M_PI * iy / Ly)
                               : (2.0f * M_PI * (iy - Ny) / Ly);
    float kz = 2.0f * M_PI * iz / Lz; // iz always in [0, Nz/2]

    k2_grid[ix + Nx * (iy + Ny * iz)] = kx*kx + ky*ky + kz*kz;
}

// ─── Spectral Poisson Host Wrapper ────────────────────────────────────────────
struct SpectralSolver {
    cufftHandle planFwd, planInv;
    cufftComplex* rho_k;
    cufftComplex* phi_k;
    float* k2_grid;
    float* phi_real;
    int Nx, Ny, Nz;
};

void spectralPoissonSolve(SpectralSolver& s,
                          float* rho_grid,
                          float* E_grid,
                          const GridParams& g,
                          cudaStream_t stream)
{
    int NzH = g.Nz / 2 + 1;
    int total_k = g.Nx * g.Ny * NzH;

    // Forward FFT: rho_real -> rho_k
    cufftExecR2C(s.planFwd, rho_grid, s.rho_k);

    // Solve in k-space
    constexpr int BLOCK = 256;
    int grid = (total_k + BLOCK - 1) / BLOCK;
    poissonSolveKspace<<<grid, BLOCK, 0, stream>>>(
        s.phi_k, s.rho_k, s.k2_grid,
        8.85418782e-12f,
        g.Nx, g.Ny, NzH);

    // Inverse FFT: phi_k -> phi_real
    cufftExecC2R(s.planInv, s.phi_k, s.phi_real);

    // Normalize (cuFFT is unnormalized)
    float norm = 1.0f / (g.Nx * g.Ny * g.Nz);
    // (apply norm inside gradPhiToE via scale, or with a cublas dscal)

    // Gradient: phi -> E
    dim3 block3(8, 8, 8);
    dim3 grid3((g.Nx + 7) / 8, (g.Ny + 7) / 8, (g.Nz + 7) / 8);
    gradPhiToE<<<grid3, block3, 0, stream>>>(
        E_grid, s.phi_real,
        0.5f * norm / g.dx,
        0.5f * norm / g.dy,
        0.5f * norm / g.dz,
        g.Nx, g.Ny, g.Nz);
}


// ═══════════════════════════════════════════════════════════════════════════════
// OPTION B: FDTD (YEE GRID) MAXWELL SOLVER
// ═══════════════════════════════════════════════════════════════════════════════
//
//  Yee grid staggering:
//    E components live at cell edges; B components at cell faces.
//    E is updated at integer timesteps; B at half-integer timesteps.
//    This ensures both fields are centered in time (leapfrog).
//
//  Update equations (SI units, no PML/boundary here):
//    B^{n+1/2} = B^{n-1/2} - dt * curl(E^n)
//    E^{n+1}   = E^n       + (dt/eps0) * (curl(B^{n+1/2}) / mu0 - J^{n+1/2})
//

// ─── B update: B^{n-1/2} + dt * curl(E) → B^{n+1/2} ─────────────────────────
__global__ void updateB(
    float* __restrict__ B,       // in/out: 3 components per node
    const float* __restrict__ E,
    float dt,
    int Nx, int Ny, int Nz,
    float inv_dx, float inv_dy, float inv_dz)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= Nx || iy >= Ny || iz >= Nz) return;

    // Periodic next-cell indices
    int ixp = (ix + 1) % Nx;
    int iyp = (iy + 1) % Ny;
    int izp = (iz + 1) % Nz;

    // E component helpers
    auto Ec = [&](int x, int y, int z, int c) -> float {
        return E[3 * (x + Nx * (y + Ny * z)) + c];
    };

    // curl(E) at this node
    float curlEx = (Ec(ix, iy, izp, 1) - Ec(ix, iy, iz, 1)) * inv_dz
                 - (Ec(ix, iyp, iz, 2) - Ec(ix, iy, iz, 2)) * inv_dy;
    float curlEy = (Ec(ixp, iy, iz, 2) - Ec(ix, iy, iz, 2)) * inv_dx
                 - (Ec(ix, iy, izp, 0) - Ec(ix, iy, iz, 0)) * inv_dz;
    float curlEz = (Ec(ix, iyp, iz, 0) - Ec(ix, iy, iz, 0)) * inv_dy
                 - (Ec(ixp, iy, iz, 1) - Ec(ix, iy, iz, 1)) * inv_dx;

    int base = 3 * (ix + Nx * (iy + Ny * iz));
    B[base + 0] -= dt * curlEx;
    B[base + 1] -= dt * curlEy;
    B[base + 2] -= dt * curlEz;
}

// ─── E update: E^n + (dt/eps0) * (curl(B)/mu0 - J) → E^{n+1} ────────────────
__global__ void updateE(
    float* __restrict__ E,       // in/out
    const float* __restrict__ B,
    const float* __restrict__ J, // current density [A/m^2]
    float dt,
    int Nx, int Ny, int Nz,
    float inv_dx, float inv_dy, float inv_dz)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= Nx || iy >= Ny || iz >= Nz) return;

    int ixm = (ix - 1 + Nx) % Nx;
    int iym = (iy - 1 + Ny) % Ny;
    int izm = (iz - 1 + Nz) % Nz;

    auto Bc = [&](int x, int y, int z, int c) -> float {
        return B[3 * (x + Nx * (y + Ny * z)) + c];
    };

    // curl(B) at this node
    float curlBx = (Bc(ix, iy, iz, 2) - Bc(ix, iym, iz, 2)) * inv_dy
                 - (Bc(ix, iy, iz, 1) - Bc(ix, iy, izm, 1)) * inv_dz;
    float curlBy = (Bc(ix, iy, iz, 0) - Bc(ix, iy, izm, 0)) * inv_dz
                 - (Bc(ix, iy, iz, 2) - Bc(ixm, iy, iz, 2)) * inv_dx;
    float curlBz = (Bc(ix, iy, iz, 1) - Bc(ixm, iy, iz, 1)) * inv_dx
                 - (Bc(ix, iy, iz, 0) - Bc(ix, iym, iz, 0)) * inv_dy;

    constexpr float mu0    = 1.25663706e-6f;
    constexpr float eps0   = 8.85418782e-12f;
    float coeff = dt / eps0;

    int base = 3 * (ix + Nx * (iy + Ny * iz));
    E[base + 0] += coeff * (curlBx / mu0 - J[base + 0]);
    E[base + 1] += coeff * (curlBy / mu0 - J[base + 1]);
    E[base + 2] += coeff * (curlBz / mu0 - J[base + 2]);
}

// ─── FDTD Host Wrapper ────────────────────────────────────────────────────────
void launchFDTD(float* E_grid, float* B_grid, float* J_grid,
                const GridParams& g, float dt,
                cudaStream_t stream)
{
    dim3 block(8, 8, 8);
    dim3 grid((g.Nx + 7) / 8, (g.Ny + 7) / 8, (g.Nz + 7) / 8);

    float idx = 1.0f / g.dx, idy = 1.0f / g.dy, idz = 1.0f / g.dz;

    // Leapfrog: B at half-step first, then E at full step
    updateB<<<grid, block, 0, stream>>>(B_grid, E_grid, dt,
                                        g.Nx, g.Ny, g.Nz,
                                        idx, idy, idz);
    updateE<<<grid, block, 0, stream>>>(E_grid, B_grid, J_grid, dt,
                                        g.Nx, g.Ny, g.Nz,
                                        idx, idy, idz);
}
