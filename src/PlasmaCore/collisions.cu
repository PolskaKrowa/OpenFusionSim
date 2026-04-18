//
// collisions.cu
// Takizuka-Abe binary Coulomb collision operator.
//
//  Particles must be pre-sorted by cell (see sorting.cu).
//  Within each cell, random pairs are selected and a random scattering
//  angle is applied, conserving energy and momentum locally.
//
//  cuRAND provides per-thread Gaussian RNG; states are initialised once
//  at simulation startup and reused every timestep.
//

#include "types.cuh"
#include <curand_kernel.h>
#include <math.h>

// ─── cuRAND state initialisation ──────────────────────────────────────────────
__global__ void initRNGStates(curandState* states, unsigned long long seed, int N)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;
    // Each thread gets a different sequence to avoid correlations
    curand_init(seed, tid, 0, &states[tid]);
}

// ─── Species mass helper ──────────────────────────────────────────────────────
__device__ __forceinline__
float speciesMass(int s)
{
    switch (s) {
        case 0:  return 9.10938370e-31f;             // electron
        case 1:  return 3.34449439e-27f;             // deuterium
        case 2:  return 5.00735588e-27f;             // tritium
        case 3:  return 4.0f * 1.67262192e-27f;      // alpha
        default: return 1.67262192e-27f;
    }
}

// ─── Scatter a pair of particles by angle (theta, phi) ───────────────────────
//  Applies random deflection in the centre-of-mass frame then back-transforms.
//  Implements the Takizuka-Abe scattering step.
__device__ void scatterPair(
    float3& va, float ma,   // velocity & mass of particle A
    float3& vb, float mb,   // velocity & mass of particle B
    float sin_theta, float cos_theta, float phi)
{
    float M  = ma + mb;
    float mu = ma * mb / M; // reduced mass

    // Centre-of-mass velocity
    float3 vcm = make_float3((ma * va.x + mb * vb.x) / M,
                             (ma * va.y + mb * vb.y) / M,
                             (ma * va.z + mb * vb.z) / M);

    // Relative velocity in CoM frame
    float3 vrel = make_float3(va.x - vb.x, va.y - vb.y, va.z - vb.z);
    float vrel_mag = sqrtf(vrel.x*vrel.x + vrel.y*vrel.y + vrel.z*vrel.z);
    if (vrel_mag < 1e-10f) return; // degenerate pair

    // Build a perpendicular basis to vrel
    // Choose an arbitrary vector not parallel to vrel
    float3 perp1;
    if (fabsf(vrel.x) < fabsf(vrel.y)) {
        perp1 = make_float3(0.0f, -vrel.z, vrel.y);
    } else {
        perp1 = make_float3(-vrel.z, 0.0f, vrel.x);
    }
    float p1mag = sqrtf(perp1.x*perp1.x + perp1.y*perp1.y + perp1.z*perp1.z);
    perp1.x /= p1mag; perp1.y /= p1mag; perp1.z /= p1mag;

    // perp2 = vrel/|vrel| x perp1
    float3 vhat = make_float3(vrel.x / vrel_mag, vrel.y / vrel_mag, vrel.z / vrel_mag);
    float3 perp2 = make_float3(vhat.y * perp1.z - vhat.z * perp1.y,
                               vhat.z * perp1.x - vhat.x * perp1.z,
                               vhat.x * perp1.y - vhat.y * perp1.x);

    float cos_phi = cosf(phi), sin_phi = sinf(phi);

    // Scattered relative velocity (same magnitude, rotated direction)
    float3 vrel_new = make_float3(
        vrel_mag * (sin_theta * cos_phi * perp1.x + sin_theta * sin_phi * perp2.x + cos_theta * vhat.x),
        vrel_mag * (sin_theta * cos_phi * perp1.y + sin_theta * sin_phi * perp2.y + cos_theta * vhat.y),
        vrel_mag * (sin_theta * cos_phi * perp1.z + sin_theta * sin_phi * perp2.z + cos_theta * vhat.z));

    // Back-transform to lab frame
    float frac_b = mb / M;
    float frac_a = ma / M;
    va = make_float3(vcm.x + frac_b * vrel_new.x,
                     vcm.y + frac_b * vrel_new.y,
                     vcm.z + frac_b * vrel_new.z);
    vb = make_float3(vcm.x - frac_a * vrel_new.x,
                     vcm.y - frac_a * vrel_new.y,
                     vcm.z - frac_a * vrel_new.z);
}

// ─── Coulomb Collision Kernel ─────────────────────────────────────────────────
//
//  Particles must be sorted by cell; cell_start[c] gives the index of the
//  first particle in cell c (CSR format), and cell_start[Ncells] = N_particles.
//
//  One thread block handles one cell.  Within the block, threads pair up
//  particles sequentially.  For cells with more particles than threads,
//  the loop repeats with stride = blockDim.x.
//
__global__ void coulombCollisions(
    float4* __restrict__ vel,           // in/out [vx,vy,vz,species]
    const float* __restrict__ weights,  // macro-particle weights
    const int* __restrict__ cell_start, // CSR cell offsets [N_cells+1]
    float dt,
    float coulomb_log,                  // ln(Lambda) ≈ 10–20 for fusion plasma
    curandState* rng_states,
    GridParams grid,
    int N_cells)
{
    int cell_id = blockIdx.x;
    if (cell_id >= N_cells) return;

    int pstart = cell_start[cell_id];
    int pend   = cell_start[cell_id + 1];
    int n_cell = pend - pstart;
    if (n_cell < 2) return;

    int tid = threadIdx.x;
    curandState local_rng = rng_states[pstart + (tid % n_cell)];

    // Cell volume for density estimation
    float V_cell = grid.dx * grid.dy * grid.dz;

    // Each thread handles one pair per iteration
    for (int i = tid; i < n_cell / 2; i += blockDim.x) {
        int pidA = pstart + 2 * i;
        int pidB = pstart + 2 * i + 1;

        float4 vA4 = vel[pidA];
        float4 vB4 = vel[pidB];

        int sA = __float_as_int(vA4.w);
        int sB = __float_as_int(vB4.w);
        float mA = speciesMass(sA);
        float mB = speciesMass(sB);
        float wA = weights[pidA];
        float wB = weights[pidB];

        float3 vA = make_float3(vA4.x, vA4.y, vA4.z);
        float3 vB = make_float3(vB4.x, vB4.y, vB4.z);

        // Relative speed
        float dvx = vA.x - vB.x, dvy = vA.y - vB.y, dvz = vA.z - vB.z;
        float vrel = sqrtf(dvx*dvx + dvy*dvy + dvz*dvz);
        if (vrel < 1.0f) continue; // skip near-zero relative velocity

        // Scattering angle variance: <delta_u^2> = Gamma * dt / (vrel)
        // where Gamma = n * q_a^2 * q_b^2 * ln(Lambda) / (4 pi eps0^2 * mu^2)
        // Approximate with particle density n ~ w/V_cell
        float n_eff = (wA + wB) * 0.5f / V_cell;
        float qA = (sA == 0) ? -PC_E : (sA == 3 ? 2.0f * PC_E : PC_E);
        float qB = (sB == 0) ? -PC_E : (sB == 3 ? 2.0f * PC_E : PC_E);
        float mu_red = mA * mB / (mA + mB);
        float eps0 = 8.85418782e-12f;
        float pi   = 3.14159265f;

        float Gamma = n_eff * qA*qA * qB*qB * coulomb_log
                    / (4.0f * pi * eps0 * eps0 * mu_red * mu_red);

        float delta_u2 = Gamma * dt / (vrel * vrel * vrel);

        // Gaussian random scattering angle
        float U1 = curand_normal(&local_rng);
        float U2 = curand_normal(&local_rng);
        float dphi = 2.0f * pi * curand_uniform(&local_rng);

        // sin/cos of scattering angle from delta_u2
        float cos_theta = 1.0f - 0.5f * delta_u2 * (1.0f - expf(-delta_u2 / 2.0f));
        cos_theta = fmaxf(-1.0f, fminf(1.0f, cos_theta));
        float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);

        scatterPair(vA, mA, vB, mB, sin_theta, cos_theta, dphi);

        vel[pidA] = make_float4(vA.x, vA.y, vA.z, vA4.w);
        vel[pidB] = make_float4(vB.x, vB.y, vB.z, vB4.w);
    }

    // Store updated RNG state
    rng_states[pstart + (tid % n_cell)] = local_rng;
}

// ─── Host Launch Wrappers ──────────────────────────────────────────────────────
void launchInitRNG(curandState* states, unsigned long long seed, int N,
                   cudaStream_t stream)
{
    constexpr int BLOCK = 256;
    int grid = (N + BLOCK - 1) / BLOCK;
    initRNGStates<<<grid, BLOCK, 0, stream>>>(states, seed, N);
}

void launchCoulombCollisions(ParticleArrays& parts,
                             const float* weights,
                             const int* cell_start,
                             float dt,
                             float coulomb_log,
                             curandState* rng_states,
                             GridParams grid,
                             int N_cells,
                             cudaStream_t stream)
{
    // One block per cell; 64 threads sufficient for most cell populations
    constexpr int THREADS_PER_CELL = 64;
    coulombCollisions<<<N_cells, THREADS_PER_CELL, 0, stream>>>(
        parts.vel, weights, cell_start,
        dt, coulomb_log, rng_states, grid, N_cells);
}
