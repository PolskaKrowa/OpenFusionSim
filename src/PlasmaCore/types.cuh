#pragma once

//
// types.cuh
// Shared types, structs, and DECLARATIONS of __constant__ variables.
//
// __constant__ variables must be DEFINED in exactly one .cu translation unit.
// All other .cu files that include this header see only `extern __constant__`
// declarations, preventing the nvlink "multiple definition" error.
//
// Definitions live in: constants.cu
//

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cufft.h>

// ─── Physical Constants (declared here, defined in constants.cu) ───────────────
extern __constant__ float PC_C;      // speed of light        [m/s]
extern __constant__ float PC_E;      // elementary charge     [C]
extern __constant__ float PC_ME;     // electron mass         [kg]
extern __constant__ float PC_MP;     // proton mass           [kg]
extern __constant__ float PC_MD;     // deuteron mass         [kg]
extern __constant__ float PC_MT;     // triton mass           [kg]
extern __constant__ float PC_EPS0;   // vacuum permittivity   [F/m]
extern __constant__ float PC_MU0;    // vacuum permeability   [H/m]

// ─── Grid Parameters (declared here, defined in constants.cu) ─────────────────
struct GridParams {
    float ox, oy, oz;       // grid origin [m]
    float dx, dy, dz;       // cell spacing [m]
    int   Nx, Ny, Nz;       // cell counts
};

extern __constant__ GridParams c_grid;

// ─── Particle Arrays (Structure of Arrays) ────────────────────────────────────
// pos.w  = macro-particle weight
// vel.w  = species ID encoded as float (0=electron,1=D,2=T,3=alpha)
struct ParticleArrays {
    float4* pos;    // [x, y, z, weight]      length N
    float4* vel;    // [vx, vy, vz, species]  length N
    int*    cell;   // sorted flat cell index  length N
    int     N;
};

// ─── Reaction Product ─────────────────────────────────────────────────────────
struct ReactionProduct {
    float4 pos;
    float4 vel_alpha;
    float4 vel_neutron;
    int    active;
};

// ─── Neutron Particle ─────────────────────────────────────────────────────────
struct NeutronParticle {
    float3 pos;
    float3 dir;
    float  energy_MeV;
    float  weight;
    int    alive;
};

// ─── Cross-Section / TBR / Heat Maps ─────────────────────────────────────────
struct XSectionTable {
    float* sigma;
    float* E_bins;
    int    n_bins;
    int    material_id;
};

struct TritiumProductionMap {
    float* tbr_voxel;
    int    Nx, Ny, Nz;
};

struct HeatDepositionMap {
    float* q_dot;
    int    Nx, Ny, Nz;
};

// ─── Sort context forward declaration ────────────────────────────────────────
struct SortContext;

// ─── Inline grid helpers ──────────────────────────────────────────────────────
__device__ __forceinline__
int flatCell(int ix, int iy, int iz, const GridParams& g) {
    return ix + g.Nx * (iy + g.Ny * iz);
}

__device__ __forceinline__
void worldToCell(float x, float y, float z, const GridParams& g,
                 int& ix, int& iy, int& iz,
                 float& fx, float& fy, float& fz)
{
    fx = (x - g.ox) / g.dx;  ix = (int)fx;  fx -= ix;
    fy = (y - g.oy) / g.dy;  iy = (int)fy;  fy -= iy;
    fz = (z - g.oz) / g.dz;  iz = (int)fz;  fz -= iz;
}
