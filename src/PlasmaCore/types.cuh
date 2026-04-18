#pragma once

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cufft.h>

// ─── Physical Constants ────────────────────────────────────────────────────────
__constant__ float PC_C        = 2.99792458e8f;   // speed of light [m/s]
__constant__ float PC_E        = 1.60217663e-19f; // elementary charge [C]
__constant__ float PC_ME       = 9.10938370e-31f; // electron mass [kg]
__constant__ float PC_MP       = 1.67262192e-27f; // proton mass [kg]
__constant__ float PC_MD       = 3.34449439e-27f; // deuteron mass [kg]
__constant__ float PC_MT       = 5.00735588e-27f; // triton mass [kg]
__constant__ float PC_EPS0     = 8.85418782e-12f; // vacuum permittivity [F/m]
__constant__ float PC_MU0      = 1.25663706e-6f;  // vacuum permeability [H/m]

// ─── Grid Parameters (stored in constant memory) ──────────────────────────────
struct GridParams {
    float ox, oy, oz;       // grid origin [m]
    float dx, dy, dz;       // cell spacing [m]
    int   Nx, Ny, Nz;       // cell counts
};

__constant__ GridParams c_grid;

// ─── Particle Arrays (Structure of Arrays) ────────────────────────────────────
// pos.w  = macro-particle weight (number of real particles represented)
// vel.w  = species ID encoded as float (0=electron, 1=deuterium, 2=tritium, 3=alpha)
// All particle arrays are SoA, NOT AoS, for coalesced access.
struct ParticleArrays {
    float4* pos;        // [x, y, z, weight]     length N
    float4* vel;        // [vx, vy, vz, species]  length N
    int*    cell;       // sorted flat cell index  length N
    int     N;
};

// ─── Reaction Product (alpha + neutron birth) ─────────────────────────────────
struct ReactionProduct {
    float4 pos;         // birth position [x,y,z,weight]
    float4 vel_alpha;   // [vx,vy,vz, species=3(alpha)]
    float4 vel_neutron; // [vx,vy,vz, 0]
    int    active;      // 1 if slot filled, 0 if empty
};

// ─── Neutron Particle ─────────────────────────────────────────────────────────
struct NeutronParticle {
    float3 pos;         // position [m]
    float3 dir;         // unit direction vector
    float  energy_MeV;  // kinetic energy
    float  weight;      // statistical weight
    int    alive;       // 0 = terminated/captured
};

// ─── Cross-Section / TBR / Heat Deposition Maps ───────────────────────────────
struct XSectionTable {
    float* sigma;       // total cross-section [m^2] vs energy bin
    float* E_bins;      // energy bin edges [MeV]
    int    n_bins;
    int    material_id;
};

struct TritiumProductionMap {
    float* tbr_voxel;   // tritium breeding rate per voxel [s^-1]
    int    Nx, Ny, Nz;
};

struct HeatDepositionMap {
    float* q_dot;       // power density [W/m^3] per voxel
    int    Nx, Ny, Nz;
};

// ─── Inline grid helpers ───────────────────────────────────────────────────────
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
