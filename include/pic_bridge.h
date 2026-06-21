//
// include/pic_bridge.h
// Host-callable interface to the CUDA PIC simulation core.
//
//  This header bridges between PlasmaCoreBridge (host C++) and the CUDA-only
//  SimBuffers / allocSimBuffers / runSimulation functions in plasmacore.cu.
//  It allows the bridge to:
//    1. Allocate / free device buffers
//    2. Initialize particles (set positions, velocities, weights, species)
//    3. Run a single PIC step
//    4. Read back per-cell diagnostics (T_e, n_e, fusion rate)
//    5. Read back particle positions for visualization
//
//  All functions are declared here; their implementations live in
//  pic_bridge.cpp (host) and pic_bridge.cu (CUDA kernels).  The host-side
//  PlasmaCoreBridge only needs the .h, not the .cu files.
//

#pragma once
#include "ReactorState.h"
#include <cstdint>
#include <vector>

// Opaque handle to the device-side SimBuffers + SortContext
struct PICState;

struct PICConfig {
    int   N_particles = 50000;        // macro-particles per species (reduced from 1M for runtime)
    int   Nx = 32, Ny = 32, Nz = 32;  // grid resolution (reduced from 64³ for runtime)
    float Lx = 20.0f, Ly = 20.0f, Lz = 20.0f;  // domain size [m]
    float dt_pic = 1e-12f;            // PIC micro-timestep [s]
    float coulomb_log = 17.0f;
    int   sort_interval = 20;
};

// Per-cell diagnostic snapshot (copied device→host)
struct PICCellDiagnostics {
    std::vector<float> T_e_keV;       // electron temperature per cell
    std::vector<float> n_e_m3;        // electron density per cell
    std::vector<float> B_T;           // magnetic field magnitude per cell
    std::vector<float> q_rad_W_m3;    // radiation power density per cell
    std::vector<float> q_alpha_W_m3;  // alpha heating per cell
    std::vector<float> tau_E_s;       // local confinement time per cell
    int N_cells;
};

// Per-particle snapshot (subset for visualization)
struct PICParticleSnapshot {
    std::vector<float> pos_x, pos_y, pos_z;   // [m]
    std::vector<float> vel_x, vel_y, vel_z;   // [m/s]
    std::vector<float> weight;
    std::vector<int>   species;               // 0=e, 1=D, 2=T, 3=α
    int N_particles;
};

// ─── Lifecycle ────────────────────────────────────────────────────────────────
PICState* picCreate(const PICConfig& cfg);
void picDestroy(PICState* state);

// ─── Initialization ──────────────────────────────────────────────────────────
//
//  Seeds the particle arrays with a Maxwellian distribution at the given
//  (T_e, n_e), with 50:50 D-T mix.  Particles are distributed uniformly in
//  the plasma ellipse (R-R_major)²/a² + (z/κ)² ≤ 1.
//
void picInitializePlasma(PICState* state,
                          float T_e_keV, float n_e_m3,
                          float R_major, float a_minor, float kappa,
                          float B_toroidal_T);

// ─── Single PIC step ─────────────────────────────────────────────────────────
//
//  Runs one complete cycle:
//    sort → deposit → solve fields → gather → push → collide → fuse → transport
//
//  Returns the wall-clock time for the step [ms] (for performance monitoring).
//
float picStep(PICState* state, float dt_pic);

// ─── Diagnostics readback ────────────────────────────────────────────────────
//
//  Copies the per-cell diagnostic arrays from device to host.  This involves
//  a device→host memcpy of ~6 × N_cells floats (for a 32³ grid, ~780 KB).
//
void picReadbackDiagnostics(PICState* state, PICCellDiagnostics& out);

// ─── Particle snapshot for visualization ─────────────────────────────────────
//
//  Copies up to max_particles positions/velocities from device to host.
//  Use this to drive the PARTICLES tab when the real PIC is running.
//
void picReadbackParticles(PICState* state, PICParticleSnapshot& out,
                           int max_particles = 4096);

// ─── Cross-section initialization (from ENDFLoader) ──────────────────────────
//
//  Calls ENDFUpload::initializeWithDefaults() — or a file loader — to
//  populate the __constant__ memory XS tables before neutron transport runs.
//
void picInitializeCrossSections(const char* ace_file_path = nullptr);

// ─── Field boundary conditions ───────────────────────────────────────────────
//
//  Sets the toroidal field B_T on the grid (uniform, aligned with y-axis).
//  Called once after picCreate, or whenever the operator changes B_T.
//
void picSetToroidalField(PICState* state, float B_T);

// ─── Material map for neutron transport ──────────────────────────────────────
//
//  Sets the per-cell material ID (0=vacuum, 1=Li6, 2=Li7, 3=steel, 4=coolant).
//  Default is a thin Li-6/Li-7 blanket shell around the plasma.
//
void picSetMaterialMap(PICState* state, const std::vector<int>& material_ids);
