//
// src/PlasmaCore/pic_bridge.cpp
// Host-side implementation of the PIC bridge interface.
//
//  This file lives in the same translation unit as PlasmaCoreBridge.cpp and
//  is compiled by the C++ compiler (not nvcc).  All actual CUDA work is
//  delegated to pic_bridge.cu via the C-linkage functions declared below.
//

#include "pic_bridge.h"
#include "ReactorState.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// ─── Forward-declare the CUDA-side functions (defined in pic_bridge.cu) ──────
extern "C" {

// Allocator: creates SimBuffers + SortContext on device, returns opaque handle
void* pic_cuda_create(int N_particles, int Nx, int Ny, int Nz,
                       float Lx, float Ly, float Lz);
void  pic_cuda_destroy(void* handle);

// Particle initialization
void pic_cuda_init_plasma(void* handle, float T_e_keV, float n_e_m3,
                           float R_major, float a_minor, float kappa,
                           float B_T);

// Single PIC step (returns wall time in ms)
float pic_cuda_step(void* handle, float dt);

// Diagnostic readback
void pic_cuda_readback_diagnostics(void* handle,
                                     float* T_e_keV, float* n_e_m3,
                                     float* B_T, float* q_rad, float* q_alpha,
                                     float* tau_E, int N_cells);

// Particle snapshot
void pic_cuda_readback_particles(void* handle,
                                   float* pos_xyz, float* vel_xyz,
                                   float* weight, int* species,
                                   int max_particles, int* N_returned);

// Field initialization
void pic_cuda_set_B_field(void* handle, float B_T);

// Material map
void pic_cuda_set_material_map(void* handle, const int* material_ids, int N_cells);

// RNG init
void pic_cuda_init_rng(void* handle, unsigned long long seed);

// Cross-section initialization (loads ENDF/B-VIII.0 XS tables to __constant__ memory)
void pic_cuda_init_cross_sections(const char* path);

}

// ─── PICState struct (host-side wrapper) ─────────────────────────────────────
struct PICState {
    void* cuda_handle = nullptr;     // opaque device-side SimBuffers + SortContext
    PICConfig cfg;

    // Cached grid dimensions (for diagnostic buffer sizing)
    int Nx, Ny, Nz, N_cells;

    // Diagnostic buffers (reused across readback calls)
    std::vector<float> T_e_keV;
    std::vector<float> n_e_m3;
    std::vector<float> B_T;
    std::vector<float> q_rad_W_m3;
    std::vector<float> q_alpha_W_m3;
    std::vector<float> tau_E_s;

    // Material map (host-side copy)
    std::vector<int> material_ids;
};

// ─── Lifecycle ────────────────────────────────────────────────────────────────
PICState* picCreate(const PICConfig& cfg)
{
    PICState* state = new PICState{};
    state->cfg = cfg;
    state->Nx = cfg.Nx;
    state->Ny = cfg.Ny;
    state->Nz = cfg.Nz;
    state->N_cells = cfg.Nx * cfg.Ny * cfg.Nz;

    state->cuda_handle = pic_cuda_create(cfg.N_particles, cfg.Nx, cfg.Ny, cfg.Nz,
                                          cfg.Lx, cfg.Ly, cfg.Lz);
    if (!state->cuda_handle) {
        delete state;
        return nullptr;
    }

    // Initialize RNG states
    pic_cuda_init_rng(state->cuda_handle, 0x5EED1234ULL);

    // Allocate diagnostic buffers
    state->T_e_keV.assign(state->N_cells, 0.0f);
    state->n_e_m3.assign(state->N_cells, 0.0f);
    state->B_T.assign(state->N_cells, 0.0f);
    state->q_rad_W_m3.assign(state->N_cells, 0.0f);
    state->q_alpha_W_m3.assign(state->N_cells, 0.0f);
    state->tau_E_s.assign(state->N_cells, 1.0f);
    state->material_ids.assign(state->N_cells, 0);

    return state;
}

void picDestroy(PICState* state)
{
    if (!state) return;
    if (state->cuda_handle) {
        pic_cuda_destroy(state->cuda_handle);
    }
    delete state;
}

// ─── Initialization ──────────────────────────────────────────────────────────
void picInitializePlasma(PICState* state,
                          float T_e_keV, float n_e_m3,
                          float R_major, float a_minor, float kappa,
                          float B_toroidal_T)
{
    if (!state || !state->cuda_handle) return;
    pic_cuda_init_plasma(state->cuda_handle, T_e_keV, n_e_m3,
                          R_major, a_minor, kappa, B_toroidal_T);
    pic_cuda_set_B_field(state->cuda_handle, B_toroidal_T);

    // Initialize a default material map: vacuum inside the plasma, Li-6/Li-7
    // blanket shell around it, steel structure outside.
    //
    // For each cell, compute its position in the (R, z) plane and assign:
    //   - Vacuum   if (R-R_major)²/a² + z²/(κa)² < 1    (inside plasma)
    //   - Li-6     if 1.0 ≤ (above) < 1.2              (inner blanket)
    //   - Li-7     if 1.2 ≤ (above) < 1.4              (outer blanket)
    //   - Steel    if 1.4 ≤ (above) < 2.0              (vacuum vessel)
    //   - Coolant  if (above) ≥ 2.0                    (thermal shield)
    for (int k = 0; k < state->Nz; k++) {
        float z = ((k + 0.5f) / state->Nz - 0.5f) * state->cfg.Lz;
        for (int j = 0; j < state->Ny; j++) {
            float y = ((j + 0.5f) / state->Ny - 0.5f) * state->cfg.Ly;
            for (int i = 0; i < state->Nx; i++) {
                float x = ((i + 0.5f) / state->Nx) * state->cfg.Lx;  // R = x
                float r_norm_sq = (x - R_major) * (x - R_major) / (a_minor * a_minor)
                                + (y * y + z * z) / (kappa * a_minor * kappa * a_minor);
                int mat;
                if (r_norm_sq < 1.0f)         mat = 0;   // vacuum
                else if (r_norm_sq < 1.44f)   mat = 1;   // Li-6
                else if (r_norm_sq < 1.96f)   mat = 2;   // Li-7
                else if (r_norm_sq < 4.0f)    mat = 3;   // steel
                else                          mat = 4;   // coolant
                int idx = i + state->Nx * (j + state->Ny * k);
                state->material_ids[idx] = mat;
            }
        }
    }
    pic_cuda_set_material_map(state->cuda_handle, state->material_ids.data(),
                                state->N_cells);
}

// ─── Single PIC step ─────────────────────────────────────────────────────────
float picStep(PICState* state, float dt_pic)
{
    if (!state || !state->cuda_handle) return 0.0f;
    return pic_cuda_step(state->cuda_handle, dt_pic);
}

// ─── Diagnostics readback ────────────────────────────────────────────────────
void picReadbackDiagnostics(PICState* state, PICCellDiagnostics& out)
{
    if (!state || !state->cuda_handle) return;
    out.N_cells = state->N_cells;
    out.T_e_keV.resize(state->N_cells);
    out.n_e_m3.resize(state->N_cells);
    out.B_T.resize(state->N_cells);
    out.q_rad_W_m3.resize(state->N_cells);
    out.q_alpha_W_m3.resize(state->N_cells);
    out.tau_E_s.resize(state->N_cells);

    pic_cuda_readback_diagnostics(state->cuda_handle,
                                    out.T_e_keV.data(), out.n_e_m3.data(),
                                    out.B_T.data(), out.q_rad_W_m3.data(),
                                    out.q_alpha_W_m3.data(), out.tau_E_s.data(),
                                    state->N_cells);
}

// ─── Particle snapshot ───────────────────────────────────────────────────────
void picReadbackParticles(PICState* state, PICParticleSnapshot& out,
                           int max_particles)
{
    if (!state || !state->cuda_handle) return;

    // Allocate flat buffers for the CUDA readback
    std::vector<float> pos_xyz(3 * max_particles);
    std::vector<float> vel_xyz(3 * max_particles);
    std::vector<float> weight(max_particles);
    std::vector<int>   species(max_particles);

    int N_returned = 0;
    pic_cuda_readback_particles(state->cuda_handle,
                                  pos_xyz.data(), vel_xyz.data(),
                                  weight.data(), species.data(),
                                  max_particles, &N_returned);

    out.N_particles = N_returned;
    out.pos_x.resize(N_returned);
    out.pos_y.resize(N_returned);
    out.pos_z.resize(N_returned);
    out.vel_x.resize(N_returned);
    out.vel_y.resize(N_returned);
    out.vel_z.resize(N_returned);
    out.weight.resize(N_returned);
    out.species.resize(N_returned);

    for (int i = 0; i < N_returned; i++) {
        out.pos_x[i] = pos_xyz[3 * i + 0];
        out.pos_y[i] = pos_xyz[3 * i + 1];
        out.pos_z[i] = pos_xyz[3 * i + 2];
        out.vel_x[i] = vel_xyz[3 * i + 0];
        out.vel_y[i] = vel_xyz[3 * i + 1];
        out.vel_z[i] = vel_xyz[3 * i + 2];
        out.weight[i] = weight[i];
        out.species[i] = species[i];
    }
}

// ─── Cross-section initialization ────────────────────────────────────────────
void picInitializeCrossSections(const char* ace_file_path)
{
    // Delegate to the CUDA-side ENDFUpload::initializeWithDefaults()
    // (or loadFromFile if path is provided).  Declared in the extern "C"
    // block above so the linker resolves it as pic_cuda_init_cross_sections
    // (C mangling) — matching the definition in pic_bridge.cu.
    pic_cuda_init_cross_sections(ace_file_path);
}

// ─── Field boundary conditions ───────────────────────────────────────────────
void picSetToroidalField(PICState* state, float B_T)
{
    if (!state || !state->cuda_handle) return;
    pic_cuda_set_B_field(state->cuda_handle, B_T);
}

// ─── Material map ────────────────────────────────────────────────────────────
void picSetMaterialMap(PICState* state, const std::vector<int>& material_ids)
{
    if (!state || !state->cuda_handle) return;
    if ((int)material_ids.size() != state->N_cells) return;
    state->material_ids = material_ids;
    pic_cuda_set_material_map(state->cuda_handle, state->material_ids.data(),
                                state->N_cells);
}
