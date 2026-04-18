#pragma once

//
// src/PlasmaCore/PlasmaCoreBridge.h
// Public interface for the PlasmaCore CUDA module.
// Only this header (plus ReactorState.h) is included by other modules.
//

#include "ReactorState.h"
#include <cuda_runtime.h>

struct PlasmaConfig {
    float pic_dt        = 1e-12f;  // PIC micro-timestep [s] (plasma freq limited)
    int   N_particles   = 1000000; // macro-particle count
    int   sort_interval = 20;      // re-sort every N PIC steps
    float coulomb_log   = 15.0f;
    // Grid dimensions (set at construction, not changed at runtime)
    int   Nx = 64, Ny = 64, Nz = 64;
    float Lx = 20.0f, Ly = 20.0f, Lz = 20.0f; // domain size [m]
};

class PlasmaCoreBridge {
public:
    explicit PlasmaCoreBridge(const PlasmaConfig& cfg);
    ~PlasmaCoreBridge();

    // Called once per game tick; reads setpoints from state, writes diagnostics back.
    void update(ReactorState& state, float dt);

    bool isInitialised() const { return initialised_; }

private:
    void initGPU();
    void shutdownGPU();
    void picStep(float dt_pic);
    void readbackDiagnostics(ReactorState& state);

    PlasmaConfig   cfg_;
    bool           initialised_ = false;
    cudaStream_t   stream_main_    = nullptr;
    cudaStream_t   stream_neutron_ = nullptr;
};