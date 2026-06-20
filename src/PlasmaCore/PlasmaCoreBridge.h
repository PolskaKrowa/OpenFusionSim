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
    int   N_particles   = 1000000; // macro-particle count (full PIC)
    int   sort_interval = 20;      // re-sort every N PIC steps
    float coulomb_log   = 17.0f;
    // Grid dimensions (set at construction, not changed at runtime)
    int   Nx = 64, Ny = 64, Nz = 64;
    float Lx = 20.0f, Ly = 20.0f, Lz = 20.0f; // domain size [m]
    bool  run_pic = true;          // if false, 0D-only mode (no PIC kernels)
};

// Forward declarations (avoid pulling CUDA headers into other modules)
struct PICState;
namespace MHDPhysics { struct MHDState; }

class PlasmaCoreBridge {
public:
    explicit PlasmaCoreBridge(const PlasmaConfig& cfg);
    ~PlasmaCoreBridge();

    // Called once per game tick; reads setpoints from state, writes diagnostics back.
    void update(ReactorState& state, float dt);

    // Cold-restart: clear the internal 0D power-balance state (stored
    // energy, ion temperature) and force the PIC plasma to be re-seeded
    // from the next non-cold tick.  The PIC grid itself is not destroyed
    // — re-creating it would cost ~50 ms of CUDA work — but the particles
    // are flagged for re-initialisation so a stale hot plasma doesn't
    // bleed across the reset.
    void reset();

    bool isInitialised() const { return initialised_; }

    // Access the MHD state for UI visualization (q-profile, island widths, VDE)
    const MHDPhysics::MHDState* mhdState() const { return mhd_state_; }
    // Access the PIC state for UI visualization (per-cell T_e, particle snapshot)
    PICState* picState() { return pic_state_; }

private:
    void initGPU();
    void shutdownGPU();
    void picStep(float dt_pic);
    void readbackDiagnostics(ReactorState& state, float dt);
    // 0D power-balance model: integrates stored energy via the IPB98(y,2)
    // confinement scaling + bremsstrahlung/synchrotron losses.
    void updatePowerBalance(ReactorState& state, float dt);
    // MHD disruption tracking: 1-D current diffusion + Rutherford island growth
    void updateMHD(ReactorState& state, float dt);

    PlasmaConfig   cfg_;
    bool           initialised_ = false;
    cudaStream_t   stream_main_    = nullptr;
    cudaStream_t   stream_neutron_ = nullptr;

    // ── PIC state (device-side particle arrays + fields) ────────────────────
    PICState*      pic_state_       = nullptr;
    bool           pic_initialized_ = false;

    // ── MHD state (1-D current profile, q-profile, tearing modes) ───────────
    MHDPhysics::MHDState* mhd_state_ = nullptr;

    // ── 0D power-balance integrator state ───────────────────────────────────
    //  These were originally `static float` locals inside updatePowerBalance,
    //  which meant they survived across SCRAM/RESET and a "cold restart"
    //  would start with the pre-SCRAM stored energy still in W_stored_MJ_.
    //  Moved to members so reset() can wipe them.
    //
    //  W_initialized_ is a sticky flag (cleared only by reset()) that tells
    //  updatePowerBalance whether the W_stored_MJ_ integrator has been
    //  initialized from the current T_e/n_e yet.  Without this flag the old
    //  `W_stored_MJ_ < 0.1f` test would re-initialize W every time it
    //  crashed during ramp-down — that re-init was the root cause of the
    //  "fusion power drops to 0 at ~30 MW then bounces back" bug.
    float W_stored_MJ_ = 0.0f;
    float T_i_keV_     = 0.0f;
    bool  W_initialized_ = false;
};
