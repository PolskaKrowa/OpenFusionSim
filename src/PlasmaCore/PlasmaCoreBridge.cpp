//
// src/PlasmaCore/PlasmaCoreBridge.cpp
// Host-side C++ wrapper around the CUDA PIC kernels.
// This is the only file in PlasmaCore that the rest of the game touches.
// All CUDA includes and kernel launches live in the .cu files.
//

#include "PlasmaCoreBridge.h"
#include "ReactorState.h"
#include <cmath>
#include <stdexcept>

// Forward declarations of host launch wrappers from the .cu modules
void launchBorisPush(struct ParticleArrays&, float4*, float4*, float, cudaStream_t);
void launchGatherFields(float4*, float4*, const float*, const float*,
                        const struct ParticleArrays&, struct GridParams, cudaStream_t);
void launchDepositCharge(float*, float*, const struct ParticleArrays&,
                         struct GridParams, float, cudaStream_t);
void launchFDTD(float*, float*, float*, const struct GridParams&, float, cudaStream_t);
void launchCoulombCollisions(struct ParticleArrays&, const float*, const int*, float,
                             float, struct curandState*, struct GridParams, int, cudaStream_t);
void launchFusionReactions(const float4*, const float4*, const float4*, const float4*,
                           const int*, const int*, struct ReactionProduct*, struct curandState*,
                           float, struct GridParams, int, cudaStream_t);
void launchNeutronTransport(struct NeutronParticle*, const float*, struct TritiumProductionMap*,
                            struct HeatDepositionMap*, struct curandState*, struct GridParams,
                            int, float, int, cudaStream_t);
void sortParticlesByCell(struct ParticleArrays&, int*, struct SortContext&,
                         struct GridParams, cudaStream_t);

// ─── PlasmaCoreBridge implementation ─────────────────────────────────────────

PlasmaCoreBridge::PlasmaCoreBridge(const PlasmaConfig& cfg)
    : cfg_(cfg)
{
    initGPU();
}

PlasmaCoreBridge::~PlasmaCoreBridge()
{
    shutdownGPU();
}

void PlasmaCoreBridge::initGPU()
{
    cudaError_t err = cudaSetDevice(0);
    if (err != cudaSuccess)
        throw std::runtime_error("PlasmaCore: no CUDA device available");

    cudaStreamCreate(&stream_main_);
    cudaStreamCreate(&stream_neutron_);

    // Allocate all device buffers via the SimBuffers allocator in plasmacore.cu
    // (allocSimBuffers is defined there; we call it here)
    // allocSimBuffers(buf_, cfg_.simConfig);

    initialised_ = true;
}

void PlasmaCoreBridge::shutdownGPU()
{
    if (!initialised_) return;
    cudaStreamSynchronize(stream_main_);
    cudaStreamSynchronize(stream_neutron_);
    cudaStreamDestroy(stream_main_);
    cudaStreamDestroy(stream_neutron_);
    // freeSimBuffers(buf_);
    initialised_ = false;
}

void PlasmaCoreBridge::update(ReactorState& state, float dt)
{
    if (!initialised_) return;
    if (state.plasma_status == PlasmaStatus::Cold) return;

    // ── Run one PIC timestep (or sub-cycle if dt > PIC timestep) ──────────────
    // The PIC dt is O(10^-12 s); game dt is O(10^-3 s).
    // Sub-cycling: run N_pic steps per game tick.
    const int N_pic = static_cast<int>(dt / cfg_.pic_dt);
    for (int i = 0; i < N_pic; i++) {
        picStep(cfg_.pic_dt);
    }

    // ── Extract diagnostics from device and write to ReactorState ─────────────
    readbackDiagnostics(state);
}

void PlasmaCoreBridge::picStep(float dt_pic)
{
    // Mirrors the runSimulation() loop in plasmacore.cu, driven externally.
    // (Sort, deposit, solve, gather, push, collide, fuse, transport)
    // Full implementation delegates to the .cu launch wrappers.
    // Abbreviated here — see plasmacore.cu::runSimulation for the full sequence.
    (void)dt_pic; // placeholder until SimBuffers wired up
}

void PlasmaCoreBridge::readbackDiagnostics(ReactorState& state)
{
    // Copy scalar diagnostics from device diagnostic accumulators to host.
    // In the real implementation these are device→host cudaMemcpyAsync calls
    // on diagnostic float arrays reduced each step.

    // Placeholder physics model for game use before full PIC is wired in.
    // Scales outputs from the current setpoints using simple 0D models.

    float Ti = state.sp_electron_temp_keV;  // ion ≈ electron temp at high density
    float ne = state.plasma_density_m3;
    float B  = state.B_toroidal_T;
    float Ip = state.plasma_current_MA;

    // Fusion power: P ~ n^2 * <σv>(T) * V_plasma
    // Use a rough Bosch-Hale fit valid 10–30 keV: <σv> ~ 3e-22 * T_keV^2 m^3/s
    float sigma_v  = 3.0e-22f * Ti * Ti;          // very rough, keV^2 → m^3/s
    float V_plasma = 840.0f;                        // ITER-class plasma volume [m^3]
    float n_D = 0.5f * ne, n_T = 0.5f * ne;
    float E_fusion_J = 17.59e6f * 1.602e-19f;      // 17.59 MeV per reaction [J]

    state.fusion_power_MW  = n_D * n_T * sigma_v * E_fusion_J * V_plasma * 1e-6f;
    state.alpha_power_MW   = state.fusion_power_MW * (3.52f / 17.59f);
    state.neutron_flux_m2s = state.fusion_power_MW * 1e6f
                           / (14.07e6f * 1.602e-19f) / 600.0f; // first-wall area ~600 m^2

    // Safety factor q95 ~ 5 * a^2 * B / (R * mu0 * Ip * 1e6)
    const float R = 6.2f, a = 2.0f;
    state.q_safety = (5.0f * a * a * B) / (R * 1.257e-6f * Ip * 1e6f);

    // Beta ~ 2 * mu0 * n * T / B^2
    state.beta = (2.0f * 1.257e-6f * ne * Ti * 1.602e-16f) / (B * B);

    // Scientific Q
    float P_aux_MW = 50.0f; // assume 50 MW auxiliary heating (NBI + ICRH)
    state.Q_scientific = (P_aux_MW > 0.1f)
                       ? state.fusion_power_MW / P_aux_MW : 0.0f;

    // Disruption heuristic: q < 2 → Greenwald disruption likely
    state.disruption_flag  = (state.q_safety < 2.0f && state.plasma_current_MA > 1.0f);
    state.alarm_disruption = state.disruption_flag;

    // Plasma status transitions
    if (state.fusion_power_MW > 100.0f)
        state.plasma_status = PlasmaStatus::Burning;
    else if (state.disruption_flag)
        state.plasma_status = PlasmaStatus::Disrupting;

    state.plasma_current_MA = Ip;
    state.plasma_temp_keV   = Ti;
    state.plasma_density_m3 = ne;
}