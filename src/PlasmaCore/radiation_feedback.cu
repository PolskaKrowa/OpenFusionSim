//
// src/PlasmaCore/radiation_feedback.cu
// Couples the radiation kernel back to the local electron temperature.
//
//  Without this feedback, the radiation kernel writes to q_dot_voxel but the
//  temperature grid is unaffected — the simulation doesn't "feel" the
//  radiation losses.  This kernel implements the electron energy equation:
//
//    (3/2) · n_e · dT_e/dt = P_alpha(r) + P_aux(r) - P_brem(r) - P_sync(r) - P_cond(r)
//
//  where the radiation terms are exactly what radiation.cu computes.  The
//  conduction term uses the IPB98(y,2) confinement scaling locally:
//
//    P_cond(r) = (3/2) · n_e(r) · T_e(r) / tau_E(r)
//
//  After this update, the temperature grid is written back to the
//  T_e_per_cell_keV array used by both the radiation kernel (next step)
//  and the visualization tabs.
//
//  References:
//    [1] ITER Physics Basis 1999, §1.5 — power balance.
//    [2] Houlberg, "Plasma transport" in Fusion Physics (IAEA, 2012).
//

#include "types.cuh"
#include <math.h>

// ─── Radiation feedback kernel ───────────────────────────────────────────────
//
//  One thread per cell.  Reads:
//    - n_e_per_cell, T_e_per_cell_keV: local plasma state
//    - B_per_cell_T: for synchrotron
//    - q_alpha_voxel: α heating deposited in this cell [W/m³]
//    - q_rad_voxel: radiation losses from radiation.cu [W/m³]
//    - tau_E_per_cell: local confinement time [s]
//  Updates:
//    - T_e_per_cell_keV: new temperature after dt
//
//  The energy equation per cell (volume V_cell):
//    dU/dt = P_alpha + P_aux - P_rad - P_cond
//    U = (3/2) n_e T_e V_cell
//  so:
//    dT_e/dt = (2/3) · (P_alpha + P_aux - P_rad - P_cond) / (n_e · V_cell)
//    P_cond = U / tau_E = (3/2) n_e T_e V_cell / tau_E
//
__global__ void radiationFeedbackKernel(
    float* __restrict__ T_e_per_cell_keV,    // in/out
    const float* __restrict__ n_e_per_cell,
    const float* __restrict__ B_per_cell_T,
    const float* __restrict__ q_alpha_voxel,  // α heating [W/m³]
    const float* __restrict__ q_rad_voxel,    // radiation loss [W/m³]
    const float* __restrict__ tau_E_per_cell, // local confinement [s]
    float P_aux_W_m3,                          // uniform aux heating [W/m³]
    float dt,
    int N_cells)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N_cells) return;

    float n_e = n_e_per_cell[tid];
    if (n_e < 1.0f) return;   // vacuum cell

    float T_e = T_e_per_cell_keV[tid];
    if (T_e < 0.01f) T_e = 0.01f;

    float P_alpha = q_alpha_voxel[tid];      // W/m³
    float P_rad   = q_rad_voxel[tid];
    float tau_E   = tau_E_per_cell[tid];
    if (tau_E < 1e-6f) tau_E = 1e-6f;

    // Energy density U = (3/2) n_e T_e [J/m³]
    // Convert T from keV to J:  T_J = T_keV × 1.602e-16
    constexpr float KEV_TO_J = 1.602176634e-16f;
    float U = 1.5f * n_e * T_e * KEV_TO_J;

    // Conduction loss: P_cond = U / tau_E
    float P_cond = U / tau_E;

    // Net heating rate [W/m³]
    float dU_dt = P_alpha + P_aux_W_m3 - P_rad - P_cond;

    // Update U
    U += dU_dt * dt;
    if (U < 0.0f) U = 0.0f;

    // Back to T_e [keV]
    float T_new = (2.0f / 3.0f) * U / (n_e * KEV_TO_J);

    // Clamp to physical range
    T_new = fmaxf(fminf(T_new, 100.0f), 0.0f);

    T_e_per_cell_keV[tid] = T_new;
}

// ─── Host launch wrapper ──────────────────────────────────────────────────────
void launchRadiationFeedback(
    float* T_e_per_cell_keV,
    const float* n_e_per_cell,
    const float* B_per_cell_T,
    const float* q_alpha_voxel,
    const float* q_rad_voxel,
    const float* tau_E_per_cell,
    float P_aux_W_m3,
    float dt,
    int N_cells,
    cudaStream_t stream)
{
    constexpr int BLOCK = 256;
    int grid = (N_cells + BLOCK - 1) / BLOCK;
    radiationFeedbackKernel<<<grid, BLOCK, 0, stream>>>(
        T_e_per_cell_keV, n_e_per_cell, B_per_cell_T,
        q_alpha_voxel, q_rad_voxel, tau_E_per_cell,
        P_aux_W_m3, dt, N_cells);
}
