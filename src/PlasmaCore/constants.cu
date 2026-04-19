//
// constants.cu
// Single definition point for all __constant__ memory variables.
//
// RULE: every __constant__ variable that is shared across .cu files must be
// DEFINED here and DECLARED extern in types.cuh.  Any .cu file that only
// #includes types.cuh sees the extern declaration and does not emit a symbol.
// nvlink therefore sees exactly one definition of each constant.
//
// To update a constant at runtime use cudaMemcpyToSymbol(), e.g.:
//   float val = 5.3f;
//   cudaMemcpyToSymbol(PC_C, &val, sizeof(float));
//

#include "types.cuh"

// ─── Physical constants ────────────────────────────────────────────────────────
__constant__ float PC_C    = 2.99792458e8f;   // speed of light  [m/s]
__constant__ float PC_E    = 1.60217663e-19f; // elementary charge [C]
__constant__ float PC_ME   = 9.10938370e-31f; // electron mass   [kg]
__constant__ float PC_MP   = 1.67262192e-27f; // proton mass     [kg]
__constant__ float PC_MD   = 3.34449439e-27f; // deuteron mass   [kg]
__constant__ float PC_MT   = 5.00735588e-27f; // triton mass     [kg]
__constant__ float PC_EPS0 = 8.85418782e-12f; // permittivity    [F/m]
__constant__ float PC_MU0  = 1.25663706e-6f;  // permeability    [H/m]

// ─── Grid parameters (set at sim init via cudaMemcpyToSymbol) ────────────────
__constant__ GridParams c_grid;

// ─── Bosch-Hale D-T coefficients (also shared across fusion/collision kernels) ─
// Declared here so fusion_reactions.cu and any future kernel can access without
// re-defining.  Add extern __constant__ BoschHaleCoeffs c_DT_BH; to types.cuh
// if other .cu files need it; for now it is only used in fusion_reactions.cu
// so it stays local there (no extern needed).
