//
// fusion_data.cu
// Single definition point for the __constant__ Bosch-Hale fusion tables.
//
//  Pattern matches constants.cu: the __constant__ symbol is declared extern
//  in fusion_data.cuh (so other .cu files can READ it) and DEFINED here
//  (exactly once across the entire binary, so nvlink sees one definition).
//

#include "fusion_data.cuh"

// ─── Verified D-T, D-D(n), D-D(p), D-³He tables ───────────────────────────────
//  Coefficients from Bosch-Hale 1992 Table VII, cross-checked vs UWFDM-1268.
//
//  D-T:   mrc2 = 1 124 656 keV  (D-T reduced mass · c² ≈ 1124.7 MeV)
//         BG   = 34.3827 keV^1/2
//  D-D:   mrc2 =   937 814 keV, BG = 31.3970 keV^1/2
//  D-³He: mrc2 = 1 124 572 keV, BG = 68.7508 keV^1/2
//
__constant__ BoschHaleTable c_fusion_tables[N_FUSION_CHANNELS] = {
    // ── D-T ──────────────────────────────────────────────────────────────────
    {
        34.3827f,                                       // BG
        1124656.0f,                                     // mrc2 (keV) — NOT 1124.6!
        1.17302e-9f,                                    // C1
        1.51361e-2f, 7.51886e-2f,                       // C2, C3
        4.60643e-3f, 1.35000e-2f,                       // C4, C5
       -1.06750e-4f, 1.36600e-5f,                       // C6, C7
        17.589f,                                        // Q
        3.521f, 14.070f,                                // E_α, E_n
        0.2f, 100.0f                                    // valid T range
    },
    // ── D(d,n)³He ────────────────────────────────────────────────────────────
    {
        31.3970f,
        937814.0f,
        5.43360e-12f,
        5.85778e-3f, 7.68222e-3f,
        0.0f,        -2.96400e-6f,
        0.0f,         0.0f,
        3.269f,
        0.820f, 2.449f,                                 // E(³He), E(n)  in CoM
        0.2f, 100.0f
    },
    // ── D(d,p)T  (produces tritium — closes the fuel cycle!) ─────────────────
    {
        31.3970f,
        937814.0f,
        5.65718e-12f,
        3.41267e-3f, 1.99167e-3f,
        0.0f,         1.05060e-5f,
        0.0f,         0.0f,
        4.033f,
        3.016f, 1.017f,                                 // E(p), E(T)
        0.2f, 100.0f
    },
    // ── ³He(d,p)⁴He ──────────────────────────────────────────────────────────
    {
        68.7508f,
        1124572.0f,
        5.51036e-10f,
        6.41918e-3f, -2.02896e-3f,
       -1.91080e-5f,  1.35776e-4f,
        0.0f,         0.0f,
        18.353f,
        14.700f, 3.653f,                                // E(p), E(α)
        0.5f, 190.0f
    }
};
