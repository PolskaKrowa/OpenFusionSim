//
// src/PlasmaCore/endf_upload.cu
// Uploads cross-section data from host (ENDFLoader::CrossSectionSet) to the
// __constant__ memory symbols used by neutron_transport.cu.
//
//  Usage:
//    ENDFLoader::CrossSectionSet xs;
//    ENDFLoader::loadDefaultCrossSections(xs);   // or loadFromTextFile()
//    uploadCrossSectionsToGPU(xs);
//
//  After this call, neutron_transport.cu's kernels will use the loaded XS
//  tables instead of the (zero-initialized) defaults.
//

#include "types.cuh"
#include "endf_loader.h"
#include <cuda_runtime.h>
#include <vector>

// Match the __constant__ declarations in neutron_transport.cu (which uses
// extern __constant__ for these symbols).  We provide the definitions here.
// Use ENDFLoader::N_MAT and ENDFLoader::N_GROUPS to avoid namespace issues.
using ENDFLoader::N_MAT;
using ENDFLoader::N_GROUPS;
using ENDFLoader::MAT_LI6;
using ENDFLoader::MAT_LI7;

extern __constant__ float c_sigma_total[N_MAT][N_GROUPS];
extern __constant__ float c_sigma_abs  [N_MAT][N_GROUPS];
extern __constant__ float c_sigma_Li6  [N_GROUPS];
extern __constant__ float c_sigma_Li7  [N_GROUPS];
extern __constant__ float c_sigma_major[N_GROUPS];
extern __constant__ float c_E_bins[N_GROUPS + 1];

namespace ENDFUpload {

// ─── Upload cross-section tables to GPU __constant__ memory ──────────────────
void uploadCrossSections(const ENDFLoader::CrossSectionSet& xs)
{
    if (!xs.loaded) return;

    // ── Energy bin boundaries ───────────────────────────────────────────────
    cudaMemcpyToSymbol(c_E_bins, xs.E_bins.data(),
                       (N_GROUPS + 1) * sizeof(float));

    // ── Per-material total / absorption cross-sections ──────────────────────
    //  The __constant__ arrays are declared as [N_MAT][N_GROUPS] (row-major),
    //  so we need to copy each material's row separately.
    for (int m = 0; m < N_MAT; m++) {
        // Total XS row
        std::vector<float> sigma_total_row(N_GROUPS);
        for (int g = 0; g < N_GROUPS; g++) {
            sigma_total_row[g] = xs.materials[m].sigma_total[g];
        }
        cudaMemcpyToSymbol(c_sigma_total, sigma_total_row.data(),
                           N_GROUPS * sizeof(float),
                           m * N_GROUPS * sizeof(float));

        // Absorption XS row
        std::vector<float> sigma_abs_row(N_GROUPS);
        for (int g = 0; g < N_GROUPS; g++) {
            sigma_abs_row[g] = xs.materials[m].sigma_abs[g];
        }
        cudaMemcpyToSymbol(c_sigma_abs, sigma_abs_row.data(),
                           N_GROUPS * sizeof(float),
                           m * N_GROUPS * sizeof(float));
    }

    // ── Li-6 (n,t)α XS (only material 1 has this) ───────────────────────────
    cudaMemcpyToSymbol(c_sigma_Li6, xs.materials[MAT_LI6].sigma_Li6_nt.data(),
                       N_GROUPS * sizeof(float));

    // ── Li-7 (n,n')t XS (only material 2 has this) ──────────────────────────
    cudaMemcpyToSymbol(c_sigma_Li7, xs.materials[MAT_LI7].sigma_Li7_nnt.data(),
                       N_GROUPS * sizeof(float));

    // ── Majorant (pre-computed in loadDefaultCrossSections) ─────────────────
    cudaMemcpyToSymbol(c_sigma_major, xs.sigma_majorant.data(),
                       N_GROUPS * sizeof(float));
}

// ─── Initialize with built-in defaults (called once at startup) ──────────────
void initializeWithDefaults()
{
    ENDFLoader::CrossSectionSet xs;
    ENDFLoader::loadDefaultCrossSections(xs);
    uploadCrossSections(xs);
}

// ─── Initialize from a text file ─────────────────────────────────────────────
bool initializeFromFile(const char* path)
{
    ENDFLoader::CrossSectionSet xs;
    if (!ENDFLoader::loadFromTextFile(path, xs)) return false;
    uploadCrossSections(xs);
    return true;
}

} // namespace ENDFUpload
