//
// plasmacore.cu
// Top-level simulation loop for FusionSim PlasmaCore.
//
//  Each timestep:
//    0. [Every N_sort steps] Sort particles by cell, compute CSR offsets.
//    1. Deposit charge and current to grid.
//    2. Solve field equations (Poisson spectral or FDTD Maxwell).
//    3. Gather fields from grid to particle positions.
//    4. Push particles (Boris).
//    5. Sample Coulomb collisions.
//    6. Sample fusion reactions; inject new alphas; hand neutrons to transport.
//    7. Transport neutrons (accumulate TBR + heat deposition).
//    8. Diagnostics (energy, momentum, particle counts).
//

#include "types.cuh"

// Forward declarations of host wrappers defined in other .cu files
void launchBorisPush(ParticleArrays&, float4*, float4*, float, cudaStream_t);
void launchGatherFields(float4*, float4*, const float*, const float*,
                        const ParticleArrays&, GridParams, cudaStream_t);
void launchDepositCharge(float*, float*, const ParticleArrays&, GridParams,
                         float, cudaStream_t);
void launchFDTD(float*, float*, float*, const GridParams&, float, cudaStream_t);
void launchCoulombCollisions(ParticleArrays&, const float*, const int*, float,
                             float, curandState*, GridParams, int, cudaStream_t);
void launchFusionReactions(const float4*, const float4*, const float4*, const float4*,
                           const int*, const int*, ReactionProduct*, curandState*,
                           float, GridParams, int, cudaStream_t);
void launchNeutronTransport(NeutronParticle*, const float*, TritiumProductionMap*,
                            HeatDepositionMap*, curandState*, GridParams,
                            int, float, int, cudaStream_t);
void sortParticlesByCell(ParticleArrays&, int*, struct SortContext&, GridParams, cudaStream_t);

// ─── Simulation Configuration ─────────────────────────────────────────────────
struct SimConfig {
    int    N_particles   = 1'000'000;
    int    N_steps       = 100'000;
    int    sort_interval = 20;       // re-sort every N steps
    float  dt            = 1e-12f;   // timestep [s]  (plasma/cyclotron frequency limited)
    float  coulomb_log   = 15.0f;    // ln(Λ) for DT plasma at ~10 keV
    float  E_n_cutoff    = 1e-7f;    // neutron energy cutoff [MeV]
    int    max_n_coll    = 1000;     // max neutron history collisions
    GridParams grid;
};

// ─── Device buffer collection ─────────────────────────────────────────────────
struct SimBuffers {
    // Particle arrays (all species in single SoA; species encoded in vel.w)
    ParticleArrays ions;        // D + T + alpha
    ParticleArrays electrons;

    // Per-species D/T subsets (views into ions, or separate arrays)
    float4* pos_D; float4* vel_D;
    float4* pos_T; float4* vel_T;

    // Field grids [3 * Nx*Ny*Nz floats]
    float* E_grid;
    float* B_grid;
    float* rho_grid;  // [Nx*Ny*Nz]
    float* J_grid;    // [3 * Nx*Ny*Nz]

    // Field at particle positions (temporary, per-step)
    float4* E_at_ion;   float4* B_at_ion;
    float4* E_at_elec;  float4* B_at_elec;

    // Particle weights
    float* ion_weights;
    float* elec_weights;

    // CSR offsets for sorted ions and electrons
    int* ion_cell_start;
    int* elec_cell_start;
    int* cell_D_start;
    int* cell_T_start;

    // Reaction products buffer
    ReactionProduct* reaction_products;

    // Neutron arrays
    NeutronParticle* neutrons;
    int N_neutrons_max;

    // TBR + heat maps
    TritiumProductionMap* tbr;
    HeatDepositionMap*    q_dot;

    // Material map for neutron transport
    float* material_map;

    // RNG states
    curandState* rng_ion;
    curandState* rng_elec;
    curandState* rng_neutron;
    curandState* rng_cell;    // one state per cell for fusion/collision kernels
};

// ─── Allocation helper (sizes inferred from SimConfig) ────────────────────────
void allocSimBuffers(SimBuffers& buf, const SimConfig& cfg)
{
    const GridParams& g = cfg.grid;
    int N_cells = g.Nx * g.Ny * g.Nz;
    int N = cfg.N_particles;

    // Ion particle SoA
    cudaMalloc(&buf.ions.pos, N * sizeof(float4));
    cudaMalloc(&buf.ions.vel, N * sizeof(float4));
    cudaMalloc(&buf.ions.cell, N * sizeof(int));
    buf.ions.N = N;

    // Electron SoA (same size as ions for quasi-neutrality)
    cudaMalloc(&buf.electrons.pos, N * sizeof(float4));
    cudaMalloc(&buf.electrons.vel, N * sizeof(float4));
    cudaMalloc(&buf.electrons.cell, N * sizeof(int));
    buf.electrons.N = N;

    // Field grids
    cudaMalloc(&buf.E_grid,   3 * N_cells * sizeof(float));
    cudaMalloc(&buf.B_grid,   3 * N_cells * sizeof(float));
    cudaMalloc(&buf.rho_grid, N_cells * sizeof(float));
    cudaMalloc(&buf.J_grid,   3 * N_cells * sizeof(float));

    // Fields at particle positions
    cudaMalloc(&buf.E_at_ion,  N * sizeof(float4));
    cudaMalloc(&buf.B_at_ion,  N * sizeof(float4));
    cudaMalloc(&buf.E_at_elec, N * sizeof(float4));
    cudaMalloc(&buf.B_at_elec, N * sizeof(float4));

    // Weights
    cudaMalloc(&buf.ion_weights,  N * sizeof(float));
    cudaMalloc(&buf.elec_weights, N * sizeof(float));

    // CSR offsets (+1 for sentinel)
    cudaMalloc(&buf.ion_cell_start,  (N_cells + 1) * sizeof(int));
    cudaMalloc(&buf.elec_cell_start, (N_cells + 1) * sizeof(int));
    cudaMalloc(&buf.cell_D_start,    (N_cells + 1) * sizeof(int));
    cudaMalloc(&buf.cell_T_start,    (N_cells + 1) * sizeof(int));

    // Reaction products: one slot per D particle (worst case all react)
    cudaMalloc(&buf.reaction_products, N * sizeof(ReactionProduct));

    // Neutrons (born from fusion; size TBD; conservatively same as ions)
    buf.N_neutrons_max = N;
    cudaMalloc(&buf.neutrons, N * sizeof(NeutronParticle));

    // Heat and TBR maps
    buf.tbr   = new TritiumProductionMap{nullptr, g.Nx, g.Ny, g.Nz};
    buf.q_dot = new HeatDepositionMap   {nullptr, g.Nx, g.Ny, g.Nz};
    cudaMalloc(&buf.tbr->tbr_voxel,   N_cells * sizeof(float));
    cudaMalloc(&buf.q_dot->q_dot,     N_cells * sizeof(float));
    cudaMemset(buf.tbr->tbr_voxel, 0, N_cells * sizeof(float));
    cudaMemset(buf.q_dot->q_dot,   0, N_cells * sizeof(float));

    // Material map
    cudaMalloc(&buf.material_map, N_cells * sizeof(float));

    // RNG states
    cudaMalloc(&buf.rng_ion,     N * sizeof(curandState));
    cudaMalloc(&buf.rng_elec,    N * sizeof(curandState));
    cudaMalloc(&buf.rng_neutron, N * sizeof(curandState));
    cudaMalloc(&buf.rng_cell,    N_cells * sizeof(curandState));
}

// ─── Main Simulation Loop ─────────────────────────────────────────────────────
void runSimulation(SimBuffers& buf, SimConfig& cfg, struct SortContext& sortCtx)
{
    cudaStream_t stream_main, stream_neutron;
    cudaStreamCreate(&stream_main);
    cudaStreamCreate(&stream_neutron);

    int N_cells = cfg.grid.Nx * cfg.grid.Ny * cfg.grid.Nz;

    for (int step = 0; step < cfg.N_steps; step++) {

        // ── 0. Sort particles by cell (every sort_interval steps) ────────────
        if (step % cfg.sort_interval == 0) {
            sortParticlesByCell(buf.ions,      buf.ion_cell_start,
                                sortCtx, cfg.grid, stream_main);
            sortParticlesByCell(buf.electrons, buf.elec_cell_start,
                                sortCtx, cfg.grid, stream_main);
            // TODO: also recompute cell_D_start and cell_T_start here
            // by filtering sorted ion array by species — left as exercise.
        }

        // ── 1. Deposit rho and J from particles to grid ───────────────────────
        launchDepositCharge(buf.rho_grid, buf.J_grid,
                            buf.ions, cfg.grid, cfg.dt, stream_main);
        // Electron contribution (opposite charge, subtract)
        // In practice: run depositCharge with sign=-1 accumulation or merge.

        // ── 2. Solve field equations ──────────────────────────────────────────
        // Option B (full EM) used here; swap for spectralPoissonSolve if electrostatic.
        launchFDTD(buf.E_grid, buf.B_grid, buf.J_grid,
                   cfg.grid, cfg.dt, stream_main);

        // ── 3. Gather fields at particle positions ────────────────────────────
        launchGatherFields(buf.E_at_ion,  buf.B_at_ion,
                           buf.E_grid, buf.B_grid,
                           buf.ions, cfg.grid, stream_main);
        launchGatherFields(buf.E_at_elec, buf.B_at_elec,
                           buf.E_grid, buf.B_grid,
                           buf.electrons, cfg.grid, stream_main);

        // ── 4. Push all particles with Boris integrator ───────────────────────
        launchBorisPush(buf.ions,      buf.E_at_ion,  buf.B_at_ion,
                        cfg.dt, stream_main);
        launchBorisPush(buf.electrons, buf.E_at_elec, buf.B_at_elec,
                        cfg.dt, stream_main);

        // ── 5. Coulomb collisions ─────────────────────────────────────────────
        launchCoulombCollisions(buf.ions, buf.ion_weights, buf.ion_cell_start,
                                cfg.dt, cfg.coulomb_log, buf.rng_cell,
                                cfg.grid, N_cells, stream_main);

        // ── 6. Fusion reaction sampling ───────────────────────────────────────
        launchFusionReactions(buf.pos_D, buf.vel_D,
                              buf.pos_T, buf.vel_T,
                              buf.cell_D_start, buf.cell_T_start,
                              buf.reaction_products,
                              buf.rng_cell,
                              cfg.dt, cfg.grid, N_cells, stream_main);

        // TODO: harvest reaction_products — inject alphas into ion array,
        //       pack neutrons into buf.neutrons[] on host or via a compact kernel.
        //       N_neutrons_this_step computed via cub::DeviceReduce::Sum on active flags.

        // ── 7. Neutron transport (can overlap ion push on separate stream) ─────
        // Assumes buf.neutrons is populated and N_neutrons_this_step is known.
        int N_neutrons_this_step = buf.N_neutrons_max; // placeholder
        launchNeutronTransport(buf.neutrons, buf.material_map,
                               buf.tbr, buf.q_dot, buf.rng_neutron,
                               cfg.grid, N_neutrons_this_step,
                               cfg.E_n_cutoff, cfg.max_n_coll,
                               stream_neutron);

        // ── 8. Diagnostics (every 100 steps) ─────────────────────────────────
        if (step % 100 == 0) {
            cudaStreamSynchronize(stream_main);
            // TODO: launch energy/momentum diagnostic kernels and copy to host
            // printf("Step %d / %d\n", step, cfg.N_steps);
        }
    }

    cudaStreamSynchronize(stream_main);
    cudaStreamSynchronize(stream_neutron);
    cudaStreamDestroy(stream_main);
    cudaStreamDestroy(stream_neutron);
}
