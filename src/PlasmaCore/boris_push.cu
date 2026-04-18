//
// boris_push.cu
// Symplectic Boris pusher — conserves phase-space volume and adiabatic invariant.
// One CUDA thread per macro-particle.
//

#include "types.cuh"

// ─── Boris Push Kernel ─────────────────────────────────────────────────────────
//
//  Algorithm (non-relativistic):
//    1. Half-accelerate with E:  v^- = v_n + (q/m)(E)(dt/2)
//    2. Rotate with B:           v^+ = rotate(v^-, B, dt)
//    3. Half-accelerate again:   v_{n+1} = v^+ + (q/m)(E)(dt/2)
//    4. Advance position:        x_{n+1} = x_n + v_{n+1} * dt
//
//  Relativistic extension: replace v^- computation with gamma-factor; omitted
//  here but straightforward to add by maintaining u = gamma*v arrays.
//
__global__ void borisPush(
    float4* __restrict__ pos,           // [x,y,z,weight]         in/out
    float4* __restrict__ vel,           // [vx,vy,vz,species_id]  in/out
    const float4* __restrict__ E_field, // E sampled at particle position
    const float4* __restrict__ B_field, // B sampled at particle position
    float dt,
    int N_particles)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N_particles) return;

    // Load particle state (128-bit aligned reads)
    float4 p = pos[pid];
    float4 v = vel[pid];
    float4 E = E_field[pid];
    float4 B = B_field[pid];

    // Determine charge-to-mass ratio from species ID
    // species: 0 = electron, 1 = deuterium, 2 = tritium, 3 = alpha
    float qm; // q/m ratio [C/kg]
    int species = __float_as_int(v.w);  // packed species ID
    switch (species) {
        case 0:  qm = -(PC_E / PC_ME); break;  // electron
        case 1:  qm =  (PC_E / PC_MD); break;  // deuterium
        case 2:  qm =  (PC_E / PC_MT); break;  // tritium
        case 3:  qm =  (2.0f * PC_E / (4.0f * PC_MP)); break; // alpha (He-4)
        default: qm =  (PC_E / PC_MP); break;
    }

    float half_dt_qm = 0.5f * dt * qm;

    // ── Step 1: half electric impulse ────────────────────────────────────────
    float vx = v.x + half_dt_qm * E.x;
    float vy = v.y + half_dt_qm * E.y;
    float vz = v.z + half_dt_qm * E.z;

    // ── Step 2: magnetic rotation (Boris) ────────────────────────────────────
    // t = (q/m)(B)(dt/2),  s = 2t / (1 + |t|^2)
    float tx = half_dt_qm * B.x;
    float ty = half_dt_qm * B.y;
    float tz = half_dt_qm * B.z;
    float t2 = tx*tx + ty*ty + tz*tz;

    float sx = 2.0f * tx / (1.0f + t2);
    float sy = 2.0f * ty / (1.0f + t2);
    float sz = 2.0f * tz / (1.0f + t2);

    // v' = v^- + v^- x t
    float vpx = vx + (vy * tz - vz * ty);
    float vpy = vy + (vz * tx - vx * tz);
    float vpz = vz + (vx * ty - vy * tx);

    // v^+ = v^- + v' x s
    vx += (vpy * sz - vpz * sy);
    vy += (vpz * sx - vpx * sz);
    vz += (vpx * sy - vpy * sx);

    // ── Step 3: second half electric impulse ─────────────────────────────────
    vx += half_dt_qm * E.x;
    vy += half_dt_qm * E.y;
    vz += half_dt_qm * E.z;

    // ── Step 4: position advance ──────────────────────────────────────────────
    p.x += vx * dt;
    p.y += vy * dt;
    p.z += vz * dt;

    // Write back (coalesced stores)
    pos[pid] = p;
    vel[pid] = make_float4(vx, vy, vz, v.w);
}

// ─── Host Launch Wrapper ───────────────────────────────────────────────────────
void launchBorisPush(ParticleArrays& parts,
                     float4* E_at_particle,
                     float4* B_at_particle,
                     float dt,
                     cudaStream_t stream)
{
    constexpr int BLOCK = 256;
    int grid = (parts.N + BLOCK - 1) / BLOCK;
    borisPush<<<grid, BLOCK, 0, stream>>>(
        parts.pos, parts.vel,
        E_at_particle, B_at_particle,
        dt, parts.N);
}
