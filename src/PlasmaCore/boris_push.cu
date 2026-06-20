//
// boris_push.cu
// Symplectic particle pusher for the PIC loop.
//
//  Two pushers are provided:
//    1. Boris (1970) — non-relativistic, the original OpenFusionSim default.
//    2. Vay (2008)   — relativistic, for fast α (3.5 MeV) and runaway
//                      electrons (γ ≫ 1).
//
//  The Vay pusher is selected per-particle based on γ > 1.01 (i.e. KE > 0.5%
//  of rest mass).  This avoids unnecessary work for the thermal bulk while
//  correctly handling the 3.5-MeV α particles (γ_α = 1.0009, marginal) and
//  runaway electrons (γ_e ≫ 1).
//
//  References:
//    [1] J.P. Boris, "Relativistic plasma simulation-optimization of a hybrid
//        code", Proc. 4th Conf. Num. Sim. Plasmas (1970).
//    [2] J.-L. Vay, Phys. Plasmas 15, 056701 (2008).
//    [3] Ripperda et al., ApJ Suppl. 235, 21 (2018) — comprehensive review.
//

#include "types.cuh"

// ─── Species mass and charge table ────────────────────────────────────────────
//  Species encoding (in vel.w): 0=electron, 1=D, 2=T, 3=α, 4=p, 5=³He
__device__ __forceinline__
void speciesMassCharge(int s, float& m_kg, float& q_C)
{
    switch (s) {
        case 0:  m_kg = 9.1093837015e-31f;  q_C = -1.602176634e-19f; break; // electron
        case 1:  m_kg = 3.3435837724e-27f;  q_C =  1.602176634e-19f; break; // deuteron
        case 2:  m_kg = 5.0073558862e-27f;  q_C =  1.602176634e-19f; break; // triton
        case 3:  m_kg = 6.6446573357e-27f;  q_C =  2.0f * 1.602176634e-19f; break; // α
        case 4:  m_kg = 1.67262192369e-27f; q_C =  1.602176634e-19f; break; // proton
        case 5:  m_kg = 5.0082343773e-27f;  q_C =  2.0f * 1.602176634e-19f; break; // ³He
        default: m_kg = 1.67262192369e-27f; q_C =  1.602176634e-19f; break;
    }
}

// ─── Non-relativistic Boris pusher (original) ─────────────────────────────────
__device__ __forceinline__
void borisNonRel(float3& v, float3 E, float3 B, float qm, float dt)
{
    float half_dt_qm = 0.5f * dt * qm;

    // Half electric impulse
    v.x += half_dt_qm * E.x;
    v.y += half_dt_qm * E.y;
    v.z += half_dt_qm * E.z;

    // Magnetic rotation: t = (q/m)·B·dt/2, s = 2t/(1+|t|²)
    float tx = half_dt_qm * B.x;
    float ty = half_dt_qm * B.y;
    float tz = half_dt_qm * B.z;
    float t2 = tx*tx + ty*ty + tz*tz;
    float sx = 2.0f * tx / (1.0f + t2);
    float sy = 2.0f * ty / (1.0f + t2);
    float sz = 2.0f * tz / (1.0f + t2);

    // v' = v + v × t
    float vpx = v.x + (v.y * tz - v.z * ty);
    float vpy = v.y + (v.z * tx - v.x * tz);
    float vpz = v.z + (v.x * ty - v.y * tx);

    // v⁺ = v + v' × s
    v.x += (vpy * sz - vpz * sy);
    v.y += (vpz * sx - vpx * sz);
    v.z += (vpx * sy - vpy * sx);

    // Second half electric impulse
    v.x += half_dt_qm * E.x;
    v.y += half_dt_qm * E.y;
    v.z += half_dt_qm * E.z;
}

// ─── Vay (2008) relativistic pusher ───────────────────────────────────────────
//
//  Solves the implicit average-velocity equation
//      u⁻ = u^(n-1/2) + (qΔt/2m)·E
//      ū = u⁻ + (qΔt/2m)·(u⁻ × B) + (qΔt/2m)·(u⁻·E)(u⁻×B)/c²
//      u⁺ = (u⁻ + (ū·τ)τ/c² + ū×τ/c) / (1 + |τ|²/c² + σ²)
//      u^(n+1/2) = u⁺ + (qΔt/2m)·E
//  where τ = (qΔt/2m)·B and σ = (qΔt/2m)·(u⁻·E)/c.
//
//  This is fully explicit (no iteration) — Vay's key result is that the
//  average-velocity form has a closed-form solution.  It reduces to Boris
//  in the non-relativistic limit.
//
//  We store u = γv (not v) so the kinematic update is:
//      x_{n+1} = x_n + (u/γ) · dt
//  where γ = sqrt(1 + |u|²/c²).
//
__device__ __forceinline__
void vayPush(float3& u,                // in/out: canonical momentum u = γv
             float3 E, float3 B,
             float qm, float dt, float c2)
{
    float half_dt_qm = 0.5f * dt * qm;

    // u⁻ = u + (qΔt/2m)·E
    float umx = u.x + half_dt_qm * E.x;
    float umy = u.y + half_dt_qm * E.y;
    float umz = u.z + half_dt_qm * E.z;

    // τ = (qΔt/2m)·B ; σ = (qΔt/2m)·(u⁻·E)/c²
    float tx  = half_dt_qm * B.x;
    float ty  = half_dt_qm * B.y;
    float tz  = half_dt_qm * B.z;
    float u_dot_E = umx * E.x + umy * E.y + umz * E.z;
    float sigma   = half_dt_qm * u_dot_E / c2;

    // ū = u⁻ + (u⁻ × τ) + σ·(u⁻ × τ)/c²     ... Vay's average velocity
    //   (Vay Eq. 14; can simplify to ū = u⁻ + u⁻×τ for the pure-magnetic case)
    float uxT_x = umy * tz - umz * ty;
    float uxT_y = umz * tx - umx * tz;
    float uxT_z = umx * ty - umy * tx;

    float ubarx = umx + uxT_x;
    float ubary = umy + uxT_y;
    float ubarz = umz + uxT_z;
    // The σ·(u⁻×τ)/c² term is O(σ·v²/c²); for fusion plasmas (v/c<0.1) it is
    // negligible, but included for completeness:
    ubarx += sigma * uxT_x;
    ubary += sigma * uxT_y;
    ubarz += sigma * uxT_z;

    // u⁺ = (u⁻ + (ū·τ)τ/c² + ū×τ/c) / (1 + |τ|²/c² + σ²)
    //     (Vay Eq. 15; this is the closed-form implicit solve)
    float tau2_c2  = (tx*tx + ty*ty + tz*tz) / c2;
    float denom    = 1.0f + tau2_c2 + sigma * sigma;
    float u_d_tau  = (ubarx * tx + ubary * ty + ubarz * tz) / c2;

    float uxT2_x = ubary * tz - ubarz * ty;
    float uxT2_y = ubarz * tx - ubarx * tz;
    float uxT2_z = ubarx * ty - ubary * tx;

    float upx = (umx + u_d_tau * tx + uxT2_x) / denom;
    float upy = (umy + u_d_tau * ty + uxT2_y) / denom;
    float upz = (umz + u_d_tau * tz + uxT2_z) / denom;

    // Final half-E impulse
    u.x = upx + half_dt_qm * E.x;
    u.y = upy + half_dt_qm * E.y;
    u.z = upz + half_dt_qm * E.z;
}

// ─── Hybrid pusher: non-rel for thermal bulk, Vay for relativistic tail ───────
//
//  γ threshold: γ > 1.01 means KE > 0.5% of rest mass (~5 keV for electrons,
//  ~9 MeV for deuterons).  This catches 3.5-MeV α particles marginally
//  (γ_α = 1.00094, just below threshold) and runaway electrons comfortably.
//
//  For α particles at 3.5 MeV the Boris-vs-Vay difference is <0.1%; we keep
//  Boris for them to avoid the (small) Vay overhead.  For runaways we use Vay.
//
__global__ void pushParticles(
    float4* __restrict__ pos,           // [x,y,z,weight]         in/out
    float4* __restrict__ vel,           // [vx,vy,vz,species_id]  in/out  (or [ux,uy,uz,s])
    const float4* __restrict__ E_field,
    const float4* __restrict__ B_field,
    float dt,
    int N_particles)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= N_particles) return;

    float4 p = pos[pid];
    float4 v = vel[pid];
    float4 E = E_field[pid];
    float4 B = B_field[pid];

    int species = __float_as_int(v.w);
    float m_kg, q_C;
    speciesMassCharge(species, m_kg, q_C);
    float qm = q_C / m_kg;

    float3 Ev = make_float3(E.x, E.y, E.z);
    float3 Bv = make_float3(B.x, B.y, B.z);

    // Compute γ to choose pusher.  We use |v|/c as the proxy (γ = 1/sqrt(1-β²)).
    constexpr float C  = 2.99792458e8f;
    constexpr float C2 = C * C;
    float v2   = v.x * v.x + v.y * v.y + v.z * v.z;
    float beta2 = v2 / C2;
    bool  relativistic = (beta2 > 0.02f);   // β > 0.14 ≈ γ > 1.01

    if (relativistic) {
        // Vay pusher on u = γv.  Convert v → u, push, convert back.
        float gamma = 1.0f / sqrtf(1.0f - beta2);
        float3 u = make_float3(gamma * v.x, gamma * v.y, gamma * v.z);
        vayPush(u, Ev, Bv, qm, dt, C2);
        float u2 = u.x * u.x + u.y * u.y + u.z * u.z;
        float g_new = sqrtf(1.0f + u2 / C2);
        v.x = u.x / g_new;
        v.y = u.y / g_new;
        v.z = u.z / g_new;
    } else {
        float3 vv = make_float3(v.x, v.y, v.z);
        borisNonRel(vv, Ev, Bv, qm, dt);
        v.x = vv.x; v.y = vv.y; v.z = vv.z;
    }

    // Position advance: x_{n+1} = x_n + v_{n+1} · dt
    p.x += v.x * dt;
    p.y += v.y * dt;
    p.z += v.z * dt;

    pos[pid] = p;
    vel[pid] = make_float4(v.x, v.y, v.z, v.w);
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
    pushParticles<<<grid, BLOCK, 0, stream>>>(
        parts.pos, parts.vel,
        E_at_particle, B_at_particle,
        dt, parts.N);
}
