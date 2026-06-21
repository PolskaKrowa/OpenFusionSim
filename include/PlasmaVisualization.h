//
// include/PlasmaVisualization.h
// Spatial diagnostic data for the plasma visualization tabs.
//
//  The PlasmaCoreBridge's CUDA PIC kernel is abbreviated in this codebase
//  (see plasmacore.cu::runSimulation — the actual particle arrays and field
//  grids are allocated but not yet wired into the host update loop).  To
//  provide realistic-looking visualization data NOW, this class synthesizes
//  physics-consistent spatial diagnostics from the 0D ReactorState using
//  standard tokamak profile assumptions:
//
//    1. Temperature profile:  T(r) = T_0 · (1 - (r/a)²)  (parabolic, peaked)
//       This is the canonical L-mode/H-mode temperature profile shape and is
//       what real tokamaks measure via Thomson scattering.
//
//    2. Density profile:  n(r) = n_0 · (1 - (r/a)²)^0.5  (flatter than T)
//       H-mode edge pedestal modelled as a step at r = 0.95·a.
//
//    3. Particle distribution:  Maxwellian samples at the local T(r), with
//       the macro-particle weight proportional to n(r)·V_cell.  The birth
//       position is sampled uniformly within the flux surface at radius r.
//
//    4. Fusion rate profile:  R_fus(r) ∝ n_D(r) · n_T(r) · <σv>(T(r))
//       — peaks sharply at the magnetic axis (T and n both peak there).
//
//  When the full PIC is wired in, this class can be replaced with a thin
//  wrapper that copies device→host arrays.  The UI tabs call only this class,
//  so no UI code needs to change.
//
//  References for the profile shapes:
//    [1] Wesson, "Tokamaks" 4th ed., Ch. 4 — equilibrium profiles.
//    [2] ITER Physics Basis 1999, §2.5 — H-mode pedestal.
//    [3] Bosch-Hale 1992 for <σv>(T).
//

#pragma once
#include "ReactorState.h"
#include "confinement_physics.h"
#include <cmath>
#include <vector>
#include <array>
#include <cstdint>

// ─── Visualization grid + particle data ──────────────────────────────────────
class PlasmaViz {
public:
    // Grid resolution for 2D temperature / fusion-rate maps (R-Z plane)
    static constexpr int GRID_N_R = 48;   // radial bins (0 → a)
    static constexpr int GRID_N_Z = 32;   // vertical bins (-κa → +κa)
    // Particle samples per species (rendering only — the real PIC has ~10⁶)
    static constexpr int N_VIZ_PARTICLES = 4096;

    // Species IDs (must match boris_push.cu's speciesMassCharge)
    enum Species : int {
        SP_ELECTRON = 0,
        SP_DEUTERIUM = 1,
        SP_TRITIUM = 2,
        SP_ALPHA = 3,
        SP_PROTON = 4,
        SP_HE3 = 5,
        N_SPECIES = 6
    };

    struct Particle {
        float x, y, z;       // position [m] in vessel coordinates
        float vx, vy, vz;    // velocity [m/s]
        float weight;        // macro-particle weight (= n·V_macro / N_macro)
        int   species;
    };

    // ─── Constructor: pre-allocate buffers ────────────────────────────────────
    PlasmaViz()
    {
        temp_grid_.resize(GRID_N_R * GRID_N_Z, 0.0f);
        dens_grid_.resize(GRID_N_R * GRID_N_Z, 0.0f);
        fusion_rate_grid_.resize(GRID_N_R * GRID_N_Z, 0.0f);
        for (int s = 0; s < N_SPECIES; s++) {
            particles_[s].resize(N_VIZ_PARTICLES);
        }
        // Pre-seed a deterministic PRNG so particle positions are stable
        // frame-to-frame (only their velocities evolve with T).
        seed_ = 0x1234abcdULL;
    }

    // ─── Regenerate spatial data from current 0D state ────────────────────────
    //
    //  Call once per game tick (or every N ticks — 1 kHz is overkill).
    //  Updates the temperature/density/fusion-rate grids AND the per-species
    //  particle samples.
    //
    void update(const ReactorState& s)
    {
        // Tokamak geometry — read from ReactorState (set by the shape
        // controller) rather than hardcoded.  Previously these were all
        // hardcoded ITER values, so the visualization didn't update when
        // the operator changed the plasma shape via the κ/δ sliders.
        const float R = s.R_out_m;          // major radius [m]
        const float a = 2.0f;               // minor radius [m] (fixed)
        const float kappa = std::clamp(s.kappa, 1.0f, 2.5f);  // elongation
        const float a_eff = a * std::sqrt(kappa);  // effective minor radius

        // 0D plasma state from the bridge
        const float T0 = s.electron_temp_keV;     // on-axis temperature
        const float n0 = s.plasma_density_m3;     // on-axis density
        const float Ip = s.plasma_current_MA;
        const bool is_burning = (s.plasma_status == PlasmaStatus::Burning ||
                                  s.plasma_status == PlasmaStatus::Initiating);
        const bool is_cold    = (s.plasma_status == PlasmaStatus::Cold);

        // ── Temperature + density grids ────────────────────────────────────────
        //  T(r) = T0 · (1 - (r/a)²)^α_T    (α_T = 1.0 for L-mode, 1.5 for H-mode)
        //  n(r) = n0 · (1 - (r/a)²)^α_n    (α_n = 0.5 flat, 1.0 peaked)
        //  H-mode edge pedestal: T jumps to 0.3·T0 at r=0.95a (ITER-like)
        //
        const float alpha_T = 1.3f;   // mild peaking (H-mode)
        const float alpha_n = 0.5f;   // flat density
        const float T_pedestal_frac = 0.30f;  // edge pedestal is 30% of core T
        const float n_pedestal_frac = 0.60f;

        for (int j = 0; j < GRID_N_Z; j++) {
            float z = (2.0f * (float)j / (float)(GRID_N_Z - 1) - 1.0f) * kappa * a;
            for (int i = 0; i < GRID_N_R; i++) {
                float r = ((float)i + 0.5f) / (float)GRID_N_R * a;
                // Account for elongation: radial coordinate scales with κ
                float r_eff = std::sqrt(r * r + (z * z) / kappa);
                if (r_eff > a) r_eff = a;

                float x = r_eff / a;   // normalized minor radius (0 → 1)
                float shape = std::max(1.0f - x * x, 0.0f);

                // Temperature
                float T_local;
                if (is_cold) {
                    T_local = 0.0f;
                } else if (x > 0.95f) {
                    // Edge pedestal
                    T_local = T0 * T_pedestal_frac;
                } else {
                    T_local = T0 * std::pow(shape, alpha_T);
                }

                // Density
                float n_local;
                if (is_cold) {
                    n_local = 0.0f;
                } else if (x > 0.95f) {
                    n_local = n0 * n_pedestal_frac;
                } else {
                    n_local = n0 * std::pow(shape, alpha_n);
                }

                // Apply He-ash and impurity dilution to the D-T fuel fraction
                float f_DT = std::max(1.0f - s.helium_fraction - s.impurity_fraction, 0.0f);

                int idx = i + GRID_N_R * j;
                temp_grid_[idx] = T_local;
                dens_grid_[idx] = n_local;

                // ── Fusion rate density ∝ n_D · n_T · <σv>(T)  [W/m³] ──────────
                //  50:50 D-T mix on the fuel fraction
                float n_D = 0.5f * n_local * f_DT;
                float n_T = 0.5f * n_local * f_DT;
                float sv = ConfinementPhysics::boschHaleSigmaV_DT(T_local);  // m³/s
                constexpr float E_fus_J = 17.59e6f * 1.602176634e-19f;
                float p_fus_Wm3 = n_D * n_T * sv * E_fus_J;
                // Add a small steady glow even when "cold" so the colormap legend
                // is visible — but make it vanishingly small for realism.
                if (is_cold) p_fus_Wm3 = 0.0f;
                fusion_rate_grid_[idx] = p_fus_Wm3;
            }
        }

        // ── Particle distributions ────────────────────────────────────────────
        //  For each species, sample N_VIZ_PARTICLES positions within the
        //  plasma boundary (ellipse r²/a² + z²/(κa)² ≤ 1) with radial
        //  distribution weighted by n(r).  Velocity is sampled from a
        //  Maxwellian at the local T(r).
        //
        //  We use a deterministic PRNG (xorshift64) seeded per-species so
        //  that particles don't teleport between frames — instead their
        //  velocities evolve smoothly with T.
        //
        const float m_kg[N_SPECIES] = {
            9.1093837015e-31f,    // electron
            3.3435837724e-27f,    // D
            5.0073558862e-27f,    // T
            6.6446573357e-27f,    // α
            1.67262192369e-27f,   // p
            5.0082343773e-27f,    // ³He
        };

        // Thermal speed: v_th = sqrt(2 k_B T / m) — but T is in keV here.
        // k_B [J/K] × T_K = T_keV × 1.602e-16 J/keV, so v_th = sqrt(2 T[J] / m).
        // For T in keV: v_th = sqrt(2 × T_keV × 1.602e-16 / m_kg).

        for (int sp = 0; sp < N_SPECIES; sp++) {
            uint64_t s_local = seed_ ^ (0x9e3779b97f4a7c15ULL * (sp + 1));
            auto next = [&]() {
                // xorshift64
                s_local ^= s_local << 13;
                s_local ^= s_local >> 7;
                s_local ^= s_local << 17;
                return (uint32_t)(s_local >> 32);
            };
            auto uniform = [&]() -> float {
                return (float)(next() >> 8) / (float)(1 << 24);
            };
            auto gauss = [&]() -> float {
                // Box-Muller
                float u1 = std::max(uniform(), 1e-7f);
                float u2 = uniform();
                return std::sqrt(-2.0f * std::log(u1)) *
                       std::cos(2.0f * 3.14159265358979f * u2);
            };

            // Species-dependent abundance (relative number density)
            //  - D and T: each 50% of fuel fraction
            //  - electrons: equal to ion charge density (~n_e)
            //  - α: proportional to local fusion rate (3.5 MeV birth)
            //  - p, ³He: side-channels (small for pure D-T)
            float abundance;
            switch (sp) {
                case SP_ELECTRON:  abundance = 1.0f;  break;  // quasi-neutral
                case SP_DEUTERIUM: abundance = 0.5f * std::max(1.0f - s.helium_fraction - s.impurity_fraction, 0.0f); break;
                case SP_TRITIUM:   abundance = 0.5f * std::max(1.0f - s.helium_fraction - s.impurity_fraction, 0.0f); break;
                case SP_ALPHA:     abundance = s.helium_fraction * 0.5f; break;
                case SP_PROTON:    abundance = 0.0f;  break;  // D-D side-channel
                case SP_HE3:       abundance = 0.0f;  break;
                default:           abundance = 0.0f;
            }

            for (int p = 0; p < N_VIZ_PARTICLES; p++) {
                Particle& P = particles_[sp][p];

                // Sample position: reject outside the plasma ellipse
                float x_norm, y_norm, z_norm;
                int tries = 0;
                do {
                    x_norm = (uniform() * 2.0f - 1.0f);
                    y_norm = (uniform() * 2.0f - 1.0f);
                    z_norm = (uniform() * 2.0f - 1.0f);
                    tries++;
                } while (x_norm*x_norm + y_norm*y_norm + z_norm*z_norm > 1.0f && tries < 8);

                // Bias sampling toward the core using the n(r) profile:
                //  multiply the radius by a power-law shrink factor.
                float shrink = std::pow(uniform(), 1.0f / (1.0f + 2.0f * alpha_n));
                x_norm *= shrink;
                y_norm *= shrink;
                z_norm *= shrink;

                // Scale to vessel coordinates (R = 6.2 m, a = 2 m, κ = 1.7)
                P.x = R + x_norm * a;
                P.y = y_norm * a;
                P.z = z_norm * kappa * a;

                // Local temperature at this position
                float r_norm = std::sqrt(x_norm*x_norm + y_norm*y_norm + z_norm*z_norm);
                float shape_local = std::max(1.0f - r_norm * r_norm, 0.0f);
                float T_local = is_cold ? 0.0f :
                                (r_norm > 0.95f ? T0 * T_pedestal_frac
                                                : T0 * std::pow(shape_local, alpha_T));

                // Maxwellian velocity (in m/s): v_th = sqrt(2 T[J] / m)
                float T_J = T_local * 1.602176634e-16f;
                float v_th = (T_J > 0.0f) ? std::sqrt(2.0f * T_J / m_kg[sp]) : 0.0f;

                // For α particles, override with birth energy (3.5 MeV CoM-frame)
                if (sp == SP_ALPHA && is_burning) {
                    constexpr float E_alpha_J = 3.521e6f * 1.602176634e-19f;
                    v_th = std::sqrt(2.0f * E_alpha_J / m_kg[sp]);
                }

                P.vx = v_th * gauss();
                P.vy = v_th * gauss();
                P.vz = v_th * gauss();

                // Macro-particle weight = local density × abundance
                P.weight = abundance * n0 * std::pow(shape_local, alpha_n);
                P.species = sp;
            }
        }
    }

    // ─── Accessors for the UI ─────────────────────────────────────────────────
    const std::vector<float>& tempGrid() const { return temp_grid_; }
    const std::vector<float>& densGrid() const { return dens_grid_; }
    const std::vector<float>& fusionRateGrid() const { return fusion_rate_grid_; }
    const std::vector<Particle>& particles(Species sp) const { return particles_[(int)sp]; }

    int gridNR() const { return GRID_N_R; }
    int gridNZ() const { return GRID_N_Z; }

    // Tokamak geometry for axis labels
    float R_major() const { return 6.2f; }
    float a_minor() const { return 2.0f; }
    float kappa() const { return 1.7f; }

    // Per-species display info
    static const char* speciesName(Species s) {
        switch (s) {
            case SP_ELECTRON:  return "Electrons (e-)";
            case SP_DEUTERIUM: return "Deuterium (D+)";
            case SP_TRITIUM:   return "Tritium (T+)";
            case SP_ALPHA:     return "Alpha (He2+)";
            case SP_PROTON:    return "Protons (p+)";
            case SP_HE3:       return "Helium-3 (3He2+)";
            default:           return "?";
        }
    }
    // Species colour (RGBA) for the particle scatter plot
    static void speciesColor(Species s, float out[4]) {
        switch (s) {
            case SP_ELECTRON:  out[0]=0.30f; out[1]=0.55f; out[2]=0.95f; out[3]=1.0f; break; // blue
            case SP_DEUTERIUM: out[0]=0.20f; out[1]=0.85f; out[2]=0.85f; out[3]=1.0f; break; // cyan
            case SP_TRITIUM:   out[0]=0.95f; out[1]=0.65f; out[2]=0.10f; out[3]=1.0f; break; // amber
            case SP_ALPHA:     out[0]=0.95f; out[2]=0.18f; out[1]=0.18f; out[3]=1.0f; break; // red
            case SP_PROTON:    out[0]=0.85f; out[1]=0.85f; out[2]=0.20f; out[3]=1.0f; break; // yellow
            case SP_HE3:       out[0]=0.75f; out[1]=0.20f; out[2]=0.85f; out[3]=1.0f; break; // purple
            default:           out[0]=out[1]=out[2]=0.5f; out[3]=1.0f;
        }
    }

    // ─── Fusion rate time-series buffer ───────────────────────────────────────
    //  Updated by recordFusionHistory() each tick — used by the "Fusion reactions
    //  over time" tab.  Keeps the last N samples.
    void recordFusionHistory(float P_fus_MW, float P_alpha_MW, float P_neutron_MW,
                              float rate_reactions_per_s)
    {
        const int CAP = 1024;
        hist_pfus_.push_back(P_fus_MW);
        hist_palpha_.push_back(P_alpha_MW);
        hist_pneut_.push_back(P_neutron_MW);
        hist_rate_.push_back(rate_reactions_per_s);
        hist_time_.push_back(hist_time_.empty() ? 0.0f : hist_time_.back() + 1.0f);
        if ((int)hist_pfus_.size() > CAP) {
            hist_pfus_.erase(hist_pfus_.begin());
            hist_palpha_.erase(hist_palpha_.begin());
            hist_pneut_.erase(hist_pneut_.begin());
            hist_rate_.erase(hist_rate_.begin());
            hist_time_.erase(hist_time_.begin());
        }
    }
    const std::vector<float>& histFusionPower() const { return hist_pfus_; }
    const std::vector<float>& histAlphaPower()  const { return hist_palpha_; }
    const std::vector<float>& histNeutronPower()const { return hist_pneut_; }
    const std::vector<float>& histReactionRate()const { return hist_rate_; }
    const std::vector<float>& histTime()        const { return hist_time_; }

private:
    std::vector<float> temp_grid_;
    std::vector<float> dens_grid_;
    std::vector<float> fusion_rate_grid_;
    std::array<std::vector<Particle>, N_SPECIES> particles_;

    std::vector<float> hist_pfus_;
    std::vector<float> hist_palpha_;
    std::vector<float> hist_pneut_;
    std::vector<float> hist_rate_;
    std::vector<float> hist_time_;

    uint64_t seed_;
};
