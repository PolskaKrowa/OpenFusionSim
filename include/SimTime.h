#pragma once

//
// include/Units.h
// Compile-time unit conversion helpers and physical constants (host-side).
// GPU-side constants live in plasmacore/types.cuh.
//

namespace units {

// ── Physical constants ────────────────────────────────────────────────────────
constexpr double C        = 2.99792458e8;    // speed of light         [m/s]
constexpr double E_CHARGE = 1.60217663e-19;  // elementary charge      [C]
constexpr double K_B      = 1.38064852e-23;  // Boltzmann constant     [J/K]
constexpr double M_E      = 9.10938370e-31;  // electron mass          [kg]
constexpr double M_P      = 1.67262192e-27;  // proton mass            [kg]
constexpr double EPS0     = 8.85418782e-12;  // vacuum permittivity    [F/m]
constexpr double MU0      = 1.25663706e-6;   // vacuum permeability    [H/m]

// ── Energy conversions ────────────────────────────────────────────────────────
constexpr double EV_TO_J  = E_CHARGE;
constexpr double J_TO_EV  = 1.0 / E_CHARGE;
constexpr double KEV_TO_J = E_CHARGE * 1e3;
constexpr double J_TO_KEV = 1.0 / KEV_TO_J;
constexpr double MEV_TO_J = E_CHARGE * 1e6;
constexpr double J_TO_MEV = 1.0 / MEV_TO_J;

// ── Power conversions ─────────────────────────────────────────────────────────
constexpr double MW_TO_W  = 1e6;
constexpr double W_TO_MW  = 1e-6;
constexpr double GW_TO_W  = 1e9;

// ── Inline conversion functions ───────────────────────────────────────────────
constexpr double keV_to_K(double keV)  { return keV * KEV_TO_J / K_B; }
constexpr double K_to_keV(double K)    { return K * K_B / KEV_TO_J; }
constexpr double MW_to_W (double MW)   { return MW * 1e6; }
constexpr double W_to_MW (double W)    { return W  * 1e-6; }

} // namespace units


#pragma once

//
// include/SimTime.h
// Lightweight simulation time manager.
// Each module's update() receives a SimTime reference.
//

struct SimTime {
    double  total_s  = 0.0;    // elapsed simulation time [s]
    float   dt_s     = 1e-3f;  // current timestep [s]
    int     tick     = 0;      // integer tick counter (wraps at INT_MAX)

    void advance() {
        total_s += dt_s;
        tick++;
    }

    // Convenience: true every N ticks
    bool every(int N) const { return (tick % N) == 0; }

    // True once, at the first tick
    bool is_first() const { return tick == 0; }
};