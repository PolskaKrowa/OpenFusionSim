#pragma once

//
// src/ThermalHydraulics/ThermalHydraulics.h
// Primary coolant loop, breeding blanket, and first-wall thermal model.
//
//  The blanket has two jobs:
//    1. Absorb 80 % of the fusion energy (the 14.1 MeV neutrons).
//    2. Breed tritium via Li-6(n,t)He-4 and Li-7(n,n')t reactions.
//
//  Coolant: Water (PWR-style) or FLiBe (molten salt) — selectable.
//  Model: lumped-parameter 1D loop (inlet → blanket → outlet → HX → back).
//

#include "ReactorState.h"
#include "SimTime.h"

enum class CoolantType { Water, FLiBe };

struct ThermalHydraulicsConfig {
    CoolantType coolant      = CoolantType::Water;

    // Blanket geometry
    float blanket_volume_m3  = 1500.0f;   // total blanket volume [m^3]
    float first_wall_area_m2 =  700.0f;   // [m^2]
    float Li6_enrichment     =    0.60f;  // 60 % Li-6 enrichment (natural: 7.5 %)

    // Coolant properties (water at ~300 °C, 15 MPa)
    float Cp_J_kgK           = 5200.0f;  // specific heat [J/(kg·K)]
    float rho_kg_m3          =  720.0f;  // density [kg/m^3]
    float max_flow_kg_s       = 18000.0f; // nominal flow rate [kg/s]

    // Heat exchanger (blanket → steam generator)
    float HX_UA_W_K           = 200e6f;  // overall heat transfer coeff × area [W/K]

    // Thermal limits
    float max_coolant_temp_K  = 900.0f;  // FLiBe limit; water ~620 K
    float max_first_wall_K    = 1800.0f; // tungsten armour limit [K]

    // Neutronics (simplified)
    float blanket_energy_mult = 1.20f;   // neutron energy amplification in Li blanket
    float TBR_target          = 1.05f;   // target tritium breeding ratio
};

class ThermalHydraulics {
public:
    explicit ThermalHydraulics(const ThermalHydraulicsConfig& cfg);

    void update(ReactorState& state, const SimTime& t);

    // Cold-restart: drop blanket/first-wall/coolant temperatures back to
    // their construction defaults.  Without this, a SCRAM fired while the
    // blanket is hot leaves blanket_temp_K_ high across the reset, and the
    // very first update() after RESET writes that hot temperature back into
    // state.first_wall_temp_K / state.thermal_runaway — re-tripping SCRAM.
    void reset();

    float outletTemp()    const { return coolant_outlet_K_; }
    float thermalPowerMW()const { return thermal_power_MW_; }

private:
    void  updateBlanket     (ReactorState& state, float dt);
    void  updateCoolantLoop (ReactorState& state, float dt);
    void  updateFirstWall   (ReactorState& state, float dt);
    float computeTBR        (float neutron_flux) const;

    ThermalHydraulicsConfig cfg_;

    float blanket_temp_K_    = 573.0f;
    float coolant_inlet_K_   = 573.0f;
    float coolant_outlet_K_  = 773.0f;
    float first_wall_temp_K_ = 300.0f;
    float thermal_power_MW_  = 0.0f;
    float flow_kg_s_         = 0.0f;
    float tbr_current_       = 0.0f;

    // ── Neutron damage accountancy (campaign-lifetime, cleared on reset) ────
    //  Integrated neutron wall load [MW·a/m²] and the resulting displacement
    //  damage in the first-wall steel [dpa].  ~10 dpa per MW·a/m² for RAFM
    //  steel under a 14 MeV DT spectrum.
    float fw_fluence_MWa_m2_ = 0.0f;
    float fw_dpa_            = 0.0f;
};