#pragma once
//
// src/TritiumPlant/TritiumPlant.h
// Tritium extraction, isotope separation, and detritiation.
//
//  Tritium is bred in the blanket (Li-6 + n → T + He-4) and extracted via
//  the Tritium Extraction System (TES).  It then goes through an Isotope
//  Separation System (ISS) to separate T from H/D/He, and is stored in
//  uranium beds.  The detritiation system cleans tritiated air from the
//  building ventilation.
//
//  Key operator actions:
//    - Start/stop TES (extracts T from blanket → adds to fuel store)
//    - Start/stop detritiation (cleans air after a tritium release)
//    - Tritium accountancy: track how much T is in the plant vs. fuel store
//

#include "ReactorState.h"
#include "SimTime.h"

struct TritiumPlantConfig {
    // TES extraction rate: 1 g/h of T per kg/s of Li blanket flow
    float tes_rate_g_per_s_per_kg_s = 0.28e-3f;  // ~1 g/h at nominal flow
    float tes_max_g = 200.f;   // storage bed capacity

    // ISS processing rate
    float iss_rate_g_per_s = 0.01f;

    // Detritiation
    float detritiation_max_m3_s = 10.f;  // building ventilation rate
    float detritiation_efficiency = 0.999f;

    // Initial tritium in plant
    float initial_T_in_plant_g = 50.f;
};

class TritiumPlantSystem {
public:
    explicit TritiumPlantSystem(const TritiumPlantConfig& cfg = {});

    void update(ReactorState& state, const SimTime& t);

    // Cold-restart: T inventory back to initial, TES/detritiation off
    void reset();

    // Operator commands
    void startTES()  { tes_on_ = true; }
    void stopTES()   { tes_on_ = false; }
    void startDetritiation(float flow_m3_s = -1.f);  // -1 = max
    void stopDetritiation() { detritiation_on_ = false; detritiation_flow_m3_s_ = 0.f; }

    float tritiumInPlant() const { return T_in_plant_g_; }
    float tesRate() const { return tes_rate_g_s_; }
    bool  tesRunning() const { return tes_on_; }
    bool  detritiationRunning() const { return detritiation_on_; }

private:
    TritiumPlantConfig cfg_;

    float T_in_plant_g_;
    float T_in_storage_g_;   // long-term storage (uranium beds)
    bool  tes_on_            = false;
    float tes_rate_g_s_      = 0.f;
    bool  detritiation_on_   = false;
    float detritiation_flow_m3_s_ = 0.f;
};
