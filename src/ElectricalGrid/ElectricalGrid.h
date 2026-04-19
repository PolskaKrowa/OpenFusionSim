#pragma once
//
// src/ElectricalGrid/ElectricalGrid.h
// Electrical grid simulation: frequency, generator synchronisation, on-site loads.
//
//  Grid model:
//    - Infinite external grid (grid holds frequency close to 50 Hz)
//    - Each generator can be synchronised and connected via a breaker
//    - Local grid frequency is a weighted mix of generator inertia and external grid
//    - On-site loads always powered first; surplus exported (or deficit imported)
//

#include "ReactorState.h"
#include "SimTime.h"
#include <array>

// ─── On-site load group ───────────────────────────────────────────────────────
struct SiteLoad {
    const char* name;
    float load_MW;
    bool  energised;
    bool  essential; // cannot be shed
};

// ─── Per-generator bus state (read from TurbineSystem, written here) ──────────
struct GeneratorBus {
    float frequency_Hz   = 0.f;
    float phase_rad      = 0.f;
    float power_MW       = 0.f;
    bool  breaker_closed = false;
    bool  synch_ok       = false;  // frequency and phase within limits
};

// ─── Grid state ───────────────────────────────────────────────────────────────
struct GridState {
    float frequency_Hz           = 50.f;
    float frequency_deviation_Hz = 0.f;   // Δf from nominal
    float voltage_pu             = 1.00f; // per-unit voltage
    float total_generation_MW    = 0.f;
    float total_site_load_MW     = 0.f;
    float export_MW              = 0.f;   // positive = to external grid
    float import_MW              = 0.f;   // positive = from external grid
    bool  underfrequency_alarm   = false;
    bool  overfrequency_alarm    = false;
    bool  underfrequency_trip    = false; // <47.5 Hz → station blackout
    bool  loss_of_offsite_power  = false; // external grid unavailable
    std::array<GeneratorBus, 4> bus;
};

class ElectricalGridSystem {
public:
    ElectricalGridSystem();

    void update(ReactorState& state, const SimTime& t,
                const struct TurbineSystem& turbines);

    const GridState& grid() const { return grid_; }
    GridState&       grid()       { return grid_; }

    // Toggle individual site loads
    void setLoadEnergised(int idx, bool on);
    int  numLoads() const { return (int)site_loads_.size(); }
    const SiteLoad& siteLoad(int i) const { return site_loads_[i]; }

    // Check if generator i is ready to sync (call before closing breaker)
    bool canSync(int gen_idx) const;

    // Force loss-of-offsite-power event
    void triggerLOOP() { grid_.loss_of_offsite_power = true; }
    void restoreOffsite() { grid_.loss_of_offsite_power = false; }

private:
    void updateFrequency(float dt, float total_gen_MW);
    void updateSynchCheck(const struct TurbineSystem& turbines);
    void updateSiteLoads();
    void updateExportImport();

    GridState grid_;
    std::array<SiteLoad, 12> site_loads_;
    float external_grid_inertia_ = 50000.f; // MW·s (infinite bus approximation)
};