#pragma once
//
// src/TurbineSystem/TurbineSystem.h
// Manages all 4 TurbineUnitControllers; writes aggregates to ReactorState.
//

#include "TurbineUnit.h"
#include "ReactorState.h"
#include "SimTime.h"
#include <array>

class TurbineSystem {
public:
    TurbineSystem();

    void update(ReactorState& state, const SimTime& t,
                float grid_frequency_Hz,
                const float salt_heat_MW[4]); // heat from MoltenSalt per SG

    // Access individual units for UI and MoltenSalt heat assignment
    TurbineUnitController& unit(int i) { return units_[i]; }
    const TurbineUnitController& unit(int i) const { return units_[i]; }

    float totalPowerMW() const;

private:
    std::array<TurbineUnitController, 4> units_;
};