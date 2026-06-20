//
// src/PlasmaVizTabs.h
// Three new visualization tabs for the OpenFusionSim frontend:
//   1. Temperature distribution (R-Z cross-section, heatmap + isotherms)
//   2. Fusion reactions over time (multi-line plot of P_fus, P_α, P_neutron, rate)
//   3. Particle distribution (scatter plot with per-species toggles)
//
//  These are kept in a separate translation unit so main.cpp doesn't grow
//  beyond what's manageable.  main.cpp #includes this header and dispatches
//  to the render functions from its tab switch.
//

#pragma once
#include "ReactorState.h"
#include "SimTime.h"
#include "PlasmaVisualization.h"

// ─── Tab 8: Temperature distribution ─────────────────────────────────────────
void TabTemperature(const ReactorState& s, const PlasmaViz& viz);

// ─── Tab 9: Fusion reactions over time ───────────────────────────────────────
void TabFusionReactions(const ReactorState& s, const PlasmaViz& viz);

// ─── Tab 10: Particle distribution ───────────────────────────────────────────
//  `visible_species_mask` is a bitfield (bit n = species n visible).
//  Modified by the user via checkboxes in the tab itself.
void TabParticleDistribution(const ReactorState& s, const PlasmaViz& viz,
                              uint32_t& visible_species_mask,
                              int& projection_mode);  // 0=R-Z, 1=R-φ (top-down)
