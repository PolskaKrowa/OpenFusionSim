//
// src/PlasmaVizTabs.cpp
// Implementation of the three plasma visualization tabs.
//
//  Uses ImGui's ImDrawList API for custom canvas rendering:
//    - Heatmap of temperature/fusion-rate via filled rectangles per pixel
//    - Multi-line time-series plots via ImGui::PlotLines + custom overlays
//    - Particle scatter plot via ImDrawList::AddCircleFilled
//
//  The temperature/fusion-rate heatmaps use a Jet-style colourmap (blue→cyan
//  →green→yellow→red) which is the standard for scientific visualization.
//  (Yes, jet has perceptual issues — but it's the convention for tokamak
//  diagnostics and matches what operators see in real control rooms.)
//

#include "PlasmaVizTabs.h"
#include "imgui.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

// ── File-scope particle size for the particle scatter plot ──────────────────
//  Previously there were TWO `static float particle_size = 1.8f;` declarations
//  in two different scopes (one in the render BeginChild, one in the controls
//  BeginChild).  The slider modified the controls-scope one but the render
//  code read the render-scope one — so moving the slider had no effect.
//  Now there's one file-scope static shared by both.
static float g_particle_size = 1.8f;

// ─── Colour palette (mirror of main.cpp's Col namespace) ─────────────────────
namespace VizCol {
    static const ImVec4 BG         = {0.04f, 0.06f, 0.04f, 1.f};
    static const ImVec4 PANEL_DARK = {0.04f, 0.05f, 0.04f, 1.f};
    static const ImVec4 GREEN      = {0.25f, 0.95f, 0.35f, 1.f};
    static const ImVec4 GREEN_DIM  = {0.14f, 0.55f, 0.20f, 1.f};
    static const ImVec4 AMBER      = {0.95f, 0.65f, 0.10f, 1.f};
    static const ImVec4 RED        = {0.95f, 0.18f, 0.18f, 1.f};
    static const ImVec4 CYAN       = {0.20f, 0.85f, 0.85f, 1.f};
    static const ImVec4 WHITE      = {0.88f, 0.92f, 0.88f, 1.f};
    static const ImVec4 GREY       = {0.35f, 0.42f, 0.35f, 1.f};
    static const ImVec4 YELLOW     = {0.85f, 0.85f, 0.20f, 1.f};
    static const ImVec4 BLUE       = {0.30f, 0.55f, 0.95f, 1.f};
    static const ImVec4 ORANGE     = {0.95f, 0.45f, 0.10f, 1.f};
}

static void Hdr(const char* t){
    ImGui::Spacing();ImGui::TextColored(VizCol::GREEN,"[ %s ]",t);
    ImGui::PushStyleColor(ImGuiCol_Separator,VizCol::GREEN_DIM);ImGui::Separator();ImGui::PopStyleColor();
}

static void Row(const char* lbl, float v, const char* fmt, const char* u,
                ImVec4 col=VizCol::GREEN){
    ImGui::TextColored(VizCol::GREY,"%-22s",lbl);ImGui::SameLine();
    char b[32];snprintf(b,32,fmt,v);ImGui::TextColored(col,"%-12s",b);
    ImGui::SameLine();ImGui::TextColored(VizCol::GREY,"%s",u);
}

// ─── Jet colourmap ───────────────────────────────────────────────────────────
//  Maps t∈[0,1] → RGBA.  Standard jet: blue → cyan → green → yellow → red.
static ImU32 jetColor(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    float r, g, b;
    if (t < 0.25f) {
        r = 0.0f;
        g = 4.0f * t;
        b = 1.0f;
    } else if (t < 0.5f) {
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - 4.0f * (t - 0.25f);
    } else if (t < 0.75f) {
        r = 4.0f * (t - 0.5f);
        g = 1.0f;
        b = 0.0f;
    } else {
        r = 1.0f;
        g = 1.0f - 4.0f * (t - 0.75f);
        b = 0.0f;
    }
    return IM_COL32((int)(r * 255), (int)(g * 255), (int)(b * 255), 255);
}

// ─── Inferno colourmap (perceptually uniform, alternative) ───────────────────
static ImU32 infernoColor(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    // Polynomial approximation of the inferno colormap
    float r = std::min(1.0f, 1.5f * t);
    float g = std::max(0.0f, std::min(1.0f, 1.2f * t - 0.2f));
    float b = std::max(0.0f, std::min(1.0f, 1.5f * t * (1.0f - t) * 1.8f));
    return IM_COL32((int)(r * 255), (int)(g * 255), (int)(b * 255), 255);
}

// ═══════════════════════════════════════════════════════════════════════════════
// TAB 8: TEMPERATURE DISTRIBUTION
// ═══════════════════════════════════════════════════════════════════════════════
//
//  Renders the plasma temperature as a 2D heatmap in the R-Z plane
//  (poloidal cross-section through the magnetic axis).  The tokamak's toroidal
//  symmetry means this cross-section is representative of the full 3D plasma.
//
//  Layout:
//    - Large heatmap (left): R on x-axis (R_major - a → R_major + a),
//                            Z on y-axis (-κa → +κa)
//    - Sidebar (right): colour-bar legend, current min/max T, profile plot
//                       T(r) along the midplane (z=0).
//
void TabTemperature(const ReactorState& s, const PlasmaViz& viz)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float ph = ImGui::GetContentRegionAvail().y;

    // ── Left: heatmap canvas ─────────────────────────────────────────────────
    float canvas_w = pw * 0.65f;
    float canvas_h = ph - 30.0f;
    if (ImGui::BeginChild("##temp_canvas", {canvas_w, canvas_h}, true)) {
        ImVec2 p0 = ImGui::GetCursorScreenPos();
        ImVec2 p1 = {p0.x + canvas_w, p0.y + canvas_h};
        ImDrawList* dl = ImGui::GetWindowDrawList();

        // Background
        dl->AddRectFilled(p0, p1, IM_COL32(15, 18, 15, 255));

        // ── Tokamak geometry ──────────────────────────────────────────────
        //  World: x = R ∈ [R_major - a, R_major + a]  → canvas_x
        //         y = Z ∈ [-κ·a, +κ·a]                → canvas_y (flipped)
        const float R_major = viz.R_major();
        const float a       = viz.a_minor();
        const float kappa   = viz.kappa();
        const int   NR = viz.gridNR();
        const int   NZ = viz.gridNZ();

        // Find max temperature for normalization (avoid divide-by-zero)
        float T_max = 1e-3f;
        for (float v : viz.tempGrid()) T_max = std::max(T_max, v);

        // ── Render heatmap as filled rectangles ───────────────────────────
        //  Each grid cell is a rectangle of size (cell_w × cell_h).
        //  Cells outside the plasma boundary (r > a) are skipped.
        float cell_w = canvas_w / NR;
        float cell_h = canvas_h / NZ;

        for (int j = 0; j < NZ; j++) {
            // World Z for this row (top = +κa, bottom = -κa)
            float z_world = (1.0f - 2.0f * (j + 0.5f) / NZ) * kappa * a;
            for (int i = 0; i < NR; i++) {
                // World R for this column (left = R-a, right = R+a)
                float r_world = R_major - a + 2.0f * a * (i + 0.5f) / NR;
                // Minor radius: r_eff = sqrt((R - R_major)² + (z/κ)²)
                float dr = r_world - R_major;
                float r_eff = std::sqrt(dr * dr + (z_world / kappa) * (z_world / kappa));
                if (r_eff > a) continue;  // outside plasma

                int idx = i + NR * j;
                float T = viz.tempGrid()[idx];
                float t_norm = T / T_max;
                ImU32 col = jetColor(t_norm);

                float x0 = p0.x + i * cell_w;
                float y0 = p0.y + j * cell_h;
                float x1 = x0 + cell_w + 1.0f;
                float y1 = y0 + cell_h + 1.0f;
                dl->AddRectFilled({x0, y0}, {x1, y1}, col);
            }
        }

        // ── Overlay: plasma boundary ellipse ──────────────────────────────
        //  Center at (R_major, 0), semi-axes (a, κ·a)
        ImVec2 center = {p0.x + canvas_w * 0.5f, p0.y + canvas_h * 0.5f};
        float ellipse_a = canvas_w * 0.5f;
        float ellipse_b = canvas_h * 0.5f;
        dl->AddCircle(center, std::min(ellipse_a, ellipse_b),
                      IM_COL32(255, 255, 255, 180), 64, 1.5f);
        // Actually use AddEllipse (or simulate with N segments)
        // ImGui's ImDrawList doesn't have AddEllipse directly — use a polygon
        const int N_ELLIPSE = 64;
        for (int k = 0; k < N_ELLIPSE; k++) {
            float a1 = 2.0f * 3.14159265f * k / N_ELLIPSE;
            float a2 = 2.0f * 3.14159265f * (k + 1) / N_ELLIPSE;
            ImVec2 p_a = {center.x + ellipse_a * std::cos(a1),
                          center.y + ellipse_b * std::sin(a1)};
            ImVec2 p_b = {center.x + ellipse_a * std::cos(a2),
                          center.y + ellipse_b * std::sin(a2)};
            dl->AddLine(p_a, p_b, IM_COL32(255, 255, 255, 200), 1.5f);
        }

        // ── Overlay: magnetic axis (R_major, 0) ───────────────────────────
        dl->AddCircleFilled(center, 4.0f, IM_COL32(255, 255, 255, 255));
        dl->AddText({center.x + 8.0f, center.y - 8.0f},
                    IM_COL32(255, 255, 255, 255), "magnetic axis");

        // ── Axis labels ───────────────────────────────────────────────────
        ImGui::SetCursorScreenPos({p0.x + 4, p0.y + canvas_h - 18});
        ImGui::TextColored(VizCol::WHITE, "R [m]: %.1f ← R_major=%.1f → %.1f",
                          R_major - a, R_major, R_major + a);
        ImGui::SetCursorScreenPos({p0.x + 4, p0.y + 4});
        ImGui::TextColored(VizCol::WHITE, "Z [m]: +%.1f (top) to -%.1f (bottom), κ=%.2f",
                          kappa * a, kappa * a, kappa);
    }
    ImGui::EndChild();

    ImGui::SameLine();

    // ── Right: legend + profile ───────────────────────────────────────────────
    if (ImGui::BeginChild("##temp_legend", {0, canvas_h}, true)) {
        Hdr("TEMPERATURE DISTRIBUTION");

        // Find T_max and T_avg
        float T_max = 0.0f, T_sum = 0.0f;
        int N_in = 0;
        for (float v : viz.tempGrid()) {
            if (v > 0) { T_max = std::max(T_max, v); T_sum += v; N_in++; }
        }
        float T_avg = N_in ? T_sum / N_in : 0.0f;

        Row("On-axis T",   s.electron_temp_keV, "%.2f", "keV", VizCol::AMBER);
        Row("Peak T (map)", T_max,              "%.2f", "keV", VizCol::RED);
        Row("Avg T (map)",  T_avg,              "%.2f", "keV");
        Row("Edge T (ped)", T_max * 0.30f,      "%.2f", "keV", VizCol::CYAN);
        Row("Profile shape", 1.3f,              "%.2f", "(α_T)");

        ImGui::Spacing();
        ImGui::TextColored(VizCol::GREY, " Colour scale (keV):");

        // Vertical colour bar
        float bar_w = 18.0f, bar_h = 140.0f;
        ImVec2 bar_p0 = ImGui::GetCursorScreenPos();
        for (int j = 0; j < 64; j++) {
            float t = 1.0f - (float)j / 63.0f;
            ImU32 col = jetColor(t);
            float y0 = bar_p0.y + j * (bar_h / 64.0f);
            float y1 = y0 + bar_h / 64.0f + 1.0f;
            ImGui::GetWindowDrawList()->AddRectFilled(
                {bar_p0.x, y0}, {bar_p0.x + bar_w, y1}, col);
        }
        // Tick labels
        ImGui::SetCursorScreenPos({bar_p0.x + bar_w + 6, bar_p0.y - 4});
        ImGui::TextColored(VizCol::WHITE, "%.2f", T_max);
        ImGui::SetCursorScreenPos({bar_p0.x + bar_w + 6, bar_p0.y + bar_h * 0.5f - 6});
        ImGui::TextColored(VizCol::WHITE, "%.2f", T_max * 0.5f);
        ImGui::SetCursorScreenPos({bar_p0.x + bar_w + 6, bar_p0.y + bar_h - 8});
        ImGui::TextColored(VizCol::WHITE, "0.00");

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();

        // Midplane T(r) profile plot
        Hdr("MIDPLANE PROFILE T(r)");
        ImGui::TextColored(VizCol::GREY, " T along z=0, r∈[0, a]");

        // Build the profile from the grid's midplane row
        const int NR = viz.gridNR();
        const int NZ = viz.gridNZ();
        int j_mid = NZ / 2;
        std::vector<float> profile(NR);
        for (int i = 0; i < NR; i++) {
            profile[i] = viz.tempGrid()[i + NR * j_mid];
        }
        float pmax = *std::max_element(profile.begin(), profile.end()) + 1e-3f;
        ImGui::PushStyleColor(ImGuiCol_PlotLines, VizCol::AMBER);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, VizCol::PANEL_DARK);
        ImGui::PlotLines("##Tprofile", profile.data(), NR, 0, nullptr,
                         0.0f, pmax, {240.0f, 80.0f});
        ImGui::PopStyleColor(2);

        ImGui::TextColored(VizCol::GREY, " ← r=0 (axis)         r=a (edge) →");
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// TAB 9: FUSION REACTIONS OVER TIME
// ═══════════════════════════════════════════════════════════════════════════════
//
//  Multi-panel time-series plot showing:
//    - P_fusion (total), P_alpha (heating), P_neutron (lost to blanket)
//    - Reaction rate (reactions/s)
//    - Cumulative reactions and cumulative energy produced
//    - Current 0D state (T, n, <σv>, Q) for context
//
//  The time-series uses the PlasmaViz history buffers (updated each tick).
//
void TabFusionReactions(const ReactorState& s, const PlasmaViz& viz)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float ph = ImGui::GetContentRegionAvail().y;

    // ── Top: time-series plots (full width) ──────────────────────────────────
    if (ImGui::BeginChild("##fr_plots", {0, ph * 0.6f}, true)) {
        Hdr("FUSION POWER & REACTION RATE — TIME SERIES");

        const auto& hp = viz.histFusionPower();
        const auto& ha = viz.histAlphaPower();
        const auto& hn = viz.histNeutronPower();
        const auto& hr = viz.histReactionRate();

        // Find max for auto-scaling
        float pmax = 1.0f;
        for (float v : hp) pmax = std::max(pmax, v);
        for (float v : ha) pmax = std::max(pmax, v);
        for (float v : hn) pmax = std::max(pmax, v);
        pmax *= 1.1f;

        float rmax = 1.0f;
        for (float v : hr) rmax = std::max(rmax, v);
        rmax *= 1.1f;

        // P_fusion (green, total)
        ImGui::TextColored(VizCol::GREEN, " ● Fusion Power [MW]  (peak=%.1f)",
                          hp.empty() ? 0.0f : *std::max_element(hp.begin(), hp.end()));
        ImGui::PushStyleColor(ImGuiCol_PlotLines, VizCol::GREEN);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, VizCol::PANEL_DARK);
        ImGui::PlotLines("##pfus", hp.data(), (int)hp.size(), 0, nullptr,
                         0.0f, pmax, {pw - 60.f, 70.f});
        ImGui::PopStyleColor(2);

        // P_alpha (amber, plasma self-heating)
        ImGui::TextColored(VizCol::AMBER, " ● Alpha Heating Power [MW]  (3.52/17.59 of P_fus)");
        ImGui::PushStyleColor(ImGuiCol_PlotLines, VizCol::AMBER);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, VizCol::PANEL_DARK);
        ImGui::PlotLines("##palpha", ha.data(), (int)ha.size(), 0, nullptr,
                         0.0f, pmax, {pw - 60.f, 70.f});
        ImGui::PopStyleColor(2);

        // P_neutron (cyan, deposited in blanket)
        ImGui::TextColored(VizCol::CYAN, " ● Neutron Power [MW]  (14.07/17.59 → blanket)");
        ImGui::PushStyleColor(ImGuiCol_PlotLines, VizCol::CYAN);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, VizCol::PANEL_DARK);
        ImGui::PlotLines("##pneut", hn.data(), (int)hn.size(), 0, nullptr,
                         0.0f, pmax, {pw - 60.f, 70.f});
        ImGui::PopStyleColor(2);

        // Reaction rate (yellow)
        ImGui::TextColored(VizCol::YELLOW, " ● Reaction Rate [reactions/s]");
        ImGui::PushStyleColor(ImGuiCol_PlotLines, VizCol::YELLOW);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, VizCol::PANEL_DARK);
        ImGui::PlotLines("##rate", hr.data(), (int)hr.size(), 0, nullptr,
                         0.0f, rmax, {pw - 60.f, 70.f});
        ImGui::PopStyleColor(2);
    }
    ImGui::EndChild();

    // ── Bottom: 0D diagnostics + cumulative stats ────────────────────────────
    if (ImGui::BeginChild("##fr_stats", {0, 0}, true)) {
        Hdr("INSTANTANEOUS FUSION DIAGNOSTICS");

        float pw2 = ImGui::GetContentRegionAvail().x;
        if (ImGui::BeginChild("##fr_left", {pw2 * 0.5f, 0.f}, false)) {
            Row("Fusion Power",    s.fusion_power_MW,   "%.2f", "MW",   VizCol::GREEN);
            Row("Alpha Heating",   s.alpha_power_MW,    "%.2f", "MW",   VizCol::AMBER);
            Row("Neutron Power",   s.fusion_power_MW * (14.070f/17.589f), "%.2f", "MW", VizCol::CYAN);
            Row("Radiated Power",  s.radiated_power_MW, "%.2f", "MW",   VizCol::GREY);

            // Reaction rate from P_fusion / E_fusion_per_reaction
            constexpr float E_fus_J = 17.59e6f * 1.602176634e-19f;  // 17.59 MeV in J
            float rate_per_s = s.fusion_power_MW * 1e6f / E_fus_J;
            Row("Reaction Rate",   rate_per_s,           "%.3e", "1/s", VizCol::YELLOW);

            Row("Q scientific",    s.Q_scientific,       "%.3f", "",    VizCol::GREEN);
            Row("Neutron Flux",    s.neutron_flux_m2s,   "%.3e", "n/m²s", VizCol::CYAN);

            ImGui::Spacing();
            // <σv> from the verified Bosch-Hale fit
            float sv = ConfinementPhysics::boschHaleSigmaV_DT(s.electron_temp_keV);
            Row("<σv>(T) Bosch-Hale", sv, "%.3e", "m³/s", VizCol::AMBER);
        }
        ImGui::EndChild();
        ImGui::SameLine();

        if (ImGui::BeginChild("##fr_right", {0, 0}, false)) {
            Hdr("CUMULATIVE TOTALS");

            // Approximate cumulative reactions from the time-series (trapezoidal)
            const auto& hp = viz.histFusionPower();
            const auto& ht = viz.histTime();
            double cumul_energy_MJ = 0.0;
            double cumul_reactions = 0.0;
            constexpr float E_fus_J = 17.59e6f * 1.602176634e-19f;
            if (hp.size() > 1) {
                // Assume each sample is 1 ms (SIM_DT_S = 1e-3)
                const float dt = 1e-3f;
                for (size_t i = 0; i < hp.size(); i++) {
                    cumul_energy_MJ  += (double)hp[i] * dt;  // MW × s = MJ
                    cumul_reactions  += (double)hp[i] * 1e6f * dt / E_fus_J;
                }
            }

            Row("Total reactions", (float)cumul_reactions, "%.3e", "count", VizCol::YELLOW);
            Row("Energy produced", (float)cumul_energy_MJ, "%.2f",  "MJ",    VizCol::GREEN);
            Row("Alpha energy",    (float)(cumul_energy_MJ * 3.521f / 17.589f), "%.2f", "MJ", VizCol::AMBER);
            Row("Neutron energy",  (float)(cumul_energy_MJ * 14.070f / 17.589f), "%.2f", "MJ", VizCol::CYAN);

            ImGui::Spacing();
            // Power-balance check.
            //  Use the actual aux heating from the H&CD systems rather than
            //  the old hardcoded 50 MW (which was left over from the original
            //  pre-H&CD power balance model).
            Hdr("POWER BALANCE");
            float P_aux = s.hcd_nbi_actual_MW + s.hcd_icrh_actual_MW
                        + s.hcd_ecrh_actual_MW + s.hcd_lhcd_actual_MW;
            float P_in  = s.alpha_power_MW + P_aux;
            float P_out = s.radiated_power_MW + s.fusion_power_MW * (14.070f/17.589f);
            Row("P_in  (α + aux)",  P_in,  "%.2f", "MW", VizCol::GREEN);
            Row("P_out (rad + n)",  P_out, "%.2f", "MW", VizCol::RED);
            Row("Net heating",      P_in - P_out, "%.2f", "MW",
                (P_in - P_out) > 0 ? VizCol::GREEN : VizCol::AMBER);
        }
        ImGui::EndChild();
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// TAB 10: PARTICLE DISTRIBUTION
// ═══════════════════════════════════════════════════════════════════════════════
//
//  Scatter plot of macro-particles in the plasma.  The user can:
//    - Toggle each species (e, D, T, α, p, ³He) on/off via checkboxes
//    - Switch projection: R-Z (poloidal cross-section) or R-φ (top-down)
//    - Adjust particle point size
//    - View per-species statistics (count, mean energy, abundance)
//
//  The particles come from PlasmaViz, which samples N_VIZ_PARTICLES per
//  species from the local Maxwellian at each radial position.
//
void TabParticleDistribution(const ReactorState& s, const PlasmaViz& viz,
                              uint32_t& visible_species_mask,
                              int& projection_mode)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float ph = ImGui::GetContentRegionAvail().y;

    // ── Left: scatter canvas ──────────────────────────────────────────────────
    float canvas_w = pw * 0.7f;
    float canvas_h = ph - 30.0f;
    if (ImGui::BeginChild("##pd_canvas", {canvas_w, canvas_h}, true)) {
        ImVec2 p0 = ImGui::GetCursorScreenPos();
        ImVec2 p1 = {p0.x + canvas_w, p0.y + canvas_h};
        ImDrawList* dl = ImGui::GetWindowDrawList();

        // Background
        dl->AddRectFilled(p0, p1, IM_COL32(8, 10, 14, 255));

        const float R_major = viz.R_major();
        const float a       = viz.a_minor();
        const float kappa   = viz.kappa();

        // ── World → canvas transform ──────────────────────────────────────
        //  R-Z mode:    x = R, y = Z       → x∈[R-a, R+a], y∈[-κa, κa]
        //  R-φ mode:    x = R, y = Y (=cartesian y of particle) → top-down
        float world_x_min, world_x_max, world_y_min, world_y_max;
        if (projection_mode == 0) {
            // R-Z poloidal cross-section
            world_x_min = R_major - a;
            world_x_max = R_major + a;
            world_y_min = -kappa * a;
            world_y_max = +kappa * a;
        } else {
            // Top-down (R-φ): particles' x-y in vessel frame
            world_x_min = R_major - a;
            world_x_max = R_major + a;
            world_y_min = -a;
            world_y_max = +a;
        }
        float wx_range = world_x_max - world_x_min;
        float wy_range = world_y_max - world_y_min;

        // Compute scale to fit (preserve aspect ratio)
        float scale_x = canvas_w / wx_range;
        float scale_y = canvas_h / wy_range;
        float scale = std::min(scale_x, scale_y) * 0.92f;
        float ox = p0.x + canvas_w * 0.5f;
        float oy = p0.y + canvas_h * 0.5f;

        auto world_to_canvas = [&](float wx, float wy) -> ImVec2 {
            return {ox + (wx - 0.5f * (world_x_min + world_x_max)) * scale,
                    oy - (wy - 0.5f * (world_y_min + world_y_max)) * scale};
        };

        // ── Plasma boundary ───────────────────────────────────────────────
        //  Draw the boundary as an ellipse for R-Z (semi-axes a, κ·a) or a
        //  circle for top-down (radius a, centred at R=R_major).
        const int N_ELLIPSE = 64;
        if (projection_mode == 0) {
            // Poloidal: ellipse centred at (R_major, 0)
            for (int k = 0; k < N_ELLIPSE; k++) {
                float a1 = 2.0f * 3.14159265f * k / N_ELLIPSE;
                float a2 = 2.0f * 3.14159265f * (k + 1) / N_ELLIPSE;
                ImVec2 p_a = world_to_canvas(R_major + a * std::cos(a1),
                                              kappa * a * std::sin(a1));
                ImVec2 p_b = world_to_canvas(R_major + a * std::cos(a2),
                                              kappa * a * std::sin(a2));
                dl->AddLine(p_a, p_b, IM_COL32(180, 220, 180, 200), 1.5f);
            }
        } else {
            // Top-down: circle at (R_major, 0)
            for (int k = 0; k < N_ELLIPSE; k++) {
                float a1 = 2.0f * 3.14159265f * k / N_ELLIPSE;
                float a2 = 2.0f * 3.14159265f * (k + 1) / N_ELLIPSE;
                ImVec2 p_a = world_to_canvas(R_major + a * std::cos(a1),
                                              a * std::sin(a1));
                ImVec2 p_b = world_to_canvas(R_major + a * std::cos(a2),
                                              a * std::sin(a2));
                dl->AddLine(p_a, p_b, IM_COL32(180, 220, 180, 200), 1.5f);
            }
        }

        // ── Magnetic axis ─────────────────────────────────────────────────
        ImVec2 axis_pt = world_to_canvas(R_major, 0.0f);
        dl->AddCircleFilled(axis_pt, 3.0f, IM_COL32(255, 255, 255, 200));
        dl->AddText({axis_pt.x + 6, axis_pt.y - 6},
                    IM_COL32(255, 255, 255, 200), "axis");

        // ── Particle scatter ──────────────────────────────────────────────
        //  Each species gets a colour; only render species whose bit is set
        //  in visible_species_mask.
        //
        //  We render species in order of abundance (electrons first, then
        //  D, T, α) so the rarer species appear on top.
        //  (particle_size is the file-scope g_particle_size, set by the
        //  "Point size" slider in the controls panel below.)
        int total_rendered = 0;

        for (int sp = 0; sp < PlasmaViz::N_SPECIES; sp++) {
            if (!((visible_species_mask >> sp) & 1)) continue;

            float col[4];
            PlasmaViz::speciesColor((PlasmaViz::Species)sp, col);
            ImU32 icol = IM_COL32((int)(col[0] * 255), (int)(col[1] * 255),
                                   (int)(col[2] * 255), 200);

            const auto& parts = viz.particles((PlasmaViz::Species)sp);
            for (const auto& P : parts) {
                if (P.weight < 1e-6f) continue;  // skip particles of absent species

                ImVec2 pt;
                if (projection_mode == 0) {
                    // R-Z: x = R = particle.x, y = Z = particle.z
                    pt = world_to_canvas(P.x, P.z);
                } else {
                    // Top-down: x = R = particle.x, y = particle.y
                    pt = world_to_canvas(P.x, P.y);
                }

                if (pt.x < p0.x || pt.x > p1.x || pt.y < p0.y || pt.y > p1.y)
                    continue;

                // Vary point size by weight (heavier = bigger)
                float size = g_particle_size * (0.5f + 1.5f * std::min(P.weight * 1e-19f, 1.0f));
                dl->AddCircleFilled(pt, size, icol);
                total_rendered++;
            }
        }

        // ── Axis labels & overlay ─────────────────────────────────────────
        ImGui::SetCursorScreenPos({p0.x + 4, p0.y + 4});
        const char* proj_name = (projection_mode == 0) ? "R-Z (poloidal)" : "R-φ (top-down)";
        ImGui::TextColored(VizCol::WHITE, "Projection: %s   |   Particles rendered: %d",
                          proj_name, total_rendered);

        ImGui::SetCursorScreenPos({p0.x + 4, p0.y + canvas_h - 18});
        ImGui::TextColored(VizCol::WHITE, "x: R∈[%.1f, %.1f] m   y: %s∈[%.1f, %.1f] m",
                          world_x_min, world_x_max,
                          projection_mode == 0 ? "Z" : "Y",
                          world_y_min, world_y_max);
    }
    ImGui::EndChild();

    ImGui::SameLine();

    // ── Right: species toggles + stats ────────────────────────────────────────
    if (ImGui::BeginChild("##pd_controls", {0, canvas_h}, true)) {
        Hdr("SPECIES VISIBILITY");

        const char* sp_names[] = {
            "Electrons (e-)",
            "Deuterium (D+)",
            "Tritium (T+)",
            "Alpha (α = He²+)",
            "Protons (p+)",
            "Helium-3 (³He²+)"
        };
        for (int sp = 0; sp < PlasmaViz::N_SPECIES; sp++) {
            bool vis = (visible_species_mask >> sp) & 1;
            float col[4];
            PlasmaViz::speciesColor((PlasmaViz::Species)sp, col);

            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(col[0], col[1], col[2], col[3]));
            if (ImGui::Checkbox(sp_names[sp], &vis)) {
                if (vis) visible_species_mask |= (1u << sp);
                else     visible_species_mask &= ~(1u << sp);
            }
            ImGui::PopStyleColor();
        }

        // Quick-select buttons
        ImGui::Spacing();
        if (ImGui::Button("All"))  visible_species_mask = 0x3F;
        ImGui::SameLine();
        if (ImGui::Button("None")) visible_species_mask = 0x00;
        ImGui::SameLine();
        if (ImGui::Button("D-T only")) visible_species_mask = (1u << 1) | (1u << 2);

        // Projection toggle
        ImGui::Spacing();
        Hdr("PROJECTION");
        ImGui::RadioButton("R-Z (poloidal)", &projection_mode, 0);
        ImGui::RadioButton("R-φ (top-down)",  &projection_mode, 1);

        // Particle size slider — controls the file-scope g_particle_size
        // which the render code above reads.  Previously this was a separate
        // local static, so the slider didn't affect rendering.
        ImGui::Spacing();
        Hdr("RENDERING");
        ImGui::SliderFloat("Point size", &g_particle_size, 0.5f, 5.0f);

        // Per-species statistics
        ImGui::Spacing();
        Hdr("SPECIES STATS");

        // Mass table (kg)
        static const float m_kg[PlasmaViz::N_SPECIES] = {
            9.1093837015e-31f, 3.3435837724e-27f, 5.0073558862e-27f,
            6.6446573357e-27f, 1.67262192369e-27f, 5.0082343773e-27f
        };

        for (int sp = 0; sp < PlasmaViz::N_SPECIES; sp++) {
            float col[4];
            PlasmaViz::speciesColor((PlasmaViz::Species)sp, col);
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(col[0], col[1], col[2], col[3]));

            const auto& parts = viz.particles((PlasmaViz::Species)sp);
            // Mean KE = (1/2) m <v²>; for Maxwellian <v²> = 3 kT/m so KE = (3/2) kT
            // but here we use the actual sampled velocities.
            float sum_KE = 0.0f;
            int N_active = 0;
            for (const auto& P : parts) {
                if (P.weight < 1e-6f) continue;
                float v2 = P.vx*P.vx + P.vy*P.vy + P.vz*P.vz;
                sum_KE += 0.5f * m_kg[sp] * v2;
                N_active++;
            }
            float mean_KE_J  = N_active ? sum_KE / N_active : 0.0f;
            float mean_KE_keV = mean_KE_J / 1.602176634e-16f;

            ImGui::TextColored(ImVec4(col[0], col[1], col[2], col[3]),
                              " %-14s N=%-5d  <KE>=%.2f keV",
                              sp_names[sp], N_active, mean_KE_keV);
            ImGui::PopStyleColor();
        }

        // Tokamak geometry reminder
        ImGui::Spacing();
        Hdr("PLASMA GEOMETRY");
        Row("Major radius R", viz.R_major(), "%.2f", "m");
        Row("Minor radius a", viz.a_minor(), "%.2f", "m");
        Row("Elongation κ",  viz.kappa(),    "%.2f", "");
        Row("Aspect ratio",  viz.R_major() / viz.a_minor(), "%.2f", "");
        Row("Plasma volume", 2.0f * 3.14159265f * viz.R_major() *
                             viz.a_minor() * viz.a_minor() * viz.kappa(),
                             "%.1f", "m³", VizCol::CYAN);
    }
    ImGui::EndChild();
}
