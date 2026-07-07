//
// src/main.cpp  —  FusionSim expanded UI
// Master branch ImGui (no docking). Tab-bar navigation over 8 views.
// Keyboard: SPACE=pause, F1=SCRAM, ESC=quit
//

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <deque>
#include <string>
#include <vector>
#include <chrono>

#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include <SDL.h>
#include <SDL_opengl.h>

// ── FusionSim ─────────────────────────────────────────────────────────────────
#include "ReactorState.h"
#include "confinement_physics.h"   // POPCON operating-space tab reuses the 0D physics
#include "SimTime.h"
#include "Control/Control.h"
#include "Magnets/Magnets.h"
#include "Fuel/Fuel.h"
#include "PlasmaCore/PlasmaCoreBridge.h"
#include "Helium/Helium.h"
#include "ThermalHydraulics/ThermalHydraulics.h"
#include "SteamPower/SteamPower.h"
// ── New modules ───────────────────────────────────────────────────────────────
#include "TurbineSystem/TurbineSystem.h"
#include "ElectricalGrid/ElectricalGrid.h"
#include "MoltenSalt/MoltenSaltSystem.h"
#include "HeliumSystem/HeliumCoolingSystem.h"
// ── Round 4 new modules ───────────────────────────────────────────────────────
#include "HCD/HCD.h"
#include "VacuumVessel/VacuumVessel.h"
#include "DisruptionMitigation/DisruptionMitigation.h"
#include "TritiumPlant/TritiumPlant.h"

// ── Plasma visualization (3 new tabs) ─────────────────────────────────────────
#include "PlasmaVisualization.h"
#include "PlasmaVizTabs.h"

// ─── Layout constants ─────────────────────────────────────────────────────────
static constexpr float SIM_DT_S        = 1e-3f;
static constexpr int   HISTORY_SAMPLES = 512;
static constexpr int   WINDOW_W        = 1680;
static constexpr int   WINDOW_H        = 980;
static constexpr float STATUS_H        = 34.f;
static constexpr float LEFT_W_FRAC     = 0.19f; // control panel width

// ═══════════════════════════════════════════════════════════════════════════════
// COLOUR PALETTE
// ═══════════════════════════════════════════════════════════════════════════════
namespace Col {
    static const ImVec4 BG         = {0.04f, 0.06f, 0.04f, 1.f};
    static const ImVec4 PANEL      = {0.07f, 0.09f, 0.07f, 1.f};
    static const ImVec4 PANEL_DARK = {0.04f, 0.05f, 0.04f, 1.f};
    static const ImVec4 BORDER     = {0.18f, 0.26f, 0.18f, 1.f};
    static const ImVec4 GREEN      = {0.25f, 0.95f, 0.35f, 1.f};
    static const ImVec4 GREEN_DIM  = {0.14f, 0.55f, 0.20f, 1.f};
    static const ImVec4 GREEN_DARK = {0.06f, 0.22f, 0.08f, 1.f};
    static const ImVec4 AMBER      = {0.95f, 0.65f, 0.10f, 1.f};
    static const ImVec4 RED        = {0.95f, 0.18f, 0.18f, 1.f};
    static const ImVec4 CYAN       = {0.20f, 0.85f, 0.85f, 1.f};
    static const ImVec4 WHITE      = {0.88f, 0.92f, 0.88f, 1.f};
    static const ImVec4 GREY       = {0.35f, 0.42f, 0.35f, 1.f};
    static const ImVec4 YELLOW     = {0.85f, 0.85f, 0.20f, 1.f};
    static const ImVec4 BLUE       = {0.30f, 0.55f, 0.95f, 1.f};
    static const ImVec4 ORANGE     = {0.95f, 0.45f, 0.10f, 1.f};
}

// ═══════════════════════════════════════════════════════════════════════════════
// UTILITY TYPES
// ═══════════════════════════════════════════════════════════════════════════════
struct ScrollBuf {
    std::deque<float> data; int cap;
    float ymin=0.f, ymax=1.f;
    explicit ScrollBuf(int c=HISTORY_SAMPLES):cap(c){}
    void push(float v){
        data.push_back(v);
        if((int)data.size()>cap)data.pop_front();
        if(!data.empty()){ymin=*std::min_element(data.begin(),data.end());
                          ymax=*std::max_element(data.begin(),data.end());}
    }
    std::vector<float> flat()const{return{data.begin(),data.end()};}
};

// ═══════════════════════════════════════════════════════════════════════════════
// ALARM SYSTEM (Round 4 — descriptive alarms with cause + remediation)
// ═══════════════════════════════════════════════════════════════════════════════
//  Each alarm now carries:
//    - msg          : short label (e.g. "Disruption risk")
//    - cause        : what triggered it (e.g. "q_95 = 1.82, below the 2.0 kink-mode limit")
//    - remediation  : what the operator should do (e.g. "Reduce plasma current to raise q_95")
//    - severity     : info / warning / critical (drives colour + ack behaviour)
//
//  The OLD AlarmSystem just had a short msg string — the operator could see
//  that "Disruption risk" was active but had no way to know *what* triggered
//  it or *what to do about it*.  The new struct gives the operator the same
//  context a real control-room alarm annunciator panel provides.
enum class AlarmSeverity : uint8_t { Info = 0, Warning = 1, Critical = 2 };

struct AlarmEntry{
    std::string msg;
    std::string cause;
    std::string remediation;
    AlarmSeverity severity = AlarmSeverity::Warning;
    double t = 0.0;
    bool active = false;
    bool acked  = false;
};

struct AlarmSystem{
    std::vector<AlarmEntry> log; bool any_unacked=false;

    // Trip an alarm with a short msg + detailed cause + remediation text.
    // If `active` is true and the alarm is new (or was previously cleared),
    // it gets added to the log with acked=false so the operator sees the
    // blinking "!! ALARM !!" indicator.
    void trip(const std::string& m, bool active, double t,
              const std::string& cause = "",
              const std::string& remediation = "",
              AlarmSeverity sev = AlarmSeverity::Warning){
        for(auto& a:log){
            if(a.msg==m){
                if(active && !a.active){
                    a.acked=false;
                    any_unacked=true;
                    // Update cause/remediation in case the trigger changed
                    if(!cause.empty()) a.cause = cause;
                    if(!remediation.empty()) a.remediation = remediation;
                    a.severity = sev;
                }
                a.active=active;
                a.t = active ? a.t : t;  // keep the original trip time
                return;
            }
        }
        if(active){
            log.push_back({m, cause, remediation, sev, t, true, false});
            any_unacked=true;
        }
    }

    void ackAll(){for(auto& a:log)a.acked=true;any_unacked=false;}
    void clearOld(){log.erase(std::remove_if(log.begin(),log.end(),
        [](const AlarmEntry&a){return!a.active&&a.acked;}),log.end());}

    int countActive() const {
        int n=0; for(const auto& a:log) if(a.active) n++; return n;
    }
    int countCritical() const {
        int n=0; for(const auto& a:log) if(a.active && a.severity==AlarmSeverity::Critical) n++;
        return n;
    }
};

// Helper: map a DisruptionCause enum to a human-readable cause string
static const char* disruptionCauseStr(DisruptionCause c){
    switch(c){
        case DisruptionCause::LowQ:              return "Safety factor q_95 below 2.0 — kink mode risk";
        case DisruptionCause::Greenwald:         return "Plasma density exceeds Greenwald limit — radiative collapse risk";
        case DisruptionCause::Troyon:            return "Normalized beta (β_N) above Troyon limit (2.5) — ideal MHD kink";
        case DisruptionCause::TearingMode:       return "Tearing mode (m/n = 2/1 or 3/2) island grew past critical width";
        case DisruptionCause::VDE:               return "Vertical displacement event — plasma lost radial position control";
        case DisruptionCause::HaloCurrent:       return "Post-disruption halo current exceeds structural design";
        case DisruptionCause::RunawayElectrons:  return "Loop electric field above Dreicer threshold — runaway electron seed";
        case DisruptionCause::DensityLimit:      return "Murphy density limit — radiative collapse from impurity radiation";
        case DisruptionCause::DisruptionImminent:return "MHD module flagged disruption_ongoing (current quench in progress)";
        case DisruptionCause::None: default:     return "None";
    }
}

// Helper: map a DisruptionCause to a recommended operator action
static const char* disruptionCauseRemediation(DisruptionCause c){
    switch(c){
        case DisruptionCause::LowQ:
            return "Reduce plasma current I_p (raises q_95 ∝ B·a²/I_p) or "
                   "increase toroidal field B_T; check for accidental I_p ramp.";
        case DisruptionCause::Greenwald:
            return "Reduce gas puffing and pellet frequency to lower n_e; "
                   "increase plasma current (raises Greenwald limit ∝ I_p/a²).";
        case DisruptionCause::Troyon:
            return "Reduce auxiliary heating (NBI/ICRH/ECRH/LHCD) to lower β; "
                   "or reduce plasma density (β ∝ n·T/B²).";
        case DisruptionCause::TearingMode:
            return "Apply localized ECRH at the q=2 or q=3/2 resonant surface "
                   "to stabilize the island; reduce neutral beam power.";
        case DisruptionCause::VDE:
            return "Vertical position controller has lost control — immediately "
                   "reduce plasma current and prepare for halo current forces. "
                   "Check PF coil currents and Z-position controller.";
        case DisruptionCause::HaloCurrent:
            return "Disruption is in progress — verify MGI/SPI fired; check "
                   "vacuum vessel structural monitoring for overstress.";
        case DisruptionCause::RunawayElectrons:
            return "Increase plasma density (collisional damping suppresses "
                   "runaways) or apply MGI/SPI for collisional suppression; "
                   "avoid high loop voltage during current quench.";
        case DisruptionCause::DensityLimit:
            return "Reduce impurity influx (check divertor, first wall); "
                   "increase pumping; consider wall boronization to getter O.";
        case DisruptionCause::DisruptionImminent:
            return "Disruption is in progress — arm and fire MGI/SPI if not "
                   "already; verify SCRAM interlock engages; prepare for "
                   "halo current forces and thermal quench.";
        case DisruptionCause::None: default:
            return "No action needed.";
    }
}

// ─── ImGui theme ─────────────────────────────────────────────────────────────
static void ApplyTheme(){
    ImGuiStyle& s=ImGui::GetStyle();
    s.WindowRounding=2.f;s.FrameRounding=2.f;s.GrabRounding=2.f;
    s.TabRounding=2.f;s.WindowBorderSize=1.f;s.FrameBorderSize=1.f;
    s.ItemSpacing={6.f,4.f};s.FramePadding={6.f,3.f};s.WindowPadding={8.f,8.f};
    auto*c=s.Colors;
    c[ImGuiCol_WindowBg]            =Col::PANEL;
    c[ImGuiCol_ChildBg]             =Col::PANEL_DARK;
    c[ImGuiCol_PopupBg]             =Col::PANEL;
    c[ImGuiCol_Border]              =Col::BORDER;
    c[ImGuiCol_FrameBg]             =Col::PANEL_DARK;
    c[ImGuiCol_FrameBgHovered]      =Col::GREEN_DARK;
    c[ImGuiCol_FrameBgActive]       ={0.10f,0.30f,0.12f,1.f};
    c[ImGuiCol_TitleBg]             =Col::PANEL_DARK;
    c[ImGuiCol_TitleBgActive]       ={0.06f,0.14f,0.07f,1.f};
    c[ImGuiCol_TitleBgCollapsed]    =Col::PANEL_DARK;
    c[ImGuiCol_ScrollbarBg]         =Col::PANEL_DARK;
    c[ImGuiCol_ScrollbarGrab]       =Col::GREEN_DIM;
    c[ImGuiCol_ScrollbarGrabHovered]=Col::GREEN;
    c[ImGuiCol_ScrollbarGrabActive] =Col::GREEN;
    c[ImGuiCol_CheckMark]           =Col::GREEN;
    c[ImGuiCol_SliderGrab]          =Col::GREEN_DIM;
    c[ImGuiCol_SliderGrabActive]    =Col::GREEN;
    c[ImGuiCol_Button]              ={0.10f,0.22f,0.11f,1.f};
    c[ImGuiCol_ButtonHovered]       ={0.15f,0.38f,0.17f,1.f};
    c[ImGuiCol_ButtonActive]        ={0.20f,0.55f,0.23f,1.f};
    c[ImGuiCol_Header]              =Col::GREEN_DARK;
    c[ImGuiCol_HeaderHovered]       ={0.12f,0.30f,0.14f,1.f};
    c[ImGuiCol_HeaderActive]        ={0.18f,0.45f,0.20f,1.f};
    c[ImGuiCol_Separator]           =Col::BORDER;
    c[ImGuiCol_Tab]                 =Col::PANEL_DARK;
    c[ImGuiCol_TabHovered]          ={0.12f,0.30f,0.14f,1.f};
    c[ImGuiCol_TabActive]           ={0.10f,0.22f,0.11f,1.f};
    c[ImGuiCol_TabUnfocused]        =Col::PANEL_DARK;
    c[ImGuiCol_TabUnfocusedActive]  =Col::PANEL_DARK;
    c[ImGuiCol_PlotLines]           =Col::GREEN;
    c[ImGuiCol_PlotLinesHovered]    =Col::AMBER;
    c[ImGuiCol_PlotHistogram]       =Col::GREEN_DIM;
    c[ImGuiCol_Text]                =Col::WHITE;
    c[ImGuiCol_TextDisabled]        =Col::GREY;
    c[ImGuiCol_NavHighlight]        =Col::GREEN;
}

// ─── UI helpers ───────────────────────────────────────────────────────────────
static void Row(const char* lbl, float v, const char* fmt, const char* u,
                ImVec4 col=Col::GREEN){
    ImGui::TextColored(Col::GREY,"%-22s",lbl);ImGui::SameLine();
    char b[32];snprintf(b,32,fmt,v);ImGui::TextColored(col,"%-12s",b);
    ImGui::SameLine();ImGui::TextColored(Col::GREY,"%s",u);
}
static void Hdr(const char* t){
    ImGui::Spacing();ImGui::TextColored(Col::GREEN,"[ %s ]",t);
    ImGui::PushStyleColor(ImGuiCol_Separator,Col::GREEN_DIM);ImGui::Separator();ImGui::PopStyleColor();
}
static void Bar(float f,float w,ImVec4 lo,ImVec4 hi,float thr=0.8f){
    ImGui::PushStyleColor(ImGuiCol_PlotHistogram,f>thr?hi:lo);
    ImGui::ProgressBar(std::clamp(f,0.f,1.f),{w,10.f},"");ImGui::PopStyleColor();
}
static void Plot(const char* id,const ScrollBuf& b,float w,float h,ImVec4 col){
    auto f=b.flat();if(f.empty())return;
    float lo=b.ymin*.9f,hi=b.ymax*1.1f+1e-4f;
    ImGui::PushStyleColor(ImGuiCol_PlotLines,col);
    ImGui::PushStyleColor(ImGuiCol_FrameBg,Col::PANEL_DARK);
    ImGui::PlotLines(id,f.data(),(int)f.size(),0,nullptr,lo,hi,{w,h});
    ImGui::PopStyleColor(2);
}
static void Lamp(const char* lbl,bool active){
    ImVec4 c=active?Col::RED:Col::GREEN_DARK;
    ImGui::PushStyleColor(ImGuiCol_Button,c);ImGui::PushStyleColor(ImGuiCol_ButtonHovered,c);ImGui::PushStyleColor(ImGuiCol_ButtonActive,c);
    ImGui::SmallButton(lbl);ImGui::PopStyleColor(3);
}
static bool GreenBtn(const char* lbl,ImVec2 sz={0.f,0.f}){
    ImGui::PushStyleColor(ImGuiCol_Button,{0.06f,0.28f,0.10f,1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered,{0.10f,0.45f,0.16f,1.f});
    bool r=ImGui::Button(lbl,sz);ImGui::PopStyleColor(2);return r;
}
static bool RedBtn(const char* lbl,ImVec2 sz={0.f,0.f}){
    ImGui::PushStyleColor(ImGuiCol_Button,{0.40f,0.05f,0.05f,1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered,{0.65f,0.08f,0.08f,1.f});
    bool r=ImGui::Button(lbl,sz);ImGui::PopStyleColor(2);return r;
}
static bool BeginTiled(const char* n,float x,float y,float w,float h,ImGuiWindowFlags ex=0){
    ImGui::SetNextWindowPos({x,y},ImGuiCond_Always);
    ImGui::SetNextWindowSize({w,h},ImGuiCond_Always);
    return ImGui::Begin(n,nullptr,
        ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoResize|
        ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoBringToFrontOnFocus|ex);
}

// ─── Helper: turbine state string + colour ────────────────────────────────────
static const char* turbStateStr(TurbineState st){
    switch(st){
        case TurbineState::Offline:      return "OFFLINE";
        case TurbineState::RollingUp:    return "ROLLING UP";
        case TurbineState::Synchronizing:return "SYNCHRONIZING";
        case TurbineState::Online:       return "ONLINE";
        case TurbineState::Runback:      return "RUNBACK";
        case TurbineState::Tripping:     return "TRIPPING";
        case TurbineState::Tripped:      return "TRIPPED";
    } return "?";
}
static ImVec4 turbStateCol(TurbineState st){
    switch(st){
        case TurbineState::Online:       return Col::GREEN;
        case TurbineState::RollingUp:
        case TurbineState::Synchronizing:return Col::AMBER;
        case TurbineState::Tripping:
        case TurbineState::Tripped:      return Col::RED;
        default:                         return Col::GREY;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// STATUS BAR
// ═══════════════════════════════════════════════════════════════════════════════
static void RenderStatusBar(const ReactorState& s,const SimTime& t,
                             bool paused,const AlarmSystem& alm,float speed){
    ImGuiIO& io=ImGui::GetIO();
    ImGui::SetNextWindowPos({0,0},ImGuiCond_Always);
    ImGui::SetNextWindowSize({io.DisplaySize.x,STATUS_H},ImGuiCond_Always);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding,0.f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize,0.f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg,Col::PANEL_DARK);
    ImGui::Begin("##sb",nullptr,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|
        ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoBringToFrontOnFocus);

    const char* ms="STARTUP"; ImVec4 mc=Col::AMBER;
    switch(s.mode){case ReactorMode::SteadyState:ms="STEADY STATE";mc=Col::GREEN;break;
                   case ReactorMode::Emergency:ms="EMERGENCY";mc=Col::RED;break;
                   case ReactorMode::Rampdown:ms="RAMP DOWN";break;default:break;}
    ImGui::TextColored(mc,"  ■ %s",ms); ImGui::SameLine(190.f);

    const char* ps="COLD"; ImVec4 pc=Col::GREY;
    switch(s.plasma_status){
        case PlasmaStatus::Initiating:ps="INITIATING";pc=Col::AMBER;break;
        case PlasmaStatus::Burning:   ps="BURNING ▲"; pc=Col::GREEN;break;
        case PlasmaStatus::Disrupting:ps="DISRUPTION";pc=Col::RED;  break;
        case PlasmaStatus::Quenched:  ps="QUENCHED";  pc=Col::RED;  break;default:break;}
    ImGui::TextColored(pc,"%-11s",ps); ImGui::SameLine(350.f);
    // Confinement regime chip — H-mode is the money regime, L-mode is amber
    // as a nudge that the operator is paying double on confinement.
    if(s.plasma_status==PlasmaStatus::Burning||s.plasma_status==PlasmaStatus::Initiating){
        ImGui::TextColored(s.h_mode?Col::GREEN:Col::AMBER, s.h_mode?"[H]":"[L]");
    } else {
        ImGui::TextColored(Col::GREY,"[-]");
    }
    ImGui::SameLine(395.f);
    ImGui::TextColored(Col::GREEN,  "P_fus %.0f MW",  s.fusion_power_MW);   ImGui::SameLine(520.f);
    ImGui::TextColored(Col::YELLOW, "P_net %.0f MW",  s.net_electric_MW);   ImGui::SameLine(650.f);
    ImGui::TextColored(Col::CYAN,   "Q=%.2f",         s.Q_scientific);      ImGui::SameLine(745.f);
    ImGui::TextColored(s.grid_frequency_Hz<49.f?Col::AMBER:Col::GREEN,
                       "f=%.2f Hz", s.grid_frequency_Hz);                    ImGui::SameLine(840.f);
    ImGui::TextColored(paused?Col::AMBER:Col::GREEN,paused?"■■ PAUSE":"▶ %.0fx",speed); ImGui::SameLine(960.f);
    ImGui::TextColored(Col::GREY,"t=%.1f s  #%d",t.total_s,t.tick);
    if(alm.any_unacked){ImGui::SameLine(1200.f);
        float bl=fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f?1.f:0.f;
        ImGui::TextColored({1.f,bl,bl,1.f},"!! ALARM !!");}

    ImGui::End();ImGui::PopStyleColor();ImGui::PopStyleVar(2);
}

// ═══════════════════════════════════════════════════════════════════════════════
// LEFT CONTROL PANEL  (always visible)
// ═══════════════════════════════════════════════════════════════════════════════
static void RenderLeftPanel(ReactorState& s,SimTime& t,
                             bool& paused,bool& req_scram,float& speed,
                             FuelSystem& fuel,
                             TurbineSystem& turbines,
                             MoltenSaltSystem& salt,
                             HeliumCoolingSystem& helium,
                             AlarmSystem& alm)
{
    ImGuiIO& io=ImGui::GetIO();
    float lw=io.DisplaySize.x*LEFT_W_FRAC;
    float ch=io.DisplaySize.y-STATUS_H;
    if(!BeginTiled("##left",0.f,STATUS_H,lw,ch,ImGuiWindowFlags_NoTitleBar))
    {ImGui::End();return;}

    float pw=ImGui::GetContentRegionAvail().x;

    // Quick plasma summary
    Hdr("PLASMA");
    ImGui::TextColored(Col::GREEN,"P_fus: %.0f MW",s.fusion_power_MW);
    ImGui::TextColored(Col::AMBER,"Ti:  %.1f keV",s.plasma_temp_keV);
    ImGui::TextColored(Col::CYAN, "ne:  %.2e m-3",s.plasma_density_m3);

    Hdr("SIMULATION");
    ImGui::SetNextItemWidth(pw);
    ImGui::SliderFloat("##spd",&speed,0.1f,100.f,"Speed: %.1fx");
    // The pause button reflects *only* the manual `paused` flag now.
    // SCRAM no longer freezes the main loop (see main loop comment for
    // why) — instead, the SCRAM popup shows what's happening and the
    // operator can either wait for the plasma to ramp down or press
    // RESET.  When a SCRAM is active we still want a visual cue on the
    // button itself, so the label switches to "■ SCRAM ACTIVE" in red
    // while cmd_scram is true.
    if (s.cmd_scram) {
        ImGui::PushStyleColor(ImGuiCol_Button,        {0.5f,0.05f,0.05f,1.f});
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, {0.65f,0.08f,0.08f,1.f});
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  {0.75f,0.10f,0.10f,1.f});
        // Disabled-style button: clicking does nothing during SCRAM,
        // because the simulation is *not* actually paused — it's running
        // the shutdown sequence.  The operator interacts via the SCRAM
        // popup instead.
        ImGui::Button("■ SCRAM ACTIVE",{pw,0.f});
        ImGui::PopStyleColor(3);
    } else {
        if(ImGui::Button(paused?"▶ RESUME":"■■ PAUSE",{pw,0.f}))paused=!paused;
    }
    ImGui::TextColored(Col::GREY," SPACE=pause  F1=scram  ESC=quit");

    Hdr("REACTOR");
    bool cold=(s.plasma_status==PlasmaStatus::Cold);
    if(cold){
        // ── Startup pre-checks ──────────────────────────────────────────────
        // Before breakdown can occur, the confining toroidal field must
        // actually be established (magnets ramp up over tens of seconds),
        // the magnet system must not be in a quench state, the cryoplant
        // must be running (superconductors need to stay cold), there has
        // to be D/T fuel in the gas system to puff in, AND — new in
        // Round 4 — the vacuum vessel must be at base pressure (< 1e-3 Pa).
        bool field_ready = s.B_toroidal_T > 0.95f * s.sp_B_toroidal_T && s.sp_B_toroidal_T > 0.5f;
        bool cryo_ready  = s.cryo_ok && !s.quench_detected;
        bool fuel_ready  = (s.fuel_D_inventory_g > 1.0f) && (s.fuel_T_inventory_g > 1.0f);
        bool vac_ready   = s.vessel_vacuum_ok;  // vacuum vessel at base pressure
        bool ready       = field_ready && cryo_ready && fuel_ready && vac_ready;

        if(ready){
            if(GreenBtn("\xe2\x96\xb6 INITIATE PLASMA",{pw,26.f})){
                s.plasma_status=PlasmaStatus::Initiating;
                s.plasma_current_MA=0.5f;s.electron_temp_keV=3.f;
                s.plasma_density_m3=5e18f;s.mode=ReactorMode::Startup;
                // Clear the previous discharge's runaway-electron latches —
                // the strike counter (re_strike_count) survives as the
                // campaign damage record, but the alarm belongs to the shot
                // that caused it.
                s.re_wall_strike=false;s.re_beam_MA=0.f;s.alarm_runaway=false;
                // ── Auto-enable ECRH for breakdown assist ──────────────────────
                //  In a real tokamak, ECRH is used to pre-ionize the gas before
                //  the loop voltage is applied — this lowers the breakdown
                //  voltage and makes startup much more reliable.  We auto-
                //  enable ECRH at 5 MW (about 20% of its max) when the plasma
                //  is initiated, so the operator doesn't have to remember to
                //  do it manually.  The operator can turn it off later if
                //  they want to run purely ohmically.
                //
                //  This also fixes the "insanely difficult to keep stable"
                //  issue: without any aux heating, the plasma only had 5 MW
                //  (now 15 MW) of ohmic heating, which was barely enough to
                //  sustain the plasma.  ECRH at 5 MW gives a total of 20 MW
                //  of heating during Initiating, which is plenty to reach
                //  fusion-relevant temperatures.
                s.hcd_ecrh_on = true;
                s.hcd_ecrh_setpoint_MW = 5.0f;
            }
        } else {
            ImGui::BeginDisabled();
            ImGui::Button("\xe2\x96\xb6 INITIATE PLASMA",{pw,26.f});
            ImGui::EndDisabled();
            if(!field_ready) ImGui::TextColored(Col::AMBER,"  Waiting on toroidal field (B_T)");
            if(!cryo_ready)  ImGui::TextColored(Col::AMBER,"  Waiting on magnet cryoplant");
            if(!fuel_ready)  ImGui::TextColored(Col::AMBER,"  Waiting on D/T fuel inventory");
            if(!vac_ready)   ImGui::TextColored(Col::AMBER,"  Waiting on vessel vacuum (< 1e-3 Pa)");
        }
    } else {
        if(RedBtn("\xe2\x96\xa0 CONTROLLED SHUTDOWN",{pw,24.f})){
            s.mode=ReactorMode::Rampdown;
            // Note: do NOT zero setpoints here — the new Rampdown mode in
            // ControlSystem::runRampdown() handles the gentle ramp-down of
            // ALL setpoints (I_p, T_e, ne, P_aux) at safe rates.  Zeroing
            // them here would re-introduce the "fusion power drops to 0 at
            // ~30 MW" bug (sudden setpoint change → power balance limit cycle).
        }
    }
    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button,{0.5f,0.05f,0.05f,1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered,{0.75f,0.08f,0.08f,1.f});
    ImGui::PushStyleColor(ImGuiCol_ButtonActive,{1.f,0.12f,0.12f,1.f});
    if(ImGui::Button("⚡ SCRAM",{pw,30.f}))req_scram=true;
    ImGui::PopStyleColor(3);

    Hdr("SETPOINTS");
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##ip",&s.sp_plasma_current_MA,0.f,20.f,"Ip: %.1f MA");
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##te",&s.sp_electron_temp_keV,0.f,40.f,"Te: %.1f keV");
    float ne20=(float)(s.sp_density_m3*1e-20);
    ImGui::SetNextItemWidth(pw);if(ImGui::SliderFloat("##ne",&ne20,0.f,2.5f,"ne: %.2f e20"))s.sp_density_m3=ne20*1e20f;
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##bt",&s.sp_B_toroidal_T,0.f,12.f,"B_T: %.2f T");
    // ── Round 4: plasma-shape setpoints ────────────────────────────────────
    //  Elongation (κ) and triangularity (δ) are now operator-adjustable.
    //  These feed the shape controller in ControlSystem::runShapeControl,
    //  which drives the PF coils.  Higher κ = better confinement but lower
    //  vertical stability; ITER targets κ ≈ 1.7.
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##kappa",&s.sp_kappa,1.0f,2.2f,"kappa: %.2f");
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##delta",&s.sp_delta,0.0f,0.6f,"delta: %.2f");

    Hdr("FUEL");
    float pct = s.sp_fuel_rate * 100.f;
    if (ImGui::SliderFloat("##gp", &pct, 0.f, 100.f, "Gas: %.0f%%"))
        s.sp_fuel_rate = pct / 100.f;
    // Pellet injector doubles as the ELM pacing actuator — range extended to
    // 20 Hz so the operator can outrun the natural Type-I ELM frequency.
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##pf",&s.pellet_frequency_Hz,0.f,20.f,"Pellets: %.1f Hz");
    float hw=(pw-4.f)*0.5f;
    if(ImGui::Button("+D2",{hw,0.f}))fuel.resupplyDeuterium(100.f);ImGui::SameLine();
    if(ImGui::Button("+T",{hw,0.f})) fuel.resupplyTritium(25.f);

    Hdr("HELIUM SYSTEM");
    if(GreenBtn("Start Cryo",{hw,0.f}))helium.startCryoplant();ImGui::SameLine();
    if(RedBtn("Stop Cryo",{hw,0.f}))  helium.stopCryoplant();
    if(GreenBtn("Start He Cooling",{pw,0.f}))helium.startReactorCooling();
    if(ImGui::Button("Pump Cryostat",{pw,0.f}))helium.startCryostatPumping();

    Hdr("ALARMS");
    int na=0,nu=0;for(auto&a:alm.log){if(a.active)na++;if(!a.acked)nu++;}
    ImGui::TextColored(na?Col::RED:Col::GREEN,"Active:%d Unacked:%d",na,nu);
    float bw=(pw-4.f)*0.5f;
    if(ImGui::SmallButton("ACK"))alm.ackAll();ImGui::SameLine();
    if(ImGui::SmallButton("CLR"))alm.clearOld();

    Hdr("TIME");
    ImGui::TextColored(Col::GREY," Sim time");ImGui::SameLine();
    ImGui::TextColored(Col::GREEN,"%.2f s",t.total_s);
    ImGui::TextColored(Col::GREY," Tick    ");ImGui::SameLine();
    ImGui::TextColored(Col::GREEN,"%d",t.tick);

    ImGui::End();
}

// ═══════════════════════════════════════════════════════════════════════════════
// OVERVIEW TAB  (plasma + alarms in one pane)
// ═══════════════════════════════════════════════════════════════════════════════
static void TabOverview(const ReactorState& s,
                         ScrollBuf& h_pfus,ScrollBuf& h_te,ScrollBuf& h_ti,
                         ScrollBuf& h_ne,ScrollBuf& h_q,
                         AlarmSystem& alm,const SimTime& t)
{
    float pw=ImGui::GetContentRegionAvail().x;
    if(ImGui::BeginChild("##ov_left",{pw*0.55f,0.f},false)){
        Hdr("CORE PLASMA");
        ImVec4 qc=(s.q_safety<2.f)?Col::RED:(s.q_safety<2.5f)?Col::AMBER:Col::GREEN;
        ImVec4 bc=(s.beta_N>2.5f)?Col::RED:(s.beta_N>2.0f)?Col::AMBER:Col::GREEN;
        Row("Plasma Current",   s.plasma_current_MA,   "%.2f","MA");
        Row("Electron Temp",    s.electron_temp_keV,   "%.1f","keV",Col::AMBER);
        Row("Ion Temp",         s.plasma_temp_keV,     "%.1f","keV",Col::AMBER);
        Row("Density",          s.plasma_density_m3,   "%.2e","m-3",Col::CYAN);
        ImVec4 gwc=(s.greenwald_frac>1.0f)?Col::RED:(s.greenwald_frac>0.85f)?Col::AMBER:Col::GREEN;
        Row("Greenwald Fraction",s.greenwald_frac,     "%.2f","",  gwc);
        Row("Safety Factor q95",s.q_safety,            "%.3f","",  qc);
        Row("Beta (normalized)",s.beta*100.f,          "%.3f","%");
        Row("Beta_N (Troyon)",  s.beta_N,              "%.2f","",  bc);
        Row("tau_E",            s.tau_E_s,             "%.2f","s", Col::CYAN);

        // ── Confinement regime (L/H) ─────────────────────────────────────────
        ImGui::TextColored(Col::GREY,"%-22s","Confinement Regime");ImGui::SameLine();
        ImGui::TextColored(s.h_mode?Col::GREEN:Col::AMBER,"%-12s",
                           s.h_mode?"H-MODE":"L-MODE");
        ImGui::SameLine();ImGui::TextColored(Col::GREY,"H98=%.2f",s.H98_factor);
        Row("L-H Threshold",    s.P_LH_MW,             "%.1f","MW",Col::CYAN);
        float P_heat=s.alpha_power_MW+s.P_aux_actual_MW+s.P_ohm_MW;
        float lh_frac=(s.P_LH_MW>0.1f)?P_heat/s.P_LH_MW:0.f;
        ImGui::TextColored(Col::GREY," P_heat/P_LH");ImGui::SameLine(150.f);
        Bar(std::min(lh_frac*0.5f,1.f),pw*0.30f,
            s.h_mode?Col::GREEN_DIM:Col::AMBER,Col::GREEN,0.5f);
        ImGui::SameLine();ImGui::TextColored(s.h_mode?Col::GREEN:Col::AMBER,"%.2f",lh_frac);

        // ── Sawtooth activity ────────────────────────────────────────────────
        if(s.sawtooth_period_s>0.f){
            ImGui::TextColored(Col::GREY,"%-22s","Sawteeth (q0<1)");ImGui::SameLine();
            ImGui::TextColored(Col::AMBER,"ACTIVE  q0=%.2f  T=%.1fs  n=%d",
                               s.q0_estimate,s.sawtooth_period_s,s.sawtooth_count);
        } else {
            Row("q0 (on-axis)", s.q0_estimate,         "%.2f","",
                s.q0_estimate<1.f?Col::AMBER:Col::GREEN);
        }

        // ── ELM activity (H-mode edge crashes) ───────────────────────────────
        //  Amber when unmitigated Type-I ELMs exceed ~2 MJ per crash (tile
        //  erosion territory); green when pellet pacing keeps them small.
        if(s.elm_freq_Hz>0.f){
            ImGui::TextColored(Col::GREY,"%-22s","ELMs (Type-I)");ImGui::SameLine();
            ImVec4 ec = s.elm_paced ? Col::GREEN
                      : (s.elm_size_MJ>2.f ? Col::AMBER : Col::GREY);
            ImGui::TextColored(ec,"%.1f Hz  dW=%.1f MJ  n=%d%s",
                               s.elm_freq_Hz,s.elm_size_MJ,s.elm_count,
                               s.elm_paced?"  [PACED]":"");
            if(!s.elm_paced && s.elm_size_MJ>2.f)
                ImGui::TextColored(Col::AMBER,
                    "  raise pellet freq above %.1f Hz to pace ELMs",s.elm_freq_Hz);
        }

        // ── Runaway electron beam (disruption aftermath) ─────────────────────
        if(s.re_beam_MA>0.1f){
            ImGui::TextColored(Col::RED,"%-22s%.2f MA  ⚡ RUNAWAY BEAM",
                               "RE Beam",s.re_beam_MA);
        }
        if(s.re_wall_strike){
            ImGui::TextColored(Col::RED,"  RE WALL STRIKE (%d this campaign) — check first wall",
                               s.re_strike_count);
        }

        Row("He Ash",           s.helium_fraction*100.f,"%.2f","%");
        Row("Impurities",       s.impurity_fraction*100.f,"%.2f","%",
            s.impurity_fraction>0.05f?Col::AMBER:Col::GREY);
        // ── Round 4: plasma shape ────────────────────────────────────────────
        Row("Elongation (kappa)",s.kappa,              "%.2f","",  Col::CYAN);
        Row("Triangularity (delta)",s.delta,           "%.2f","",  Col::CYAN);
        Row("Volume",           s.plasma_volume_m3,    "%.0f","m3",Col::GREY);

        Hdr("IN-VESSEL INVENTORY");
        Row("Deuterium",  s.vessel_D_g,        "%.4f","g",Col::CYAN);
        Row("Tritium",    s.vessel_T_g,        "%.4f","g",Col::AMBER);
        Row("He-4 Ash",   s.vessel_He_g,       "%.4f","g",Col::GREY);
        Row("Junk Mass",  s.vessel_impurity_g, "%.4f","g",
            s.vessel_impurity_g>1.f?Col::RED:Col::GREY);
        Row("Exhaust Pump Speed", s.exhaust_pumping_Ls, "%.1f","L/s");

        Hdr("POWER BALANCE");
        Row("Fusion Power",   s.fusion_power_MW,  "%.1f","MW");
        Row("Alpha Heating",  s.alpha_power_MW,   "%.1f","MW",Col::AMBER);
        Row("Aux Heating",    s.P_aux_actual_MW,  "%.1f","MW",Col::AMBER);
        Row("Ohmic Heating",  s.P_ohm_MW,         "%.2f","MW",Col::AMBER);
        Row("Loop Voltage",   s.loop_voltage_V,   "%.2f","V", Col::CYAN);
        Row("Bremsstrahlung", s.P_brem_MW,        "%.1f","MW",Col::GREY);
        Row("Synchrotron",    s.P_sync_MW,        "%.2f","MW",Col::GREY);
        Row("Line Radiation", s.P_line_MW,        "%.1f","MW",
            s.P_line_MW>0.5f*(s.alpha_power_MW+s.P_aux_actual_MW+s.P_ohm_MW+0.1f)
            ?Col::RED:Col::GREY);
        Row("Conduction",     s.P_cond_MW,        "%.1f","MW",Col::GREY);
        Row("dW/dt",          s.dW_dt_MW,         "%+.1f","MW",
            s.dW_dt_MW>=0.f?Col::GREEN:Col::AMBER);
        Row("Q scientific",   s.Q_scientific,     "%.3f","");
        Row("Blanket Heat",   s.blanket_heat_MW,  "%.1f","MW",Col::AMBER);
        Row("Net Electric",   s.net_electric_MW,  "%.1f","MW",Col::YELLOW);
        Row("Grid Frequency", s.grid_frequency_Hz,"%.3f","Hz",
            s.grid_frequency_Hz<49.f?Col::AMBER:Col::GREEN);

        Hdr("MAGNETS");
        Row("B_T",          s.B_toroidal_T,    "%.3f","T");
        Row("Coil Current", s.coil_current_kA, "%.1f","kA");
        Row("CS Current",   s.cs_current_kA,   "%.1f","kA", Col::CYAN);
        // CS volt-second budget — the pulse-length clock of the discharge.
        Row("CS Flux Left", s.cs_flux_remaining_Wb, "%.1f","Wb",
            s.cs_flux_exhausted?Col::RED
            :(s.cs_flux_remaining_Wb<20.f?Col::AMBER:Col::GREEN));
        if(s.cs_flux_exhausted)
            ImGui::TextColored(Col::RED,"  CS EXHAUSTED — non-inductive drive only");
        Row("Magnet Temp",  s.magnet_temp_K,   "%.2f","K",
            s.magnet_temp_K>6.f?Col::AMBER:Col::GREEN);
        Lamp(" QUENCH ",s.quench_detected);

        Hdr("VACUUM & TRITIUM");
        Row("Vessel Pressure", s.vessel_pressure_Pa, "%.2e", "Pa",
            s.vessel_vacuum_ok ? Col::GREEN : (s.vessel_pressure_Pa > 1.f ? Col::RED : Col::AMBER));
        Row("T in Plant",      s.tritium_in_plant_g, "%.1f", "g", Col::AMBER);
        Row("TBR",             s.tbr_current,         "%.3f", "",
            s.tbr_current >= 1.0f ? Col::GREEN : Col::AMBER);

        Hdr("FIRST WALL");
        Row("Wall Temp",     s.first_wall_temp_K,        "%.0f","K",
            s.first_wall_temp_K>1200.f?Col::RED
            :(s.first_wall_temp_K>800.f?Col::AMBER:Col::GREEN));
        Row("Neutron Load",  s.neutron_wall_load_MW_m2,  "%.3f","MW/m2",Col::CYAN);
        // Displacement damage — the campaign-lifetime wear meter.  ITER's
        // first wall is rated for ~3 dpa before replacement.
        Row("Damage",        s.fw_dpa,                   "%.4f","dpa",
            s.fw_dpa>3.f?Col::RED:(s.fw_dpa>1.f?Col::AMBER:Col::GREY));

        Hdr("PLASMA TRENDS");
        ImGui::TextColored(Col::GREY," Fusion Power [MW]");Plot("##op",h_pfus,pw*.55f,40.f,Col::GREEN);
        ImGui::TextColored(Col::GREY," Te [keV]");         Plot("##ote",h_te,pw*.55f,40.f,Col::AMBER);
        ImGui::TextColored(Col::GREY," Ti [keV]  (watch for sawteeth)");
                                                            Plot("##oti",h_ti,pw*.55f,40.f,Col::ORANGE);
        ImGui::TextColored(Col::GREY," ne [e20 m-3]");     Plot("##one",h_ne,pw*.55f,40.f,Col::CYAN);
        ImGui::TextColored(Col::GREY," q95");              Plot("##oq",h_q, pw*.55f,40.f,qc);
    }
    ImGui::EndChild();
    ImGui::SameLine();

    if(ImGui::BeginChild("##ov_right",{0.f,0.f},false)){
        Hdr("SALT SYSTEM");
        Row("Hot Tank Temp",  s.hot_tank_temp_K,   "%.0f","K",Col::ORANGE);
        Row("Cold Tank Temp", s.cold_tank_temp_K,  "%.0f","K",Col::CYAN);
        Row("Hot Level",      s.hot_tank_level_m,  "%.1f","m");
        Row("Salt Flow",      s.salt_flow_total_kg_s,"%.0f","kg/s");

        Hdr("HELIUM");
        Row("Magnet He",     s.magnet_he_temp_K,   "%.2f","K",
            s.magnet_he_temp_K>5.5f?Col::RED:Col::GREEN);
        Row("Cryo OK",       s.cryo_ok?1.f:0.f,   s.cryo_ok?"YES":"NO","",
            s.cryo_ok?Col::GREEN:Col::RED);
        Row("He Outlet",     s.reactor_he_outlet_K, "%.0f","K",Col::AMBER);

        Hdr("ALARMS");
        float bl=fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f?1.f:0.f;
        bool any=false;
        for(auto& a:alm.log){
            if(!a.active)continue;
            any=true;
            // Critical alarms blink red, warnings are amber, info is cyan
            ImVec4 ac = (a.severity == AlarmSeverity::Critical && !a.acked)
                       ? ImVec4{1.f, bl, bl, 1.f}
                       : (a.severity == AlarmSeverity::Critical ? Col::RED
                          : (a.severity == AlarmSeverity::Warning ? Col::AMBER : Col::CYAN));
            const char* icon = (a.severity == AlarmSeverity::Critical) ? "⛔"
                              : (a.severity == AlarmSeverity::Warning)  ? "⚠"
                              : "ℹ";
            ImGui::TextColored(ac," %s [t=%.1fs] %s%s",
                icon, a.t, a.msg.c_str(), a.acked?" (ack)":"");
            // Show cause on the next line (truncated to fit)
            ImGui::TextColored(Col::WHITE, "      cause: %s",
                a.cause.length() > 90 ? (a.cause.substr(0,87) + "...").c_str()
                                      : a.cause.c_str());
        }
        if(!any)ImGui::TextColored(Col::GREEN,"  No active alarms");
        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,"  See ALARMS tab for full details + remediation");
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// OP SPACE TAB — live POPCON operating-space diagram
// ═══════════════════════════════════════════════════════════════════════════════
//
//  The classic reactor-design plot: density (y) vs ion temperature (x), with
//  the required auxiliary power P_aux(n,T) shaded at every point, the
//  operating limits overlaid (Greenwald density, Troyon beta, L-H power
//  threshold), and the LIVE operating point tracing a breadcrumb trail as
//  the operator drives the discharge.  Computed with exactly the same
//  ConfinementPhysics functions as the simulation itself, so what you see
//  is what the plasma feels.
//
//    GREEN region   : ignition (P_aux ≤ 0 — the plasma self-heats)
//    BLUE → AMBER   : increasing auxiliary power required for steady state
//    RED line       : Greenwald density limit
//    MAGENTA curve  : Troyon β_N = 2.5 limit
//    CYAN dots      : L→H transition boundary (Martin 2008)
//
static void TabOpSpace(const ReactorState& s)
{
    using namespace ConfinementPhysics;

    constexpr int   NX = 56, NY = 44;
    constexpr float T_MIN = 1.f,    T_MAX = 40.f;     // ion temp [keV]
    constexpr float N_MIN = 0.05f,  N_MAX = 2.5f;     // density [1e20 m-3]

    // ── Cached P_aux(n,T) map — recomputed when the machine point moves ────
    static float paux_map[NY][NX];
    static float lh_map  [NY][NX];      // P_heat − P_LH (sign → boundary)
    static float q10_map [NY][NX];      // P_fus − 10·P_aux (Q=10 contour)
    static float c_Ip=-1.f, c_B=-1.f, c_k=-1.f, c_d=-1.f, c_Zeff=-1.f, c_imp=-1.f;
    static int   frame_cnt=0;

    TokamakGeometry g;
    g.R_major_m=s.R_out_m; g.a_minor_m=2.0f;
    g.kappa=std::clamp(s.kappa,1.0f,2.5f);
    g.B_toroidal_T=std::max(s.B_toroidal_T,0.01f);
    g.I_plasma_MA=std::max(s.plasma_current_MA,0.5f);
    float Z_eff = 1.0f + 4.0f*s.helium_fraction
                + (s.impurity_fraction>0?10.0f*s.impurity_fraction:0.0f);
    float f_DT  = std::max(1.0f-s.helium_fraction-s.impurity_fraction,0.0f);

    bool stale = (fabsf(c_Ip-g.I_plasma_MA)>0.05f)||(fabsf(c_B-g.B_toroidal_T)>0.02f)
               ||(fabsf(c_k-g.kappa)>0.01f)||(fabsf(c_d-s.delta)>0.01f)
               ||(fabsf(c_Zeff-Z_eff)>0.02f)||(fabsf(c_imp-s.impurity_fraction)>0.002f)
               ||(++frame_cnt>240);
    if(stale){
        frame_cnt=0; c_Ip=g.I_plasma_MA; c_B=g.B_toroidal_T; c_k=g.kappa;
        c_d=s.delta; c_Zeff=Z_eff; c_imp=s.impurity_fraction;
        float V=plasmaVolume(g), S=plasmaSurfaceArea(g);
        for(int iy=0;iy<NY;iy++){
            float n20=N_MIN+(N_MAX-N_MIN)*(iy+0.5f)/NY;
            float ne=n20*1e20f;
            float nD=0.5f*ne*f_DT, nT=0.5f*ne*f_DT;     // 50:50 D-T map
            float P_LH=martinLHThreshold(ne,g.B_toroidal_T,S);
            for(int ix=0;ix<NX;ix++){
                float T=T_MIN+(T_MAX-T_MIN)*(ix+0.5f)/NX;   // assume Te=Ti on map
                float W_MJ=1.5f*ne*(2.f*T)*1.602176634e-16f*V*1e-6f;
                //  Same profile-peaking factors as solvePowerBalance so the
                //  map matches what the live plasma actually experiences.
                float P_fus=fusionPowerDensity(nD,nT,T)*2.5f*V*1e-6f;
                float P_a  =P_fus*(3.521f/17.589f);
                float P_ohm=ohmicHeating(g,g.I_plasma_MA,T,Z_eff).P_ohm_MW;
                float P_rad=(bremSStrahlung(ne,T,Z_eff)*1.4f
                            +synchrotron(ne,T,g.B_toroidal_T)
                            +lineRadiationDensity(ne,T,s.impurity_fraction,
                                                  0.001f*s.impurity_fraction))*V*1e-6f;
                // Fixed-point iterate P_aux: τ_E depends on total heating
                float P_aux=10.f;
                for(int it=0;it<5;it++){
                    float P_loss=std::max(P_a+P_ohm+std::max(P_aux,0.f),0.5f);
                    float tau=ipb98y2(g,P_loss,ne*1e-19f,2.5f,1.0f); // H-mode map
                    float P_cond=(tau>1e-4f)?W_MJ/tau:1e6f;
                    P_aux=P_cond+P_rad-P_a-P_ohm;
                }
                paux_map[iy][ix]=P_aux;
                lh_map  [iy][ix]=(P_a+P_ohm+std::max(P_aux,0.f))-P_LH;
                q10_map [iy][ix]=P_fus-10.f*std::max(P_aux,0.01f);
            }
        }
    }

    // ── Layout: canvas left, regime/ledger panel right ──────────────────────
    float pw=ImGui::GetContentRegionAvail().x;
    float ph=ImGui::GetContentRegionAvail().y;

    if(ImGui::BeginChild("##ops_canvas",{pw*0.62f,0.f},false)){
        Hdr("OPERATING SPACE  (POPCON)");
        ImGui::TextColored(Col::GREY,
            " Steady-state P_aux required at each (T_i, n_e) for the CURRENT machine point");
        ImGui::TextColored(Col::GREY,
            " I_p=%.1f MA  B_T=%.2f T  kappa=%.2f  delta=%.2f  Z_eff=%.2f",
            g.I_plasma_MA,g.B_toroidal_T,g.kappa,s.delta,Z_eff);

        ImVec2 p0=ImGui::GetCursorScreenPos();
        ImVec2 avail=ImGui::GetContentRegionAvail();
        float ML=52.f, MB=30.f, MT=8.f, MR=14.f;
        float cw=std::max(avail.x-ML-MR,120.f);
        float chh=std::max(avail.y-MT-MB-4.f,120.f);
        ImVec2 o={p0.x+ML,p0.y+MT};                  // plot origin (top-left)
        auto X=[&](float T){return o.x+cw*(T-T_MIN)/(T_MAX-T_MIN);};
        auto Y=[&](float n20){return o.y+chh*(1.f-(n20-N_MIN)/(N_MAX-N_MIN));};
        ImDrawList* dl=ImGui::GetWindowDrawList();

        // Heatmap cells
        float dx=cw/NX, dy=chh/NY;
        for(int iy=0;iy<NY;iy++)for(int ix=0;ix<NX;ix++){
            float pa=paux_map[iy][ix];
            ImU32 col;
            if(pa<=0.f){                                   // ignition
                float d=std::clamp(-pa/100.f,0.f,1.f);
                col=ImGui::ColorConvertFloat4ToU32({0.05f+0.15f*d,0.45f+0.45f*d,0.12f+0.2f*d,1.f});
            } else if(pa<150.f){
                float d=pa/150.f;
                col=ImGui::ColorConvertFloat4ToU32({0.06f+0.5f*d,0.10f+0.30f*d,0.28f*(1.f-d)+0.04f,1.f});
            } else {
                col=ImGui::ColorConvertFloat4ToU32({0.10f,0.05f,0.05f,1.f});
            }
            float x=o.x+ix*dx, y=o.y+chh-(iy+1)*dy;
            dl->AddRectFilled({x,y},{x+dx+0.5f,y+dy+0.5f},col);
        }

        // Q=10 contour (yellow dots at sign changes) & L-H boundary (cyan)
        for(int iy=0;iy<NY;iy++)for(int ix=1;ix<NX;ix++){
            float xm=o.x+ix*dx, ym=o.y+chh-(iy+0.5f)*dy;
            if(q10_map[iy][ix-1]*q10_map[iy][ix]<0.f)
                dl->AddCircleFilled({xm,ym},1.6f,IM_COL32(230,230,60,255));
            if(lh_map[iy][ix-1]*lh_map[iy][ix]<0.f)
                dl->AddCircleFilled({xm,ym},1.6f,IM_COL32(60,220,220,255));
        }
        for(int ix=0;ix<NX;ix++)for(int iy=1;iy<NY;iy++){
            float xm=o.x+(ix+0.5f)*dx, ym=o.y+chh-iy*dy;
            if(q10_map[iy-1][ix]*q10_map[iy][ix]<0.f)
                dl->AddCircleFilled({xm,ym},1.6f,IM_COL32(230,230,60,255));
            if(lh_map[iy-1][ix]*lh_map[iy][ix]<0.f)
                dl->AddCircleFilled({xm,ym},1.6f,IM_COL32(60,220,220,255));
        }

        // Greenwald limit (horizontal red line)
        float nG=g.I_plasma_MA/(3.14159265f*g.a_minor_m*g.a_minor_m);
        if(nG>N_MIN&&nG<N_MAX)
            dl->AddLine({o.x,Y(nG)},{o.x+cw,Y(nG)},IM_COL32(235,50,50,255),2.f);

        // Troyon beta_N = 2.5 curve: n20(T) with T_sum = 2T
        {
            constexpr float MU0=1.256637e-6f;
            ImVec2 prev{-1,-1};
            for(int ix=0;ix<=64;ix++){
                float T=T_MIN+(T_MAX-T_MIN)*ix/64.f;
                float beta_lim=2.5f*g.I_plasma_MA/(100.f*g.a_minor_m*g.B_toroidal_T);
                float n=beta_lim*g.B_toroidal_T*g.B_toroidal_T
                        /(MU0*(2.f*T)*1.602e-16f)*1e-20f;
                if(n>N_MIN&&n<N_MAX){
                    ImVec2 pt{X(T),Y(n)};
                    if(prev.x>0)dl->AddLine(prev,pt,IM_COL32(220,70,220,255),2.f);
                    prev=pt;
                }else prev={-1,-1};
            }
        }

        // Breadcrumb trail of the live operating point (sampled ~4 Hz)
        static std::deque<ImVec2> trail;   // data coords (T_i, n20)
        static double last_t=-1.0;
        if(s.time_s<last_t) last_t=-1.0;   // sim clock reset (COLD RESTART)
        if(s.time_s-last_t>0.25){
            last_t=s.time_s;
            if(s.plasma_status!=PlasmaStatus::Cold
             &&s.plasma_status!=PlasmaStatus::Quenched){
                //  Clamp into the plotted range so a cold/edge operating
                //  point renders ON the frame border rather than outside
                //  it (an unclamped just-initiated plasma at T≈0, n≈0
                //  mapped to a spot beyond the bottom-left corner).
                float Tc=std::clamp(s.ion_temp_keV,T_MIN,T_MAX);
                float nc=std::clamp((float)(s.plasma_density_m3*1e-20),N_MIN,N_MAX);
                trail.push_back({Tc,nc});
            } else trail.clear();
            while(trail.size()>240)trail.pop_front();
        }
        for(size_t i=1;i<trail.size();i++){
            float a=0.15f+0.85f*(float)i/trail.size();
            dl->AddLine({X(trail[i-1].x),Y(trail[i-1].y)},
                        {X(trail[i].x),  Y(trail[i].y)},
                        ImGui::ColorConvertFloat4ToU32({0.95f,0.95f,0.95f,a}),1.5f);
        }
        if(!trail.empty()){
            ImVec2 c{X(trail.back().x),Y(trail.back().y)};
            dl->AddCircle(c,6.f,IM_COL32(255,255,255,255),12,2.f);
            dl->AddLine({c.x-10,c.y},{c.x+10,c.y},IM_COL32(255,255,255,180),1.f);
            dl->AddLine({c.x,c.y-10},{c.x,c.y+10},IM_COL32(255,255,255,180),1.f);
        }

        // Axes + frame + tick labels
        dl->AddRect({o.x,o.y},{o.x+cw,o.y+chh},IM_COL32(70,110,70,255));
        char lbl[32];
        for(int T=5;T<=40;T+=5){
            snprintf(lbl,32,"%d",T);
            dl->AddText({X((float)T)-6,o.y+chh+4},IM_COL32(140,160,140,255),lbl);
        }
        dl->AddText({o.x+cw*0.5f-40,o.y+chh+16},IM_COL32(140,160,140,255),"T_i [keV]");
        for(float n=0.5f;n<N_MAX;n+=0.5f){
            snprintf(lbl,32,"%.1f",n);
            dl->AddText({o.x-34,Y(n)-7},IM_COL32(140,160,140,255),lbl);
        }
        dl->AddText({p0.x,o.y-2},IM_COL32(140,160,140,255),"n_e");
        dl->AddText({p0.x,o.y+12},IM_COL32(140,160,140,255),"e20");

        ImGui::Dummy({avail.x,chh+MB+MT});
        ImGui::TextColored(Col::GREEN, " GREEN = ignition");ImGui::SameLine();
        ImGui::TextColored({0.9f,0.9f,0.3f,1.f}," YELLOW = Q=10");ImGui::SameLine();
        ImGui::TextColored(Col::CYAN,  " CYAN = L-H");ImGui::SameLine();
        ImGui::TextColored(Col::RED,   " RED = Greenwald");ImGui::SameLine();
        ImGui::TextColored({0.86f,0.3f,0.86f,1.f}," MAGENTA = Troyon");
    }
    ImGui::EndChild();
    ImGui::SameLine();

    // ── Right column: regime + live power ledger ────────────────────────────
    if(ImGui::BeginChild("##ops_right",{0.f,0.f},false)){
        (void)ph;
        Hdr("CONFINEMENT REGIME");
        ImGui::TextColored(s.h_mode?Col::GREEN:Col::AMBER,"  %s",
                           s.h_mode?"H-MODE (H98=1.00)":"L-MODE (H98=0.55)");
        Row("P_LH (Martin 08)",s.P_LH_MW,"%.1f","MW",Col::CYAN);
        float P_heat=s.alpha_power_MW+s.P_aux_actual_MW+s.P_ohm_MW;
        Row("P_heat total",    P_heat,   "%.1f","MW",Col::AMBER);
        ImGui::TextColored(Col::GREY,
            s.h_mode?"  Back-transition below %.1f MW":"  Transition above %.1f MW",
            s.h_mode?0.5f*s.P_LH_MW:s.P_LH_MW);

        Hdr("POWER LEDGER [MW]");
        float rw=ImGui::GetContentRegionAvail().x;
        float maxP=std::max({s.alpha_power_MW,s.P_aux_actual_MW,s.P_ohm_MW,
                             s.P_brem_MW,s.P_sync_MW,s.P_line_MW,s.P_cond_MW,10.f});
        auto Ledger=[&](const char* n,float v,ImVec4 c){
            ImGui::TextColored(Col::GREY,"%-10s",n);ImGui::SameLine(84.f);
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram,c);
            ImGui::ProgressBar(std::clamp(v/maxP,0.f,1.f),{rw*0.52f,11.f},"");
            ImGui::PopStyleColor();
            ImGui::SameLine();ImGui::TextColored(c,"%6.1f",v);
        };
        ImGui::TextColored(Col::GREEN," SOURCES");
        Ledger("alpha",   s.alpha_power_MW,  Col::GREEN);
        Ledger("aux H&CD",s.P_aux_actual_MW, Col::GREEN_DIM);
        Ledger("ohmic",   s.P_ohm_MW,        Col::GREEN_DIM);
        ImGui::TextColored(Col::RED," SINKS");
        Ledger("brems",   s.P_brem_MW, Col::ORANGE);
        Ledger("synchro", s.P_sync_MW, Col::ORANGE);
        Ledger("line rad",s.P_line_MW, s.P_line_MW>P_heat?Col::RED:Col::ORANGE);
        Ledger("conduct", s.P_cond_MW, Col::AMBER);
        ImGui::Spacing();
        Row("dW/dt", s.dW_dt_MW, "%+.1f","MW",
            s.dW_dt_MW>=0.f?Col::GREEN:Col::AMBER);
        Row("Stored W", s.stored_energy_GJ*1000.f, "%.0f","MJ",Col::CYAN);

        Hdr("OPERATING LIMITS");
        auto LimitBar=[&](const char* n,float frac,const char* txt){
            ImGui::TextColored(Col::GREY,"%-12s",n);ImGui::SameLine(96.f);
            ImVec4 c=frac>1.f?Col::RED:frac>0.85f?Col::AMBER:Col::GREEN;
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram,c);
            ImGui::ProgressBar(std::clamp(frac,0.f,1.f),{rw*0.45f,11.f},"");
            ImGui::PopStyleColor();
            ImGui::SameLine();ImGui::TextColored(c,"%s",txt);
        };
        char b1[24],b2[24],b3[24];
        snprintf(b1,24,"%.2f/1.00",s.greenwald_frac);
        snprintf(b2,24,"%.2f/2.50",s.beta_N);
        float qfrac=(s.q_safety>0.01f)?2.0f/s.q_safety:0.f;
        snprintf(b3,24,"q95=%.2f",s.q_safety);
        LimitBar("Greenwald",s.greenwald_frac,b1);
        LimitBar("Troyon bN",s.beta_N/2.5f,b2);
        LimitBar("Kink (q>2)",qfrac,b3);

        Hdr("SAWTOOTH / CORE MHD");
        Row("q0 estimate",s.q0_estimate,"%.2f","",
            s.q0_estimate<1.f?Col::AMBER:Col::GREEN);
        if(s.sawtooth_period_s>0.f){
            ImGui::TextColored(Col::AMBER,"  SAWTEETH ACTIVE  period %.1f s",
                               s.sawtooth_period_s);
            Row("Crash count",(float)s.sawtooth_count,"%.0f","");
        } else {
            ImGui::TextColored(Col::GREY,"  quiescent core (q0 > 1)");
        }
        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,
            " Tip: raise I_p to lift the Greenwald line and\n"
            " widen the map; raise B_T to relax Troyon; pay\n"
            " the L-H toll (cyan) before pushing for burn.");
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// TURBINE TAB  (one turbine unit's full panel)
// ═══════════════════════════════════════════════════════════════════════════════
static void TabTurbine(TurbineUnitController& ctrl,
                        ElectricalGridSystem& egrid,
                        int idx)
{
    TurbineUnit& u=ctrl.s;
    float pw=ImGui::GetContentRegionAvail().x;
    float col_w=pw*0.33f;

    // ── Column 1: steam generator + MSIV ─────────────────────────────────────
    if(ImGui::BeginChild("##tc1",{col_w,0.f},false)){
        ImVec4 stc=turbStateCol(u.state);
        Hdr("TURBINE STATUS");
        ImGui::TextColored(stc,"  %s",turbStateStr(u.state));
        Row("Rotor Speed",    u.rpm,              "%.1f","RPM",
            u.rpm>3100.f?Col::RED:(u.rpm>2900.f&&u.rpm<3100.f?Col::GREEN:Col::AMBER));
        Row("Shaft Power",    u.shaft_power_MW,   "%.1f","MW");
        Row("Generator MW",   u.generator.power_MW,"%.1f","MW",Col::YELLOW);
        Row("Gen Frequency",  u.generator.frequency_Hz,"%.3f","Hz",
            std::abs(u.generator.frequency_Hz-50.f)<0.1f?Col::GREEN:Col::AMBER);

        Hdr("STEAM GENERATOR");
        Row("SG Pressure",    u.sg_pressure_MPa,  "%.2f","MPa",
            u.sg_pressure_MPa>19.f?Col::RED:Col::GREEN);
        Row("Steam Temp",     u.sg_steam_temp_K,  "%.0f","K",Col::AMBER);
        Row("SG Level",       u.sg_level_m,       "%.1f","m",
            u.sg_level_m<3.f?Col::RED:Col::GREEN);
        Row("Steam Flow",     u.sg_steam_flow_kg_s,"%.0f","kg/s");
        Row("Heat Input",     u.sg_heat_input_MW, "%.0f","MW",Col::ORANGE);
        Lamp(" HI PRESSURE ",u.alarm_hi_sg_pressure);
        ImGui::SameLine(); Lamp(" LO LEVEL ",u.alarm_lo_sg_level);

        Hdr("FEEDWATER PUMPS");
        for(int i=0;i<2;i++){
            auto& p=u.fw_pump[i];
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY,"  FW Pump %d",i+1);
            ImGui::SameLine(100.f);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY,
                               p.running?"%-6s %-5.0f kg/s  %.2f MW":"STOPPED",
                               p.trip?"TRIP":">",p.flow_kg_s,p.power_MW);
            float bw=(pw*.33f-8.f)*.5f;
            ImGui::PushID(i*100+idx*10+1);
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f})){p.running=false;}
            ImGui::SetNextItemWidth(pw*.33f-8.f);
            float pcta = p.speed_frac * 100.f;
            if (ImGui::SliderFloat("##sp", &pcta, 0.f, 100.f, "Speed: %.0f%%"))
                p.speed_frac = pcta / 100.f;
            ImGui::PopID();
        }
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Column 2: condenser, hotwell, preheaters ──────────────────────────────
    if(ImGui::BeginChild("##tc2",{col_w,0.f},false)){
        Hdr("MSIV & STEAM PATH");
        ImGui::TextColored(u.msiv_open?Col::GREEN:Col::RED,"  MSIV: %s (pos %.0f%%)",
                           u.msiv_open?"OPEN":"CLOSED",u.msiv_position*100.f);
        ImGui::TextColored(Col::GREY,"  Setpoint:");ImGui::SameLine();
        ImGui::SetNextItemWidth(pw*.33f-80.f);
        float pctb = u.msiv_setpoint * 100.f;
        if (ImGui::SliderFloat("##msiv", &pctb, 0.f, 100.f, "%.0f%%"))
            u.msiv_setpoint = pctb / 100.f;

        float bw=(pw*.33f-8.f)*.5f;
        if(GreenBtn("Open MSIV",{bw,0.f})){u.msiv_open=true;u.msiv_trip_latch=false;}
        ImGui::SameLine();
        if(RedBtn("Close/Trip",{bw,0.f}))ctrl.cmdTrip();

        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,"  Bypass Valve:");ImGui::SameLine();
        ImGui::SetNextItemWidth(pw*.33f-100.f);
        float pctc = u.bypass_valve_pos * 100.f;
        if (ImGui::SliderFloat("##bvp", &pctc, 0.f, 100.f, "%.0f%%"))
            u.bypass_valve_pos = pctc / 100.f;
        ImGui::Text("  Relief Valve: %s  (setpt %.1f MPa)",
                    u.relief_valve_open?"OPEN":"CLOSED",u.relief_setpoint_MPa);
        ImGui::SetNextItemWidth(pw*.33f-8.f);
        ImGui::SliderFloat("##rv",&u.relief_setpoint_MPa,15.f,22.f,"Relief: %.1f MPa");

        Hdr("PREHEATERS");
        for(int i=0;i<4;i++){
            auto& ph=u.preheater[i];
            ImGui::PushID(i*200+idx*20+5);
            ImGui::TextColored(ph.enabled?Col::GREEN:Col::GREY," PH-%d",i+1);
            ImGui::SameLine(60.f);
            if(ph.enabled)
                ImGui::TextColored(Col::AMBER,"%.0f K  %.1f MW",ph.fw_outlet_temp_K,ph.heat_transferred_MW);
            else ImGui::TextColored(Col::GREY,"bypassed");
            ImGui::SameLine(pw*.33f-70.f);
            ImGui::Checkbox(ph.enabled?"Enabled##ph":"Disabled##ph",&ph.enabled);
            ImGui::PopID();
        }
        Row("FW Temp post-PH",u.fw_temp_after_ph_K,"%.0f","K",Col::AMBER);

        Hdr("CONDENSER");
        auto& c=u.condenser;
        Row("Cond Pressure",  c.pressure_kPa,"%.2f","kPa",
            c.pressure_kPa>12.f?Col::RED:(c.pressure_kPa>8.f?Col::AMBER:Col::GREEN));
        Row("Cond Temp",      c.temp_K,      "%.1f","K");
        Row("Air Fraction",   c.air_fraction*100.f,"%.3f","%",
            c.air_fraction>0.01f?Col::AMBER:Col::GREEN);
        Row("Condensate Flow",c.condensate_flow_kg_s,"%.0f","kg/s");
        Lamp(" LO VAC ",u.alarm_lo_condenser_vac);

        ImGui::Spacing();ImGui::Text("  CAR Pumps:");
        for(int i=0;i<4;i++){
            ImGui::SameLine();
            ImGui::PushID(i*300+idx*30);
            bool on=c.car_pump[i];
            ImGui::PushStyleColor(ImGuiCol_Button,on?Col::GREEN:ImVec4{0.2f,0.2f,0.2f,1.f});
            char lbl[8];snprintf(lbl,8,"CAR%d",i+1);
            if(ImGui::SmallButton(lbl))c.car_pump[i]=!c.car_pump[i];
            ImGui::PopStyleColor();ImGui::PopID();
        }
        ImGui::Text("  SJAE:"); ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol_Button,c.sjae_running?Col::GREEN:ImVec4{0.2f,0.2f,0.2f,1.f});
        if(ImGui::SmallButton(c.sjae_running?"SJAE:ON":"SJAE:OFF"))c.sjae_running=!c.sjae_running;
        ImGui::PopStyleColor();
        ImGui::SetNextItemWidth(pw*.33f-8.f);
        float pctd = c.condensate_pump_speed * 100.f;
        if (ImGui::SliderFloat("##cp", &pctd, 0.f, 100.f, "Cond Pump: %.0f%%"))
            c.condensate_pump_speed = pctd / 100.f;
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Column 3: hotwell + generator + controls ──────────────────────────────
    if(ImGui::BeginChild("##tc3",{0.f,0.f},false)){
        Hdr("HOTWELL");
        auto& hw=u.hotwell;
        Row("Hotwell Level", hw.level_m,"%.2f","m",
            hw.lo_level_alarm?Col::RED:(hw.hi_level_alarm?Col::AMBER:Col::GREEN));
        Bar(hw.level_m/3.5f,pw*.38f,Col::GREEN,Col::RED,0.9f);
        Lamp(" LO LEVEL ",hw.lo_level_alarm);ImGui::SameLine();
        Lamp(" HI LEVEL ",hw.hi_level_alarm);
        ImGui::Text("  Makeup: %s   Drain: %s",
                    hw.makeup_valve?"OPEN":"closed",hw.drain_valve?"OPEN":"closed");
        bool mk=hw.makeup_valve,dr=hw.drain_valve;
        if(ImGui::Checkbox("Manual Makeup",&mk))hw.makeup_valve=mk;ImGui::SameLine();
        if(ImGui::Checkbox("Manual Drain",&dr))hw.drain_valve=dr;
        Row("CEP Running",   hw.cep_running?1.f:0.f, hw.cep_running?"YES":"NO","",
            hw.cep_running?Col::GREEN:Col::RED);
        Row("FW Suction",    hw.suction_avail_frac*100.f,"%.0f","%",
            hw.suction_avail_frac<0.3f?Col::RED:(hw.suction_avail_frac<0.8f?Col::AMBER:Col::GREEN));

        Hdr("GENERATOR");
        auto& g=u.generator;
        Row("Frequency",     g.frequency_Hz,    "%.3f","Hz",
            std::abs(g.frequency_Hz-50.f)<0.1f?Col::GREEN:Col::AMBER);
        Row("Active Power",  g.power_MW,        "%.1f","MW",Col::YELLOW);
        Row("Reactive Power",g.reactive_MVAR,   "%.1f","MVAR");
        Row("Voltage",       g.terminal_voltage_kV,"%.2f","kV");

        bool synchOk=egrid.canSync(idx);
        ImGui::TextColored(synchOk?Col::GREEN:Col::GREY,"  Synch check: %s",
                           synchOk?"READY":"NOT READY");
        ImGui::TextColored(g.breaker_closed?Col::GREEN:Col::GREY,"  Breaker: %s",
                           g.breaker_closed?"CLOSED":"OPEN");
        Lamp(" OVERSPEED ",u.overspeed_trip);ImGui::SameLine();
        Lamp(" DIFF TRIP ",g.diff_trip);

        float bw=(pw*.38f-4.f)*.5f;
        if(GreenBtn("Close Breaker",{bw,0.f}))ctrl.cmdCloseBreakerRequest();ImGui::SameLine();
        if(RedBtn("Open Breaker",{bw,0.f})) ctrl.cmdOpenBreaker();

        Hdr("UNIT CONTROLS");
        if(GreenBtn("▶ START UNIT",{pw*.38f,26.f})) ctrl.cmdStart();
        if(RedBtn("■ STOP UNIT",{pw*.38f,24.f}))    ctrl.cmdStop();
        float bw2=(pw*.38f-4.f)*.5f;
        if(RedBtn("TRIP",{bw2,0.f})) ctrl.cmdTrip();ImGui::SameLine();
        if(GreenBtn("RESET",{bw2,0.f}))ctrl.cmdReset();

        // Governor demand
        Hdr("GOVERNOR");
        ImGui::SetNextItemWidth(pw*.38f);
        float pcte = u.governor_demand * 100.f;
        if (ImGui::SliderFloat("##gd", &pcte, 0.f, 100.f, "Demand: %.0f%%"))
            u.governor_demand = pcte / 100.f;
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// ELECTRICAL GRID TAB
// ═══════════════════════════════════════════════════════════════════════════════
static void TabGrid(ElectricalGridSystem& eg,TurbineSystem& turbines)
{
    const GridState& g=eg.grid();
    float pw=ImGui::GetContentRegionAvail().x;
    float cw=pw*0.45f;

    if(ImGui::BeginChild("##grid_l",{cw,0.f},false)){
        Hdr("GRID STATUS");
        ImVec4 fc=(g.frequency_Hz<49.f||g.frequency_Hz>51.f)?Col::RED:Col::GREEN;
        Row("Grid Frequency",  g.frequency_Hz,       "%.4f","Hz",fc);
        Row("Freq Deviation",  g.frequency_deviation_Hz,"%.4f","Hz",
            std::abs(g.frequency_deviation_Hz)>0.5f?Col::RED:Col::AMBER);
        Row("Total Generation",g.total_generation_MW,"%.1f","MW",Col::YELLOW);
        Row("Site Load",       g.total_site_load_MW, "%.1f","MW",Col::GREY);
        Row("Export to Grid",  g.export_MW,          "%.1f","MW",Col::GREEN);
        Row("Import from Grid",g.import_MW,          "%.1f","MW",Col::AMBER);
        Lamp(" UNDERFREQ ",g.underfrequency_alarm);ImGui::SameLine();
        Lamp(" OVERFREQ ",g.overfrequency_alarm);  ImGui::SameLine();
        Lamp(" LOOP ",g.loss_of_offsite_power);

        ImGui::Spacing();
        if(ImGui::Button("Trigger LOOP",{cw-16.f,0.f}))eg.triggerLOOP();
        if(ImGui::Button("Restore Offsite",{cw-16.f,0.f}))eg.restoreOffsite();

        Hdr("GENERATOR BUSES");
        for(int i=0;i<4;i++){
            const auto& b=g.bus[i];
            ImGui::TextColored(b.breaker_closed?Col::GREEN:Col::GREY,
                "  Gen %d: %s  %6.3f Hz  %6.1f MW  Synch:%s",
                i+1, b.breaker_closed?"ON GRID":"offline",
                b.frequency_Hz, b.power_MW,
                b.synch_ok?"READY":"    -");
        }
    }
    ImGui::EndChild();ImGui::SameLine();

    if(ImGui::BeginChild("##grid_r",{0.f,0.f},false)){
        Hdr("SITE LOADS");
        for(int i=0;i<eg.numLoads();i++){
            const SiteLoad& l=eg.siteLoad(i);
            ImGui::TextColored(l.energised?Col::GREEN:Col::GREY,
                " %-28s %5.1f MW  %s",l.name,l.load_MW,
                l.essential?"[ESS]":"     ");
            if(!l.essential){
                ImGui::SameLine(pw*.55f-10.f);
                ImGui::PushID(i+1000);
                bool on=l.energised;
                if(ImGui::Checkbox("##load",&on))eg.setLoadEnergised(i,on);
                ImGui::PopID();
            }
        }
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// MOLTEN SALT TAB
// ═══════════════════════════════════════════════════════════════════════════════
static void TabMoltenSalt(MoltenSaltSystem& salt, TurbineSystem& turbines)
{
    MoltenSaltState& s=salt.saltState();
    float pw=ImGui::GetContentRegionAvail().x;
    float cw=pw*0.33f;

    // ── Tanks ─────────────────────────────────────────────────────────────────
    if(ImGui::BeginChild("##salt_tanks",{cw,0.f},false)){
        Hdr("HOT TANK");
        Row("Temperature",s.hot_tank.temp_K,"%.0f","K",
            s.hot_tank.temp_K>950.f?Col::RED:Col::ORANGE);
        Row("Level",      s.hot_tank.level_m,"%.2f","m",
            s.hot_tank.lo_level_alarm?Col::RED:Col::GREEN);
        Bar(s.hot_tank.level_m/16.f,cw-16.f,Col::GREEN,Col::RED,0.9f);
        Lamp(" HI TEMP ",s.hot_tank.hi_temp_alarm);ImGui::SameLine();
        Lamp(" LO TEMP ",s.hot_tank.lo_temp_alarm);ImGui::SameLine();
        Lamp(" LO LEVEL",s.hot_tank.lo_level_alarm);

        Hdr("COLD TANK");
        Row("Temperature",s.cold_tank.temp_K,"%.0f","K",
            s.cold_tank.lo_temp_alarm?Col::RED:Col::CYAN);
        Row("Level",      s.cold_tank.level_m,"%.2f","m",
            s.cold_tank.lo_level_alarm?Col::RED:Col::GREEN);
        Bar(s.cold_tank.level_m/16.f,cw-16.f,Col::CYAN,Col::RED,0.9f);
        Lamp(" LO TEMP ",s.cold_tank.lo_temp_alarm);ImGui::SameLine();
        Lamp(" LO LEVEL",s.cold_tank.lo_level_alarm);

        Hdr("BLANKET CIRC PUMPS");
        ImGui::TextColored(Col::GREY," Blanket outlet: %.0f K  (target ~900 K)",
                           s.blanket_outlet_K);
        for(int i=0;i<2;i++){
            auto& p=s.blanket_circ[i];
            ImGui::PushID(i+2000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," BCP-%d: %s  %.0f kg/s",
                               i+1,p.running?"RUN":"STOP",p.flow_kg_s);
            float bw=(cw-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;p.auto_start=false;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f})){p.running=false;p.auto_start=false;}
            ImGui::Checkbox("AUTO (hold 900 K outlet)",&p.auto_start);
            ImGui::SetNextItemWidth(cw-16.f);
            //  The slider edits the 0-1 fraction but must DISPLAY percent —
            //  the old format string printed the raw fraction, so 100%
            //  speed showed as "Speed: 1%".
            float bcp_pct = p.speed_frac * 100.f;
            if(ImGui::SliderFloat("##bcps",&bcp_pct,0.f,100.f,"Speed: %.0f%%"))
                { p.speed_frac = bcp_pct / 100.f; p.auto_start = false; }
            ImGui::PopID();
        }
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Distribution ──────────────────────────────────────────────────────────
    if(ImGui::BeginChild("##salt_dist",{cw,0.f},false)){
        Hdr("DISTRIBUTION — 1&2 GROUP");
        ImGui::Checkbox("Enable 1&2 Flow",&s.dist_12_enabled);
        ImGui::SetNextItemWidth(cw-16.f);
        float pct1 = s.dist_1_frac * 100.f;
        if (ImGui::SliderFloat("##d1f", &pct1, 0.f, 100.f, "T1 split: %.0f%%"))
            s.dist_1_frac = pct1 / 100.f;
        ImGui::TextColored(Col::GREY, "  T1: %.0f%%  T2: %.0f%%",
                        s.dist_1_frac * 100.f, (1.f - s.dist_1_frac) * 100.f);
        // Show computed flows
        float f12=(s.hotleg[0].flow_kg_s+s.hotleg[1].flow_kg_s)*(s.dist_12_enabled?1.f:0.f);
        ImGui::TextColored(Col::ORANGE,"  T1 flow: %.0f kg/s  →  %.0f MW",
                           f12*s.dist_1_frac, s.sg_heat_MW[0]);
        ImGui::TextColored(Col::ORANGE,"  T2 flow: %.0f kg/s  →  %.0f MW",
                           f12*(1.f-s.dist_1_frac), s.sg_heat_MW[1]);

        Hdr("DISTRIBUTION — 3&4 GROUP");
        ImGui::Checkbox("Enable 3&4 Flow",&s.dist_34_enabled);
        ImGui::SetNextItemWidth(cw-16.f);
        float pct3 = s.dist_3_frac * 100.f;
        if (ImGui::SliderFloat("##d2f", &pct3, 0.f, 100.f, "T3 split: %.0f%%"))
            s.dist_3_frac = pct3 / 100.f;
        ImGui::TextColored(Col::GREY, "  T3: %.0f%%  T4: %.0f%%",
                        s.dist_3_frac * 100.f, (1.f - s.dist_3_frac) * 100.f);
        float f34=(s.hotleg[2].flow_kg_s+s.hotleg[3].flow_kg_s)*(s.dist_34_enabled?1.f:0.f);
        ImGui::TextColored(Col::ORANGE,"  T3 flow: %.0f kg/s  →  %.0f MW",
                           f34*s.dist_3_frac, s.sg_heat_MW[2]);
        ImGui::TextColored(Col::ORANGE,"  T4 flow: %.0f kg/s  →  %.0f MW",
                           f34*(1.f-s.dist_3_frac), s.sg_heat_MW[3]);

        Hdr("SG SALT CONDITIONS");
        for(int i=0;i<4;i++){
            float demand = turbines.unit(i).sgDemandFactor();
            float t_sat  = turbines.unit(i).s.sg_steam_temp_K;
            ImGui::TextColored(Col::ORANGE," SG%d: %.0f→%.0f K  Q=%.0f MW",
                i+1,s.sg_salt_inlet_K[i],s.sg_salt_outlet_K[i],s.sg_heat_MW[i]);
            ImGui::TextColored(Col::GREY,"      steam side T_sat=%.0f K  demand=%.0f%%",
                t_sat, demand*100.f);
        }
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Pumps ─────────────────────────────────────────────────────────────────
    if(ImGui::BeginChild("##salt_pumps",{0.f,0.f},false)){
        Hdr("HOTLEG PUMPS");
        const char* hnames[4]={"HL-1A (Grp 1&2)","HL-1B (Grp 1&2)","HL-2A (Grp 3&4)","HL-2B (Grp 3&4)"};
        for(int i=0;i<4;i++){
            auto& p=s.hotleg[i];
            ImGui::PushID(i+3000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," %s",hnames[i]);
            ImGui::TextColored(Col::GREY,"   %.0f kg/s  %.2f MW",p.flow_kg_s,p.power_MW);
            float bw=(pw*.34f-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;p.auto_start=false;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f})){p.running=false;p.auto_start=false;}
            ImGui::SameLine();ImGui::Checkbox("AUTO",&p.auto_start);
            ImGui::SetNextItemWidth(pw*.34f-16.f);
            float hl_pct = p.speed_frac * 100.f;   // display %, store fraction
            if(ImGui::SliderFloat("##hls",&hl_pct,0.f,100.f,"%.0f%%"))
                { p.speed_frac = hl_pct / 100.f; p.auto_start = false; }
            ImGui::PopID();
        }

        Hdr("COLDLEG PUMPS");
        for(int i=0;i<4;i++){
            auto& p=s.coldleg[i];
            ImGui::PushID(i+4000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," CL-%d (SG%d return)",i+1,i+1);
            ImGui::TextColored(Col::GREY,"   %.0f kg/s",p.flow_kg_s);
            float bw=(pw*.34f-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;p.auto_start=false;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f})){p.running=false;p.auto_start=false;}
            ImGui::SameLine();ImGui::Checkbox("AUTO",&p.auto_start);
            ImGui::SetNextItemWidth(pw*.34f-16.f);
            float cl_pct = p.speed_frac * 100.f;   // display %, store fraction
            if(ImGui::SliderFloat("##cls",&cl_pct,0.f,100.f,"%.0f%%"))
                { p.speed_frac = cl_pct / 100.f; p.auto_start = false; }
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// HELIUM SYSTEM TAB
// ═══════════════════════════════════════════════════════════════════════════════
static void TabHelium(HeliumCoolingSystem& he)
{
    HeliumSystemState& s=he.heState();
    float pw=ImGui::GetContentRegionAvail().x;
    float cw=pw*0.34f;

    // ── Reactor He circuit ────────────────────────────────────────────────────
    if(ImGui::BeginChild("##he_reactor",{cw,0.f},false)){
        Hdr("REACTOR He CIRCUIT (WARM)");
        auto& rc=s.reactor_circuit;
        Row("Inlet Temp",  rc.inlet_temp_K,  "%.0f","K",Col::CYAN);
        Row("Outlet Temp", rc.outlet_temp_K, "%.0f","K",
            rc.hi_outlet_temp_alarm?Col::RED:Col::AMBER);
        Row("Pressure",    rc.pressure_MPa,  "%.1f","MPa");
        Row("Heat Removed",rc.heat_removed_MW,"%.1f","MW",Col::AMBER);
        Lamp(" HI OUTLET ",rc.hi_outlet_temp_alarm);ImGui::SameLine();
        Lamp(" LO FLOW ",  rc.lo_flow_alarm);

        Hdr("REACTOR He PUMPS");
        for(int i=0;i<2;i++){
            auto& p=rc.pump[i];
            ImGui::PushID(i+5000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," He Pump %d: %s  %.0f kg/s  %.2f MW",
                               i+1,p.running?"RUN":"STOP",p.flow_kg_s,p.power_MW);
            float bw=(cw-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f}))p.running=false;
            ImGui::SetNextItemWidth(cw-16.f);
            //  Edit the 0-1 fraction, display percent (was "Speed: 1%" at full)
            float he_pct = p.speed_frac * 100.f;
            if(ImGui::SliderFloat("##rps",&he_pct,0.f,100.f,"Speed: %.0f%%"))
                p.speed_frac = he_pct / 100.f;
            ImGui::PopID();
        }
        Hdr("AUX HEAT EXCHANGERS");
        ImGui::TextColored(Col::GREY,"  Air cooler duty: ");ImGui::SameLine();
        ImGui::TextColored(Col::ORANGE,"%.1f MW",rc.aux_HX_duty_MW);
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Cryostat ──────────────────────────────────────────────────────────────
    if(ImGui::BeginChild("##he_cryo",{cw,0.f},false)){
        Hdr("CRYOSTAT");
        auto& cv=s.cryostat;
        ImVec4 vc=cv.vacuum_ok?Col::GREEN:Col::RED;
        Row("Vacuum Pressure", cv.vacuum_pressure_Pa,"%.2e","Pa",vc);
        Row("Thermal Shield",  cv.thermal_shield_K,  "%.0f","K",Col::CYAN);
        Row("Heat Leak",       cv.heat_leak_W,       "%.0f","W",Col::AMBER);
        Row("Inner Diameter",  cv.inner_diameter_m,  "%.1f","m",Col::GREY);
        Lamp(" VACUUM OK ",cv.vacuum_ok);ImGui::SameLine();
        Lamp(" VACUUM FAIL ",!cv.vacuum_ok);

        Hdr("VACUUM PUMPING");
        ImGui::TextColored(cv.roughing_pump_on?Col::GREEN:Col::GREY,
                           "  Roughing Pump: %s",cv.roughing_pump_on?"ON":"off");
        ImGui::TextColored(cv.turbo_pump_on?Col::GREEN:Col::GREY,
                           "  Turbo-molecular Pump: %s",cv.turbo_pump_on?"ON":"off");
        float bw=(cw-16.f)*.5f;
        if(GreenBtn("Start Roughing",{bw,0.f}))cv.roughing_pump_on=true;ImGui::SameLine();
        if(RedBtn("Stop Roughing",{bw,0.f}))   cv.roughing_pump_on=false;
        if(GreenBtn("Start Turbo",{bw,0.f}))   cv.turbo_pump_on=true;ImGui::SameLine();
        if(RedBtn("Stop Turbo",{bw,0.f}))      cv.turbo_pump_on=false;

        // Visual vacuum scale (log)
        ImGui::Spacing();
        float log_pres=std::clamp(
            (std::log10(cv.vacuum_pressure_Pa)+8.f)/8.f,0.f,1.f);
        ImGui::TextColored(Col::GREY,"  Vacuum quality:");ImGui::SameLine();
        Bar(1.f-log_pres,cw-130.f,Col::GREEN,Col::RED,0.1f);
    }
    ImGui::EndChild();ImGui::SameLine();

    // ── Magnet cryo circuit ───────────────────────────────────────────────────
    if(ImGui::BeginChild("##he_magnet",{0.f,0.f},false)){
        Hdr("MAGNET CRYO CIRCUIT (4.5 K)");
        auto& mc=s.magnet_circuit;
        ImVec4 tc=(mc.cold_mass_temp_K>5.5f)?Col::RED:
                  (mc.cold_mass_temp_K>5.0f)?Col::AMBER:Col::GREEN;
        Row("Cold Mass Temp", mc.cold_mass_temp_K, "%.3f","K",tc);
        Row("He Supply",      mc.supply_temp_K,    "%.3f","K",Col::CYAN);
        Row("He Return",      mc.return_temp_K,    "%.3f","K",Col::CYAN);
        Row("Pressure",       mc.pressure_MPa,     "%.3f","MPa");
        Row("Heat Load",      mc.total_heat_load_W,"%.0f","W",Col::AMBER);
        Row("Cryo Power",     mc.cryo_refrigerator_MW,"%.2f","MW",Col::YELLOW);
        Lamp(" LO TEMP ",mc.lo_temp_alarm);ImGui::SameLine();
        Lamp(" HI TEMP ",mc.hi_temp_alarm);ImGui::SameLine();
        Lamp(" LO PRES ",mc.lo_pressure_alarm);

        Hdr("CRYOGENIC REFRIGERATOR");
        ImGui::TextColored(mc.refrigerator_on?Col::GREEN:Col::GREY,
                           "  Cryo Refrigerator: %s",mc.refrigerator_on?"RUNNING":"stopped");
        ImGui::TextColored(Col::GREY,"  COP at 4.5K: ~1/250");
        float bw=(pw*.32f-16.f)*.5f;
        if(GreenBtn("Start Cryoplant",{bw,0.f}))he.startCryoplant();ImGui::SameLine();
        if(RedBtn("Stop Cryoplant",{bw,0.f}))   he.stopCryoplant();

        Hdr("CRYO CIRCULATORS");
        for(int i=0;i<2;i++){
            auto& p=mc.cryo_pump[i];
            ImGui::PushID(i+6000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY,
                " Circ %d: %s  %.2f kg/s",i+1,p.running?"RUN":"STOP",p.flow_kg_s);
            if(GreenBtn("On",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Off",{bw,0.f}))p.running=false;
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// H&CD TAB — Heating and Current Drive actuators (NBI / ICRH / ECRH / LHCD)
// ═══════════════════════════════════════════════════════════════════════════════
//  Gives the operator independent control over the four H&CD systems.
//  Each system has an enable checkbox, a power setpoint slider, the actual
//  delivered power, and a fault indicator.  The total aux heat (sum of all
//  four actuals) feeds back into the plasma power balance via sp_aux_heat_MW.
static void TabHCD(ReactorState& s, HCDSystem& hcd)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float cw = pw * 0.24f;

    // ── Summary column ─────────────────────────────────────────────────────────
    if(ImGui::BeginChild("##hcd_summary",{cw,0.f},false)){
        Hdr("TOTAL AUX HEAT");
        Row("Setpoint (PID)", s.sp_aux_heat_MW, "%.1f", "MW", Col::AMBER);
        Row("Actual", s.hcd_nbi_actual_MW + s.hcd_icrh_actual_MW
                       + s.hcd_ecrh_actual_MW + s.hcd_lhcd_actual_MW,
                       "%.1f", "MW", Col::GREEN);
        Row("Q_scientific", s.Q_scientific, "%.2f", "");
        ImGui::Spacing();

        // ── Heating presets ──────────────────────────────────────────────────
        //  One-click buttons to set up common heating configurations.  These
        //  are a convenience so the operator doesn't have to individually
        //  enable and set each system's power.  The "Standard Heating" preset
        //  gives the typical ITER ramp-up mix (~32 MW total) which is enough
        //  to reach Burning status from a cold start.
        Hdr("PRESETS");
        if(GreenBtn("Standard Heating",{cw-16.f,0.f})){
            // ITER-class ramp-up: NBI full + ICRH half + ECRH quarter
            s.hcd_nbi_on  = true;  s.hcd_nbi_setpoint_MW  = 16.5f;
            s.hcd_icrh_on = true;  s.hcd_icrh_setpoint_MW = 10.0f;
            s.hcd_ecrh_on = true;  s.hcd_ecrh_setpoint_MW = 6.0f;
            s.hcd_lhcd_on = false; s.hcd_lhcd_setpoint_MW = 0.0f;
        }
        if(GreenBtn("ECRH Startup Only",{cw-16.f,0.f})){
            // Minimal: just ECRH for breakdown + early ramp
            s.hcd_nbi_on  = false; s.hcd_nbi_setpoint_MW  = 0.0f;
            s.hcd_icrh_on = false; s.hcd_icrh_setpoint_MW = 0.0f;
            s.hcd_ecrh_on = true;  s.hcd_ecrh_setpoint_MW = 5.0f;
            s.hcd_lhcd_on = false; s.hcd_lhcd_setpoint_MW = 0.0f;
        }
        if(GreenBtn("Full Power",{cw-16.f,0.f})){
            // Everything at max — for Q>10 scenarios
            s.hcd_nbi_on  = true;  s.hcd_nbi_setpoint_MW  = 16.5f;
            s.hcd_icrh_on = true;  s.hcd_icrh_setpoint_MW = 20.0f;
            s.hcd_ecrh_on = true;  s.hcd_ecrh_setpoint_MW = 24.0f;
            s.hcd_lhcd_on = true;  s.hcd_lhcd_setpoint_MW = 20.0f;
        }
        if(RedBtn("All Off",{cw-16.f,0.f})){
            s.hcd_nbi_on  = false; s.hcd_nbi_setpoint_MW  = 0.0f;
            s.hcd_icrh_on = false; s.hcd_icrh_setpoint_MW = 0.0f;
            s.hcd_ecrh_on = false; s.hcd_ecrh_setpoint_MW = 0.0f;
            s.hcd_lhcd_on = false; s.hcd_lhcd_setpoint_MW = 0.0f;
        }
        ImGui::Spacing();

        Hdr("CURRENT DRIVE");
        Row("NBI CD",  s.hcd_nbi_current_drive_MA,  "%.2f", "MA", Col::CYAN);
        Row("LHCD CD", s.hcd_lhcd_current_drive_MA, "%.2f", "MA", Col::CYAN);
        Row("Bootstrap", s.hcd_bootstrap_current_MA, "%.2f", "MA", Col::CYAN);
        float total_cd = s.hcd_nbi_current_drive_MA + s.hcd_lhcd_current_drive_MA
                       + s.hcd_bootstrap_current_MA;
        Row("Total non-inductive", total_cd, "%.2f", "MA",
            total_cd / std::max(0.01f, s.plasma_current_MA) > 0.5f ? Col::GREEN : Col::AMBER);
        ImGui::TextColored(Col::GREY,"  (Fraction of I_p: %.0f%%)",
            100.f * total_cd / std::max(0.01f, s.plasma_current_MA));
        // ── CS flux budget: why non-inductive current matters ────────────────
        ImVec4 fxc = s.cs_flux_exhausted ? Col::RED
                   : (s.cs_flux_remaining_Wb < 20.f ? Col::AMBER : Col::GREEN);
        Row("CS Flux Left", s.cs_flux_remaining_Wb, "%.1f", "Wb", fxc);
        if(s.cs_flux_exhausted)
            ImGui::TextColored(Col::RED, "  FLUX EXHAUSTED — I_p held by CD only");

        if(s.nbi_shinethrough_block && s.hcd_nbi_on){
            ImGui::Spacing();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(Col::AMBER,
                "  ⚠ NBI SHINE-THROUGH INTERLOCK: n_e < 1.5e19 m-3 — beams "
                "blocked until the plasma is dense enough to stop them.");
            ImGui::PopTextWrapPos();
        }

        Hdr("FAULTS");
        Lamp(" NBI ",  hcd.nbiFault());
        ImGui::SameLine(); Lamp(" ICRH ", hcd.icrhFault());
        ImGui::SameLine(); Lamp(" ECRH ", hcd.ecrhFault());
        ImGui::SameLine(); Lamp(" LHCD ", hcd.lhcdFault());
        if(s.alarm_aux_heat_fault){
            ImGui::Spacing();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(Col::RED, "  Active fault — see below for reason");
            ImGui::PopTextWrapPos();
            if(hcd.nbiFault())  ImGui::TextColored(Col::AMBER,"  NBI:  %s", hcd.nbiFaultReason().c_str());
            if(hcd.icrhFault()) ImGui::TextColored(Col::AMBER,"  ICRH: %s", hcd.icrhFaultReason().c_str());
            if(hcd.ecrhFault()) ImGui::TextColored(Col::AMBER,"  ECRH: %s", hcd.ecrhFaultReason().c_str());
            if(hcd.lhcdFault()) ImGui::TextColored(Col::AMBER,"  LHCD: %s", hcd.lhcdFaultReason().c_str());
        }
        ImGui::Spacing();
        if(GreenBtn("Clear Faults",{cw-16.f,0.f})) hcd.clearFaults();

        // Inject fault buttons (for training/testing)
        Hdr("INJECT FAULT (TEST)");
        float bw = (cw-20.f)*0.5f;
        if(ImGui::SmallButton("NBI##f"))  hcd.faultNBI("Beamline valve closed (operator test)");
        ImGui::SameLine();
        if(ImGui::SmallButton("ICRH##f")) hcd.faultICRH("Antenna VSWR > 3.0 (operator test)");
        if(ImGui::SmallButton("ECRH##f")) hcd.faultECRH("Gyrotron collector overtemp (operator test)");
        ImGui::SameLine();
        if(ImGui::SmallButton("LHCD##f")) hcd.faultLHCD("Klystron HV trip (operator test)");
    }
    ImGui::EndChild(); ImGui::SameLine();

    // ── Per-system controls ────────────────────────────────────────────────────
    auto renderSystem = [&](const char* name, bool& on, float& setpoint, float actual,
                            float max_MW, const char* description, bool faulted,
                            bool ready, float warmup_remaining) {
        if(ImGui::BeginChild(name,{cw,0.f},false)){
            Hdr(name);
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(Col::GREY, "%s", description);
            ImGui::PopTextWrapPos();
            ImGui::Spacing();

            ImGui::Checkbox("Enabled", &on);
            ImGui::SetNextItemWidth(cw-16.f);
            float pct = setpoint / max_MW * 100.f;
            if(ImGui::SliderFloat("##sp", &pct, 0.f, 100.f, "Setpoint: %.0f%%"))
                setpoint = pct / 100.f * max_MW;

            Row("Setpoint", setpoint, "%.1f", "MW", Col::AMBER);
            Row("Actual",   actual,   "%.1f", "MW",
                faulted ? Col::RED : (actual > 0.1f ? Col::GREEN : Col::GREY));
            Row("Max",      max_MW,   "%.1f", "MW", Col::GREY);

            if(warmup_remaining > 0.f){
                ImGui::TextColored(Col::AMBER, "  WARMUP: %.1fs remaining", warmup_remaining);
            } else if(ready){
                ImGui::TextColored(Col::GREEN, "  READY");
            } else {
                ImGui::TextColored(Col::RED, "  NOT READY");
            }

            if(faulted){
                ImGui::Spacing();
                ImGui::TextColored(Col::RED, "  ⚠ FAULT");
            }
        }
        ImGui::EndChild(); ImGui::SameLine();
    };

    renderSystem("NBI",  s.hcd_nbi_on,  s.hcd_nbi_setpoint_MW,
                 s.hcd_nbi_actual_MW, 16.5f,
                 "Neutral Beam Injection. 1 MeV deuterium neutrals, 16.5 MW "
                 "(3 injectors × 5.5 MW). Heats ions; co-current drive ~1 MA. "
                 "Interlocked on n_e > 1.5e19 (shine-through protection).",
                 hcd.nbiFault(),
                 hcd.nbiReady() && !s.nbi_shinethrough_block,
                 hcd.nbiWarmupRemaining());

    renderSystem("ICRH", s.hcd_icrh_on, s.hcd_icrh_setpoint_MW,
                 s.hcd_icrh_actual_MW, 20.0f,
                 "Ion Cyclotron Resonance Heating. 40–55 MHz, 20 MW. "
                 "Minority heating (3He in D-T). No net current drive.",
                 hcd.icrhFault(), !hcd.icrhFault(), 0.f);

    renderSystem("ECRH", s.hcd_ecrh_on, s.hcd_ecrh_setpoint_MW,
                 s.hcd_ecrh_actual_MW, 24.0f,
                 "Electron Cyclotron Resonance Heating. 170 GHz, 24 MW. "
                 "Local electron heating; used for startup assist + NTM control.",
                 hcd.ecrhFault(), !hcd.ecrhFault(), 0.f);

    renderSystem("LHCD", s.hcd_lhcd_on, s.hcd_lhcd_setpoint_MW,
                 s.hcd_lhcd_actual_MW, 20.0f,
                 "Lower Hybrid Current Drive. 5 GHz, 20 MW. Off-axis "
                 "non-inductive current drive for steady-state scenarios.",
                 hcd.lhcdFault(), hcd.lhcdReady(), hcd.lhcdWarmupRemaining());
}

// ═══════════════════════════════════════════════════════════════════════════════
// VACUUM & DM TAB — Vessel pumping + Disruption Mitigation
// ═══════════════════════════════════════════════════════════════════════════════
static void TabVacuumDM(ReactorState& s, VacuumVesselSystem& vac,
                        DisruptionMitigationSystem& dm)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float cw = pw * 0.5f;

    // ── Left column: vacuum vessel ─────────────────────────────────────────────
    if(ImGui::BeginChild("##vac",{cw,0.f},false)){
        Hdr("VESSEL PRESSURE");
        ImVec4 pc = vac.vacuumOK() ? Col::GREEN
                   : (s.vessel_pressure_Pa > 1.f ? Col::RED : Col::AMBER);
        Row("Pressure", s.vessel_pressure_Pa, "%.3e", "Pa", pc);
        Row("Wall Temp", s.vessel_wall_temp_K, "%.0f", "K", Col::AMBER);
        Row("Vacuum OK", vac.vacuumOK()?1.f:0.f,
            vac.vacuumOK()?"YES":"NO", "", vac.vacuumOK()?Col::GREEN:Col::RED);

        // Logarithmic pressure gauge
        ImGui::Spacing();
        float log_pres = std::clamp(
            (std::log10(std::max(s.vessel_pressure_Pa, 1e-8f)) + 8.f) / 13.f, 0.f, 1.f);
        // 1e-8 Pa → 0, 1e5 Pa → 1
        ImGui::TextColored(Col::GREY, "  Vacuum quality (log scale):");
        Bar(1.f - log_pres, cw-130.f, Col::GREEN, Col::RED, 0.05f);
        ImGui::TextColored(Col::GREY, "  Initiation threshold: 1e-3 Pa");

        Hdr("PUMPS");
        ImGui::Checkbox("Roughing Pump", &s.vessel_roughing_on);
        ImGui::Checkbox("Turbo Pump",    &s.vessel_turbo_on);
        // Auto-stop turbo if pressure too high (handled by VacuumVessel update)
        if(s.vessel_turbo_on && s.vessel_pressure_Pa > 10.f){
            ImGui::TextColored(Col::RED, "  ⚠ TURBO INTERLOCK: pressure too high");
        }

        float bw = (cw-20.f)*0.5f;
        if(GreenBtn("Start Roughing",{bw,0.f})) s.vessel_roughing_on = true;
        ImGui::SameLine();
        if(RedBtn("Stop Roughing",{bw,0.f})) s.vessel_roughing_on = false;
        if(GreenBtn("Start Turbo",{bw,0.f})){
            if(s.vessel_pressure_Pa < 10.f) s.vessel_turbo_on = true;
        }
        ImGui::SameLine();
        if(RedBtn("Stop Turbo",{bw,0.f})) s.vessel_turbo_on = false;

        Hdr("BAKEOUT & CONDITIONING");
        ImGui::Checkbox("Wall Bakeout (423 K, 5 min)", &s.vessel_bakeout_on);
        Row("Bakeout Progress", vac.bakeoutProgress()*100.f, "%.1f", "%",
            vac.bakeoutComplete() ? Col::GREEN : Col::AMBER);
        Row("Boronization Thickness", s.boronization_thickness_nm, "%.0f", "nm", Col::CYAN);
        ImGui::Spacing();
        if(GreenBtn("Start Boronization",{cw-16.f,0.f})){
            // Boronization requires good vacuum and cold plasma
            if(s.vessel_vacuum_ok && s.plasma_status == PlasmaStatus::Cold){
                s.boronization_in_progress = true;
            }
        }
        if(s.boronization_in_progress){
            ImGui::TextColored(Col::AMBER, "  Boronization in progress — vacuum poor");
        }

        Hdr("VESSEL BREACH (LOVA)");
        if(s.alarm_vessel_breach){
            ImGui::TextColored(Col::RED, "  ⚠ VESSEL BREACH DETECTED");
            if(GreenBtn("Clear Breach (test)",{cw-16.f,0.f})) vac.clearBreach();
        } else {
            if(RedBtn("Inject Breach (test)",{cw-16.f,0.f})) vac.triggerVesselBreach();
        }
    }
    ImGui::EndChild(); ImGui::SameLine();

    // ── Right column: disruption mitigation ────────────────────────────────────
    if(ImGui::BeginChild("##dm",{0.f,0.f},false)){
        Hdr("DISRUPTION STATUS");
        Row("Disruption Flag", s.disruption_flag?1.f:0.f,
            s.disruption_flag?"ACTIVE":"clear", "",
            s.disruption_flag?Col::RED:Col::GREEN);
        Row("Cause", 0.f, disruptionCauseStr(s.disruption_cause), "",
            s.disruption_flag?Col::RED:Col::GREY);
        Row("q_95", s.q_safety, "%.2f", "",
            s.q_safety < 2.f ? Col::RED : (s.q_safety < 2.5f ? Col::AMBER : Col::GREEN));
        Row("β_N", s.beta_N, "%.2f", "",
            s.beta_N > 2.5f ? Col::RED : (s.beta_N > 2.0f ? Col::AMBER : Col::GREEN));
        Row("Plasma Status", 0.f,
            (s.plasma_status==PlasmaStatus::Burning) ? "BURNING" :
            (s.plasma_status==PlasmaStatus::Disrupting) ? "DISRUPTING" :
            (s.plasma_status==PlasmaStatus::Quenched) ? "QUENCHED" :
            (s.plasma_status==PlasmaStatus::Initiating) ? "INITIATING" : "COLD",
            "", s.disruption_flag ? Col::RED : Col::GREEN);

        Hdr("MGI — Massive Gas Injection");
        ImGui::Checkbox("Arm MGI", &s.mgi_armed);
        Lamp(" ARMED ", s.mgi_armed);
        ImGui::SameLine(); Lamp(" FIRED ", s.mgi_fired);
        ImGui::TextColored(Col::GREY,
            "  Injects Ne/Ar gas for radiative thermal quench.");
        ImGui::TextColored(Col::GREY,
            "  Reduces halo current forces by ~66%%.");
        ImGui::TextColored(Col::GREY,
            "  Vessel pressure rises to ~10 Pa after firing.");
        if(s.mgi_armed && !s.mgi_fired){
            if(RedBtn("⚡ FIRE MGI",{cw*0.95f-16.f,28.f})){
                s.cmd_disrupt_mitigation = true;
            }
        }

        Hdr("SPI — Shattered Pellet Injection");
        ImGui::Checkbox("Arm SPI", &s.spi_armed);
        Lamp(" ARMED ", s.spi_armed);
        ImGui::SameLine(); Lamp(" FIRED ", s.spi_fired);
        ImGui::TextColored(Col::GREY,
            "  Fires 500 mg D2 pellet at 250 m/s into shatter ring.");
        ImGui::TextColored(Col::GREY,
            "  Reduces halo current forces by ~75%% (better than MGI).");
        ImGui::TextColored(Col::GREY,
            "  Preferred for high-current (ITER-class) plasmas.");
        if(s.spi_armed && !s.spi_fired){
            if(RedBtn("⚡ FIRE SPI",{cw*0.95f-16.f,28.f})){
                s.cmd_disrupt_mitigation = true;
            }
        }

        Hdr("AUTO-MGI");
        //  The checkbox binds to state.dm_auto_mgi (persistent).  It used
        //  to bind to a local `bool auto_mgi = false` recreated every
        //  frame, so it rendered unticked and immediately disabled
        //  auto-MGI again — enabling it was impossible.
        ImGui::Checkbox("Enable auto-MGI (fires on disruption or q_95 < 1.8)",
                        &s.dm_auto_mgi);
        if(s.dm_auto_mgi && !s.mgi_armed)
            ImGui::TextColored(Col::AMBER, "  ⚠ Auto-MGI needs MGI ARMED to fire");
        if(s.dm_injected_impurity > 0.001f)
            ImGui::TextColored(Col::RED, "  RADIATOR IN VESSEL: %.1f%% (pumping out)",
                               s.dm_injected_impurity*100.f);
        if(s.dm_mitigated)
            ImGui::TextColored(Col::CYAN, "  Current quench under mitigation control");

        Hdr("MITIGATION STATUS");
        Row("Active", dm.active()?1.f:0.f, dm.active()?"YES":"no", "",
            dm.active()?Col::RED:Col::GREEN);
        Row("Last Fire Type", 0.f, dm.lastFireType().c_str(), "", Col::GREY);
        Row("Time Since Fire", dm.timeSinceLastFire_s(), "%.1f", "s", Col::GREY);
        Row("Halo Force Reduction", s.mitigation_force_N/1e6f, "%.1f", "MN", Col::CYAN);
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// TRITIUM PLANT TAB — TES / Detritiation / Accountancy
// ═══════════════════════════════════════════════════════════════════════════════
static void TabTritiumPlant(ReactorState& s, TritiumPlantSystem& tp)
{
    float pw = ImGui::GetContentRegionAvail().x;
    float cw = pw * 0.5f;

    if(ImGui::BeginChild("##tp_l",{cw,0.f},false)){
        Hdr("TRITIUM EXTRACTION SYSTEM (TES)");
        Row("T in Plant",  s.tritium_in_plant_g,   "%.1f", "g", Col::AMBER);
        Row("T in Fuel Store", s.fuel_T_inventory_g, "%.1f", "g", Col::GREEN);
        Row("T Recovery Rate", s.tritium_recovery_rate_g_s, "%.4f", "g/s", Col::CYAN);
        Row("TBR (current)", s.tbr_current, "%.3f", "",
            s.tbr_current >= 1.0f ? Col::GREEN : Col::AMBER);
        ImGui::Spacing();
        if(GreenBtn("Start TES",{cw*0.48f,0.f})) tp.startTES();
        ImGui::SameLine();
        if(RedBtn("Stop TES",{cw*0.48f,0.f}))   tp.stopTES();
        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,
            "  TES extracts T bred in the Li blanket via He purge gas.");
        ImGui::TextColored(Col::GREY,
            "  Throughput scales with blanket coolant flow.");
        ImGui::TextColored(Col::GREY,
            "  TBR > 1.05 needed for self-sufficient fuel cycle.");
    }
    ImGui::EndChild(); ImGui::SameLine();

    if(ImGui::BeginChild("##tp_r",{0.f,0.f},false)){
        Hdr("DETRITIATION");
        Row("Detritiation Flow", s.detritiation_flow_m3_s, "%.1f", "m³/s", Col::CYAN);
        Row("Detritiation On", tp.detritiationRunning()?1.f:0.f,
            tp.detritiationRunning()?"YES":"NO", "",
            tp.detritiationRunning()?Col::GREEN:Col::GREY);
        ImGui::Spacing();
        if(GreenBtn("Start Detritiation",{cw*0.48f,0.f})) tp.startDetritiation();
        ImGui::SameLine();
        if(RedBtn("Stop Detritiation",{cw*0.48f,0.f})) tp.stopDetritiation();
        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,
            "  Detritiation cleans tritiated air from the reactor building.");
        ImGui::TextColored(Col::GREY,
            "  Activated after a tritium release or during maintenance.");

        Hdr("TRITIUM ACCOUNTANCY");
        Row("Plant Inventory",  s.tritium_in_plant_g,    "%.1f", "g", Col::AMBER);
        Row("Fuel Store",       s.fuel_T_inventory_g,    "%.1f", "g", Col::GREEN);
        Row("In Vessel",        s.vessel_T_g,            "%.4f", "g", Col::CYAN);
        float total = s.tritium_in_plant_g + s.fuel_T_inventory_g + s.vessel_T_g;
        Row("TOTAL (inventory)", total, "%.1f", "g",
            total > 150.f ? Col::RED : Col::GREEN);
        ImGui::Spacing();
        ImGui::TextColored(Col::GREY,
            "  Regulatory limit: 200 g on-site (ITER-class).");
        ImGui::TextColored(Col::GREY,
            "  Accountancy tracks all T: plant + fuel store + in-vessel.");
    }
    ImGui::EndChild();
}

// ═══════════════════════════════════════════════════════════════════════════════
// OPERATIONS TAB — Step-by-step startup / shutdown procedures with checklists
// ═══════════════════════════════════════════════════════════════════════════════
//  This is the "feel of operating the whole reactor" piece — the operator
//  gets a real procedure to follow, with each step's prerequisites and
//  completion criteria shown.  Steps auto-check when their criteria are met.
static void TabOperations(const ReactorState& s, const SimTime& t)
{
    float pw = ImGui::GetContentRegionAvail().x;

    Hdr("STARTUP PROCEDURE");
    ImGui::TextColored(Col::GREY,
        "  Follow these steps in order to initiate a plasma discharge.");
    ImGui::Spacing();

    auto stepLine = [](int n, const char* title, const char* detail,
                       bool done, bool blocked) {
        ImVec4 iconCol = done ? Col::GREEN : (blocked ? Col::RED : Col::AMBER);
        const char* icon = done ? "[OK]" : (blocked ? "[!!]" : "[  ]");
        ImGui::TextColored(iconCol, "%s", icon); ImGui::SameLine();
        ImGui::TextColored(done ? Col::GREEN : Col::WHITE, "Step %d: %s", n, title);
        ImGui::TextColored(Col::GREY, "      %s", detail);
        ImGui::Spacing();
    };

    // Step 1: Cryoplant running (magnet He temp < 5 K)
    bool s1 = s.cryo_ok && s.magnet_temp_K < 5.5f;
    stepLine(1, "Start cryoplant",
             "Cryoplant must be running and magnet cold mass < 5.5 K.  "
             "Go to HELIUM tab → Start Cryoplant.  Wait for cold mass to reach 4.5 K.",
             s1, !s1 && s.magnet_temp_K > 8.f);

    // Step 2: Magnets ramped (B_T > 95% setpoint)
    bool s2 = s.B_toroidal_T > 0.95f * s.sp_B_toroidal_T && s.sp_B_toroidal_T > 0.5f;
    stepLine(2, "Ramp toroidal field",
             "Set B_T setpoint (left panel slider) to 5.3 T.  Wait for coil "
             "current to ramp (limited to 2 kA/s).  Field must be > 95%% of setpoint.",
             s2, s1 && !s2);

    // Step 3: Vacuum vessel at base pressure
    bool s3 = s.vessel_vacuum_ok;
    stepLine(3, "Pump down vessel",
             "Go to VACUUM&DM tab → Start Roughing Pump.  Wait for pressure < 1 Pa, "
             "then Start Turbo Pump.  Wait for pressure < 1e-3 Pa.  "
             "If pressure plateaus, run Wall Bakeout (5 min at 423 K).",
             s3, s2 && !s3);

    // Step 4: Fuel inventory available
    bool s4 = (s.fuel_D_inventory_g > 1.0f) && (s.fuel_T_inventory_g > 1.0f);
    stepLine(4, "Verify fuel inventory",
             "Deuterium and tritium inventories must be > 1 g.  Use +D2 and +T "
             "buttons on left panel if low.  Consider starting TES on TRITIUM tab "
             "to extract T from blanket (only meaningful after first discharge).",
             s4, s3 && !s4);

    // Step 5: Auxiliary heating armed
    bool s5 = s.hcd_ecrh_on || s.hcd_nbi_on;
    stepLine(5, "Arm ECRH for breakdown assist",
             "Go to H&CD tab → click 'ECRH Startup Only' preset (or manually "
             "enable ECRH at ~5 MW).  ECRH pre-ionizes the gas, lowering the "
             "loop voltage needed for breakdown.  NOTE: ECRH is auto-enabled "
             "at 5 MW when you press INITIATE PLASMA, so this step is optional "
             "— but it's good practice to verify it's armed before initiation.",
             s5, s4 && !s5);

    // Step 6: INITIATE PLASMA button
    bool s6 = s.plasma_status != PlasmaStatus::Cold;
    stepLine(6, "Initiate plasma",
             "All preconditions met — press INITIATE PLASMA on the left panel.  "
             "The central solenoid induces a loop voltage, breakdown occurs, "
             "and plasma current begins ramping toward setpoint (15 MA).  "
             "ECRH is auto-enabled at 5 MW for breakdown assist.  Ohmic heating "
             "provides 15 MW during Initiating to ramp T_e up.",
             s6, s5 && !s6);

    // Step 7: Ramp to burning
    bool s7 = s.plasma_status == PlasmaStatus::Burning;
    stepLine(7, "Ramp to burning",
             "Once I_p > 13 MA and T_e > 10 keV, go to H&CD tab → click "
             "'Standard Heating' preset (NBI 16.5 + ICRH 10 + ECRH 6 = 32 MW).  "
             "This drives T_e to 20 keV.  Once fusion power > 10 MW, alpha "
             "heating self-sustains and status flips to BURNING.  "
             "Then increase ne toward 1e20 m⁻³ for full power.",
             s7, s6 && !s7);

    Hdr("STEADY-STATE OPERATION");
    ImGui::TextColored(Col::GREY,
        "  Monitor: q_95 > 2.0, β_N < 2.5, n_e < Greenwald limit,");
    ImGui::TextColored(Col::GREY,
        "  divertor temp < 2500 K, first wall < 2000 K, magnet temp < 5.5 K.");
    ImGui::TextColored(Col::GREY,
        "  Adjust H&CD mix to control T_e.  Pellet frequency for density.");
    ImGui::TextColored(Col::GREY,
        "  Start TES once TBR > 1 to keep tritium inventory topped up.");
    ImGui::Spacing();

    Hdr("CONTROLLED SHUTDOWN PROCEDURE");
    ImGui::TextColored(Col::GREY,
        "  Option A (preferred): CONTROLLED SHUTDOWN button (left panel).");
    ImGui::TextColored(Col::GREY,
        "    - Ramps down I_p, T_e, ne, P_aux simultaneously at safe rates.");
    ImGui::TextColored(Col::GREY,
        "    - Plasma decays smoothly to Quenched.  No thermal quench, no halo currents.");
    ImGui::Spacing();
    ImGui::TextColored(Col::GREY,
        "  Option B (emergency): SCRAM button or F1 key.");
    ImGui::TextColored(Col::GREY,
        "    - Immediately zeros all setpoints.  Magnet quench dump fires.");
    ImGui::TextColored(Col::GREY,
        "    - Use only when the plant is in an unsafe state (overtemp, quench, etc.).");
    ImGui::Spacing();

    Hdr("DISRUPTION MITIGATION PROCEDURE");
    ImGui::TextColored(Col::GREY,
        "  If a disruption is detected (alarm fires or q_95 < 1.5):");
    ImGui::TextColored(Col::GREY,
        "  1. Arm MGI or SPI on VACUUM&DM tab (if not already armed).");
    ImGui::TextColored(Col::GREY,
        "  2. Press FIRE to inject the mitigation pellet/gas.");
    ImGui::TextColored(Col::GREY,
        "  3. Verify halo current forces are reduced (mitigation_force_N > 0).");
    ImGui::TextColored(Col::GREY,
        "  4. Wait for plasma to Quench, then RESET for cold restart.");
    ImGui::Spacing();

    Hdr("CURRENT PROCEDURE STATUS");
    const char* proc = "UNKNOWN";
    ImVec4 pc = Col::GREY;
    if (s.plasma_status == PlasmaStatus::Cold) { proc = "PRE-STARTUP"; pc = Col::GREY; }
    else if (s.plasma_status == PlasmaStatus::Initiating) { proc = "STARTUP"; pc = Col::AMBER; }
    else if (s.plasma_status == PlasmaStatus::Burning) { proc = "OPERATING"; pc = Col::GREEN; }
    else if (s.plasma_status == PlasmaStatus::Disrupting) { proc = "DISRUPTION"; pc = Col::RED; }
    else if (s.plasma_status == PlasmaStatus::Quenched) { proc = "SHUTDOWN"; pc = Col::GREY; }

    if (s.mode == ReactorMode::Rampdown) { proc = "CONTROLLED SHUTDOWN"; pc = Col::AMBER; }
    if (s.mode == ReactorMode::Emergency) { proc = "EMERGENCY"; pc = Col::RED; }

    ImGui::TextColored(pc, "  Current procedure: %s", proc);
    ImGui::TextColored(Col::GREY, "  Sim time: %.1f s  (tick %d)", t.total_s, t.tick);
}

// ═══════════════════════════════════════════════════════════════════════════════
// ALARMS TAB — Full alarm list with cause + remediation
// ═══════════════════════════════════════════════════════════════════════════════
static void TabAlarms(AlarmSystem& alm)
{
    float pw = ImGui::GetContentRegionAvail().x;

    Hdr("ALARM SUMMARY");
    int n_active    = alm.countActive();
    int n_critical  = alm.countCritical();
    int n_total     = (int)alm.log.size();
    ImGui::TextColored(n_active ? Col::RED : Col::GREEN, "  Active: %d", n_active);
    ImGui::SameLine(160.f);
    ImGui::TextColored(n_critical ? Col::RED : Col::GREEN, "  Critical: %d", n_critical);
    ImGui::SameLine(320.f);
    ImGui::TextColored(Col::GREY, "  Total logged: %d", n_total);
    ImGui::Spacing();
    float bw = (pw-20.f)*0.33f;
    if(GreenBtn("ACK ALL",{bw,0.f})) alm.ackAll();
    ImGui::SameLine();
    if(ImGui::Button("CLEAR ACKED",{bw,0.f})) alm.clearOld();
    ImGui::SameLine();
    ImGui::TextColored(Col::GREY, "  (ACK = acknowledge; CLR = remove acked+inactive)");

    Hdr("ACTIVE ALARMS");
    if(n_active == 0){
        ImGui::TextColored(Col::GREEN, "  No active alarms — plant is nominal.");
    } else {
        // Sort: critical first, then by trip time
        std::vector<const AlarmEntry*> active;
        for(const auto& a : alm.log) if(a.active) active.push_back(&a);
        std::sort(active.begin(), active.end(),
            [](const AlarmEntry* a, const AlarmEntry* b){
                if(a->severity != b->severity)
                    return (int)a->severity > (int)b->severity;
                return a->t < b->t;
            });

        for(const auto* a : active){
            ImGui::Spacing();
            // Severity colour: critical=red, warning=amber, info=cyan
            ImVec4 sevCol = (a->severity == AlarmSeverity::Critical) ? Col::RED
                          : (a->severity == AlarmSeverity::Warning)  ? Col::AMBER
                          : Col::CYAN;
            // Header line with severity icon + msg + time
            const char* sevIcon = (a->severity == AlarmSeverity::Critical) ? "⛔"
                                : (a->severity == AlarmSeverity::Warning)  ? "⚠"
                                : "ℹ";
            float bl = (!a->acked && a->severity == AlarmSeverity::Critical
                        && fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f) ? 1.f : 0.f;
            ImGui::TextColored({sevCol.x, sevCol.y*bl + sevCol.y*(1-bl),
                                sevCol.z*bl + sevCol.z*(1-bl), 1.f},
                "  %s  [t=%.1fs]  %s  %s",
                sevIcon, a->t, a->msg.c_str(),
                a->acked ? "(ACK)" : "(UNACK)");
            // Cause
            ImGui::TextColored(Col::WHITE, "    CAUSE: %s", a->cause.c_str());
            // Remediation
            ImGui::TextColored(Col::CYAN, "    ACTION: %s", a->remediation.c_str());
        }
    }

    Hdr("ALARM HISTORY (acked / cleared)");
    bool any_history = false;
    for(const auto& a : alm.log){
        if(a.active) continue;  // only show inactive (history)
        any_history = true;
        ImGui::TextColored(Col::GREY, "  [t=%.1fs]  %s  (cleared)",
            a.t, a.msg.c_str());
    }
    if(!any_history){
        ImGui::TextColored(Col::GREY, "  No alarm history.");
    }
}


// ═══════════════════════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════════════════════
int main(int,char**)
{
    // ── SDL2 + OpenGL setup ───────────────────────────────────────────────────
    if(SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER)){
        fprintf(stderr,"SDL_Init: %s\n",SDL_GetError());return 1;
    }
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS,0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK,SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION,3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION,3);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,24);

    SDL_Window* win=SDL_CreateWindow("FusionSim — Tokamak Power Plant",
        SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,WINDOW_W,WINDOW_H,
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
    if(!win){fprintf(stderr,"Window: %s\n",SDL_GetError());return 1;}

    SDL_GLContext gl=SDL_GL_CreateContext(win);
    SDL_GL_MakeCurrent(win,gl);SDL_GL_SetSwapInterval(1);

    IMGUI_CHECKVERSION();ImGui::CreateContext();
    ImGuiIO& io=ImGui::GetIO();
    io.ConfigFlags|=ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename="fusionsim.ini";
    if(!io.Fonts->AddFontFromFileTTF("assets/fonts/JetBrainsMono-Regular.ttf",13.f))
        io.Fonts->AddFontDefault();
    ApplyTheme();
    ImGui_ImplSDL2_InitForOpenGL(win,gl);
    ImGui_ImplOpenGL3_Init("#version 330");

    // ── Physics modules ───────────────────────────────────────────────────────
    ReactorState state; SimTime sim; sim.dt_s=SIM_DT_S;
    state.sp_plasma_current_MA=15.f; state.sp_electron_temp_keV=20.f;
    state.sp_density_m3=1e20f; state.sp_B_toroidal_T=5.3f;
    state.sp_fuel_rate=0.5f; state.sp_coolant_flow=0.8f; state.D_T_ratio=1.f;
    state.mode=ReactorMode::Startup; state.plasma_status=PlasmaStatus::Cold;
    // Initialize plasma-shape setpoints to ITER nominal
    state.sp_kappa = 1.7f; state.sp_delta = 0.33f;

    ControlConfig cc; ControlSystem  control(cc);
    MagnetConfig  mc; MagnetSystem   magnets(mc);
    FuelConfig    fc; FuelSystem     fuel(fc);
    PlasmaConfig  pc; pc.Nx=64;pc.Ny=64;pc.Nz=64;pc.pic_dt=1e-12f;
    // ── PIC mode ─────────────────────────────────────────────────────────────
    //  run_pic = false → "arcade mode": the sim uses the 0D power-balance
    //  model only (IPB98(y,2) confinement + brem + sync + line radiation),
    //  which is plenty for gameplay and runs at full framerate even on
    //  modest GPUs.  Set this to true if you want the full CUDA PIC loop
    //  (32^3 grid, 50k particles/species) running in lockstep with the
    //  0D model — note that even at 1 step/tick this drops the framerate
    //  to ~20 fps on a Tesla T4.
    pc.run_pic = true;
    PlasmaCoreBridge plasmacore(pc);
    HeliumConfig             hc;  HeliumSystem         heliumAsh(hc);
    ThermalHydraulicsConfig  tc2; tc2.coolant=CoolantType::FLiBe;
    ThermalHydraulics        thermalhydraulics(tc2);

    TurbineSystem       turbines;
    ElectricalGridSystem egrid;
    MoltenSaltSystem    salt;
    HeliumCoolingSystem heSystem;

    // ── Round 4 new modules ───────────────────────────────────────────────────
    HCDSystem                  hcd;            // NBI/ICRH/ECRH/LHCD actuators
    VacuumVesselSystem         vacuum;          // roughing+turbo+bakeout+boronization
    DisruptionMitigationSystem dm;             // MGI/SPI
    TritiumPlantSystem         tritiumplant;   // TES + detritiation

    // Initialise molten salt pumps
    for(auto& p:salt.saltState().hotleg)    {p.running=true;p.speed_frac=1.f;}
    for(auto& p:salt.saltState().coldleg)   {p.running=true;p.speed_frac=1.f;}
    for(auto& p:salt.saltState().blanket_circ){p.running=true;p.speed_frac=1.f;}

    // History buffers
    ScrollBuf h_pfus,h_te,h_ti,h_ne,h_q,h_elec,h_freq;

    // Plasma visualization (spatial diagnostics for the 3 new tabs)
    PlasmaViz plasmaviz;
    uint32_t viz_species_mask = (1u << PlasmaViz::SP_ELECTRON)
                              | (1u << PlasmaViz::SP_DEUTERIUM)
                              | (1u << PlasmaViz::SP_TRITIUM)
                              | (1u << PlasmaViz::SP_ALPHA);
    int viz_projection = 0;  // 0=R-Z, 1=top-down

    // UI state
    AlarmSystem alarms;
    bool paused=false,req_scram=false;
    float speed=1.f;
    int active_tab=0; // 0=overview, 1-4=turbines, 5=grid, 6=salt, 7=helium,
                      // 8=temperature, 9=fusion, 10=particles
    auto last_wall=std::chrono::steady_clock::now();

    // ── Main loop ─────────────────────────────────────────────────────────────
    for(bool running=true;running;){
        SDL_Event ev;
        while(SDL_PollEvent(&ev)){
            ImGui_ImplSDL2_ProcessEvent(&ev);
            if(ev.type==SDL_QUIT)running=false;
            if(ev.type==SDL_KEYDOWN){
                switch(ev.key.keysym.sym){
                    case SDLK_SPACE: paused=!paused;   break;
                    case SDLK_F1:    req_scram=true;   break;
                    case SDLK_ESCAPE:running=false;    break;
                    default: break;
                }
            }
        }

        // Physics
        //
        //  IMPORTANT: the main physics loop runs whenever the operator
        //  hasn't manually paused — *including* during a SCRAM.  Previously
        //  this was gated on `!state.cmd_scram` as well, which froze the
        //  simulation the moment a SCRAM tripped.  That had two nasty
        //  consequences:
        //
        //    1. The plasma never actually ramped down (the controller sets
        //       sp_plasma_current_MA = 0 on SCRAM, but the bridge needs to
        //       run to integrate that), so the operator couldn't see the
        //       shutdown happening — the sim was just frozen.
        //
        //    2. It was the root cause of the soft-lock: after the operator
        //       pressed RESET in the SCRAM popup, the ReactorState was
        //       wiped but the per-module state (MagnetSystem::dump_triggered_,
        //       ThermalHydraulics blanket temps, etc.) wasn't, and the very
        //       first tick after RESET re-tripped SCRAM.  With physics now
        //       running continuously, RESET can call reset() on every
        //       module and the next tick starts cleanly.
        //
        //  The SCRAM popup is now informational only — it shows while
        //  state.cmd_scram is true, but the simulation keeps running so
        //  the operator can watch the plasma ramp down.
        if(!paused){
            auto   now=std::chrono::steady_clock::now();
            double wdt=std::chrono::duration<double>(now-last_wall).count();
            last_wall=now;
            int ticks=std::clamp((int)(wdt*speed/SIM_DT_S),1,500);

            for(int i=0;i<ticks;i++){
                if(req_scram){state.cmd_scram=true;state.mode=ReactorMode::Emergency;req_scram=false;}
                control          .update(state,sim);
                magnets          .update(state,sim);
                fuel             .update(state,sim);
                // ── H&CD must run BEFORE plasmacore ───────────────────────────
                //  The plasmacore reads state.hcd_*_actual_MW to compute
                //  P_aux in the power balance.  Without running H&CD first,
                //  the bridge would see the previous tick's H&CD outputs
                //  (one-tick lag) — not a problem in steady state, but it
                //  causes a one-tick spike on H&CD on/off transitions.
                hcd              .update(state,sim);
                plasmacore       .update(state,SIM_DT_S);
                heliumAsh        .update(state,sim);
                thermalhydraulics.update(state,sim);
                // ── Vacuum vessel + tritium plant + DM ────────────────────────
                //  Vacuum runs every tick so the operator can see pressure
                //  evolve during pumpdown.  Tritium plant extracts T from
                //  the blanket.  Disruption mitigation consumes the
                //  cmd_disrupt_mitigation flag if the operator clicked FIRE.
                vacuum           .update(state,sim);
                tritiumplant     .update(state,sim);
                dm               .update(state,sim);

                // New module chain
                // Salt/steam coupling: read each SG's steam-side saturation
                // temperature and how much heat the secondary side can
                // currently absorb (from the end of the *previous* tick —
                // turbines.update() runs after salt.update() below). This
                // is what limits how much the molten salt actually cools in
                // each SG, rather than a fixed ΔT regardless of turbine state.
                float sg_tsat[4], sg_demand[4];
                for(int u=0;u<4;u++){
                    sg_tsat[u]   = turbines.unit(u).s.sg_steam_temp_K;
                    sg_demand[u] = turbines.unit(u).sgDemandFactor();
                }
                salt    .update(state,sim,state.blanket_heat_MW,sg_tsat,sg_demand);
                turbines.update(state,sim,state.grid_frequency_Hz,salt.sgHeatMW());
                egrid   .update(state,sim,turbines);
                heSystem.update(state,sim);

                sim.advance();
                //  Mirror the simulation clock into ReactorState.  The
                //  time_s field existed but was never advanced, so anything
                //  keyed off it (the OP SPACE breadcrumb trail sampler) saw
                //  a clock frozen at 0.0 and only ever fired once — the
                //  crosshair sat welded to its very first sample.
                state.time_s = sim.total_s;
            }
            h_pfus.push(state.fusion_power_MW);
            h_te  .push(state.electron_temp_keV);
            h_ti  .push(state.ion_temp_keV);
            h_ne  .push((float)(state.plasma_density_m3*1e-20));
            h_q   .push(state.q_safety);
            h_elec.push(state.net_electric_MW);
            h_freq.push(state.grid_frequency_Hz);

            // Update the spatial visualization data every tick.  This is
            // the only call site that touches the PlasmaViz — the three new
            // tabs are pure renderers that read from it.
            plasmaviz.update(state);

            // Record fusion history for the time-series tab
            constexpr float E_fus_J = 17.59e6f * 1.602176634e-19f;
            float rate_per_s = state.fusion_power_MW * 1e6f / E_fus_J;
            plasmaviz.recordFusionHistory(
                state.fusion_power_MW,
                state.alpha_power_MW,
                state.fusion_power_MW * (14.070f / 17.589f),
                rate_per_s);
        } else { last_wall=std::chrono::steady_clock::now(); }

        // ── Alarms (Round 4: descriptive cause + remediation) ───────────────────
        //  Each alarm now carries a "what went wrong" cause string and a "what
        //  to do" remediation string, plus a severity level.  The OVERVIEW tab
        //  and the new ALARMS tab display these in a structured layout
        //  instead of just a one-line label.
        {
            // Build the disruption cause string dynamically from the current
            // state so the operator sees the EXACT trigger (with values)
            // rather than a generic message.
            char cause_buf[256];
            char remed_buf[512];
            if (state.alarm_disruption) {
                const char* base_cause = disruptionCauseStr(state.disruption_cause);
                const char* base_remed = disruptionCauseRemediation(state.disruption_cause);
                snprintf(cause_buf, sizeof(cause_buf),
                    "%s.  Current: q_95=%.2f, β_N=%.2f, n_e=%.2e (Greenwald frac=%.2f), I_p=%.2f MA.",
                    base_cause, state.q_safety, state.beta_N,
                    state.plasma_density_m3,
                    state.plasma_density_m3 * 1e-20f /
                        std::max(0.001f, state.plasma_current_MA / (3.14159265f * 4.0f)),
                    state.plasma_current_MA);
                snprintf(remed_buf, sizeof(remed_buf), "%s", base_remed);
                alarms.trip("Disruption risk", state.alarm_disruption, sim.total_s,
                    cause_buf, remed_buf, AlarmSeverity::Critical);
            } else {
                alarms.trip("Disruption risk", false, sim.total_s);
            }

            // Magnet quench — coil exceeded critical temperature.
            snprintf(cause_buf, sizeof(cause_buf),
                "Superconducting coil temperature %.2f K exceeded the critical "
                "temperature %.1f K. The coil transitioned from superconducting "
                "to resistive state, dumping stored magnetic energy (%.1f GJ) "
                "into the dump resistor.  Coil current is decaying with τ = %.1f s.",
                state.magnet_temp_K, 18.0f, state.stored_energy_GJ, 10.0f);
            alarms.trip("Magnet quench", state.alarm_quench, sim.total_s,
                cause_buf,
                "Verify the QPS (Quench Protection System) fired correctly and "
                "the dump resistor absorbed the energy.  Inspect the coil for "
                "hotspots.  Do NOT re-energize the magnet until coil temp < 5 K "
                "and a full QPS diagnostic has been completed.",
                AlarmSeverity::Critical);

            // Overtemperature — first wall, divertor, or coolant.
            snprintf(cause_buf, sizeof(cause_buf),
                "First wall temp %.0f K (limit %d K), divertor temp %.0f K "
                "(limit %d K), coolant outlet %.0f K (limit %d K).  "
                "Excessive heat flux is exceeding the cooling capacity.",
                state.first_wall_temp_K, 2000,
                state.divertor_temp_K, 2200,
                state.coolant_outlet_temp_K, 900);
            alarms.trip("Overtemperature", state.alarm_overtemp, sim.total_s,
                cause_buf,
                "Reduce fusion power (lower I_p or T_e setpoint) to drop the "
                "heat flux.  Verify coolant flow rate (state.coolant_flow_kg_s) "
                "is at setpoint — if not, check molten salt pumps.  If divertor "
                "overtemp, increase impurity seeding (Ar/Ne) for radiative "
                "cooling of the scrape-off layer.",
                AlarmSeverity::Critical);

            // Low tritium — fuel store running out.
            snprintf(cause_buf, sizeof(cause_buf),
                "Tritium fuel inventory %.1f g is below the 10 g minimum.  "
                "The blanket TBR is currently %.3f (target 1.05) and TES is %s.  "
                "Tritium consumption rate is %.3f g/s at the current fusion power.",
                state.fuel_T_inventory_g, state.tbr_current,
                state.tritium_recovery_on ? "RUNNING" : "OFF",
                state.fusion_power_MW * 1e6f / (17.59e6f * 1.602176634e-19f)
                    * 5.00735588e-27f * 1e3f);
            alarms.trip("Low tritium", state.alarm_low_tritium, sim.total_s,
                cause_buf,
                "Start the Tritium Extraction System (TES) on the TRITIUM "
                "PLANT tab to recover T from the blanket.  If TES is already "
                "running, reduce fusion power to slow T consumption.  Request "
                "an external T resupply (+T button on the left panel) for "
                "immediate relief.",
                AlarmSeverity::Warning);

            // Grid underfrequency — turbine/grid frequency droop.
            snprintf(cause_buf, sizeof(cause_buf),
                "Grid frequency %.3f Hz is below the 49.0 Hz lower operating "
                "limit (nominal 50.0 Hz).  Total generation %.1f MW vs site "
                "load %.1f MW — deficit %.1f MW is being imported from the "
                "external grid.  If LOOP occurs the plant will lose offsite power.",
                egrid.grid().frequency_Hz,
                egrid.grid().total_generation_MW,
                egrid.grid().total_site_load_MW,
                egrid.grid().total_site_load_MW - egrid.grid().total_generation_MW);
            alarms.trip("Grid underfrequency",
                egrid.grid().underfrequency_alarm, sim.total_s,
                cause_buf,
                "Increase generator governor demand on the turbine tabs to "
                "raise generation.  Verify all 4 turbine breakers are closed "
                "and the units are online.  If a turbine has tripped, restart "
                "it.  If load exceeds generation capacity, shed non-essential "
                "site loads (toggle off on the ELEC GRID tab).",
                AlarmSeverity::Warning);

            // Cryostat vacuum loss — thermal insulation compromised.
            snprintf(cause_buf, sizeof(cause_buf),
                "Cryostat vacuum pressure %.2e Pa exceeds the 1e-2 Pa upper "
                "limit.  Heat leak to the 4.5 K cold mass has increased, "
                "raising the cryoplant load beyond its capacity.  Magnet "
                "temperature will rise toward the quench threshold.",
                heSystem.heState().cryostat.vacuum_pressure_Pa);
            alarms.trip("Cryostat vacuum loss",
                !heSystem.heState().cryostat.vacuum_ok, sim.total_s,
                cause_buf,
                "Start the cryostat roughing pump, then the turbo pump, on "
                "the HELIUM tab to recover vacuum.  Check for O-ring leaks "
                "on any recently opened cryostat penetration.  If the magnet "
                "temp rises above 6 K, prepare for a controlled SCRAM.",
                AlarmSeverity::Critical);

            // Magnet He high temp — cryoplant can't keep up.
            snprintf(cause_buf, sizeof(cause_buf),
                "Magnet cold mass temperature %.3f K exceeds the 5.5 K "
                "operating limit (target 4.5 K).  Cryoplant heat load is "
                "%.0f W vs refrigerator capacity.  Approaching the quench "
                "threshold of %.0f K.",
                heSystem.heState().magnet_circuit.cold_mass_temp_K,
                heSystem.heState().magnet_circuit.total_heat_load_W,
                15.0f);
            alarms.trip("Magnet He hi temp",
                heSystem.heState().magnet_circuit.hi_temp_alarm, sim.total_s,
                cause_buf,
                "Verify the cryoplant refrigerator is running (HELIUM tab).  "
                "If running, the heat load is too high — reduce magnet current "
                "ramp rate.  Check for abnormal eddy-current heating (rapid "
                "current changes).  If temp > 6 K, prepare for SCRAM.",
                AlarmSeverity::Critical);

            // Hot tank low level — salt inventory insufficient.
            snprintf(cause_buf, sizeof(cause_buf),
                "Hot tank level %.2f m is below the 2.0 m minimum (capacity "
                "16 m).  Salt flow to the steam generators will be interrupted "
                "if level drops further, causing turbine trips.",
                salt.saltState().hot_tank.level_m);
            alarms.trip("Hot tank lo level",
                salt.saltState().hot_tank.lo_level_alarm, sim.total_s,
                cause_buf,
                "Verify blanket circulation pumps (BCP-1, BCP-2) are running "
                "and at full speed (MOLTEN SALT tab).  Check for salt leaks "
                "in the hot leg.  Reduce turbine governor demand if salt "
                "flow can't keep up with steam demand.",
                AlarmSeverity::Warning);

            // Hot tank freeze — FLiBe has frozen in the tank.
            snprintf(cause_buf, sizeof(cause_buf),
                "Hot tank temperature %.0f K is below the FLiBe freeze point "
                "733 K.  Solidified salt will block flow paths and damage "
                "pumps.  Trace heating must be activated to recover.",
                salt.saltState().hot_tank.temp_K);
            alarms.trip("Hot tank freeze",
                salt.saltState().hot_tank.lo_temp_alarm, sim.total_s,
                cause_buf,
                "Activate trace heating on the affected tank (MOLTEN SALT "
                "tab) and wait for temperature to recover above 800 K before "
                "starting any salt pumps.  Do NOT pump solid salt — wait "
                "for it to melt.  Reduce fusion power to allow the salt to "
                "reheat from blanket flow.",
                AlarmSeverity::Critical);

            // Vessel vacuum loss — can't initiate plasma.
            snprintf(cause_buf, sizeof(cause_buf),
                "Vacuum vessel pressure %.2e Pa exceeds the 1e-3 Pa initiation "
                "threshold.  Plasma cannot be (re)initiated until base vacuum "
                "is recovered.  Wall outgassing factor = %.3f (1 = unbaked).",
                state.vessel_pressure_Pa,
                1.0f);  // simplification — VacuumVessel tracks this internally
            alarms.trip("Vessel vacuum loss", state.alarm_vacuum_loss, sim.total_s,
                cause_buf,
                "Start the vessel roughing pump, wait for pressure < 1 Pa, "
                "then start the turbo pump (VACUUM & DM tab).  If pressure "
                "plateaus above 1e-3 Pa, run wall bakeout (423 K for 5 min) "
                "to desorb water.  Check for vessel breach (loss-of-vacuum "
                "accident) if pressure rises rapidly with pumps on.",
                AlarmSeverity::Warning);

            // Vessel breach — LOVA (loss-of-vacuum accident).
            snprintf(cause_buf, sizeof(cause_buf),
                "Vacuum vessel breach detected!  Pressure is rising rapidly "
                "toward atmospheric (%.0f Pa and climbing).  This is a "
                "loss-of-vacuum accident (LOVA) — air ingress will quench "
                "the plasma and may mobilize tritium.",
                state.vessel_pressure_Pa);
            alarms.trip("Vessel breach", state.alarm_vessel_breach, sim.total_s,
                cause_buf,
                "SCRAM the reactor immediately if not already (F1).  Isolate "
                "the tritium plant and start detritiation.  Evacuate the "
                "reactor building.  Do NOT attempt to recover vacuum until "
                "the breach location has been identified and isolated.",
                AlarmSeverity::Critical);

            // H&CD fault — one of the four H&CD systems has tripped.
            snprintf(cause_buf, sizeof(cause_buf),
                "An H&CD system fault is active (NBI: %s, ICRH: %s, "
                "ECRH: %s, LHCD: %s).  Total aux heating is now %.1f MW "
                "(was %.1f MW setpoint).",
                hcd.nbiFault()  ? "FAULT" : "ok",
                hcd.icrhFault() ? "FAULT" : "ok",
                hcd.ecrhFault() ? "FAULT" : "ok",
                hcd.lhcdFault() ? "FAULT" : "ok",
                state.hcd_nbi_actual_MW + state.hcd_icrh_actual_MW
                    + state.hcd_ecrh_actual_MW + state.hcd_lhcd_actual_MW,
                state.sp_aux_heat_MW);
            alarms.trip("H&CD system fault", state.alarm_aux_heat_fault, sim.total_s,
                cause_buf,
                "Open the H&CD tab to identify which system faulted.  Common "
                "causes: gyrotron overtemp (ECRH), antenna VSWR (ICRH), "
                "beamline valve closed (NBI), klystron HV trip (LHCD).  "
                "Reset the fault on the H&CD tab once the root cause is "
                "fixed.  Plasma temperature may droop while the fault is "
                "active — consider reducing T_e setpoint.",
                AlarmSeverity::Warning);

            // Runaway electron beam / wall strike.
            snprintf(cause_buf, sizeof(cause_buf),
                "Runaway electron activity: beam current %.2f MA, wall "
                "strikes this campaign: %d.  An unmitigated fast current "
                "quench avalanche-converted plasma current into a "
                "relativistic electron beam.",
                state.re_beam_MA, state.re_strike_count);
            alarms.trip("Runaway electrons", state.alarm_runaway, sim.total_s,
                cause_buf,
                "The RE beam forms during UNMITIGATED disruptions above "
                "~4 MA.  Prevention: keep MGI/SPI armed (VACUUM & DM tab) "
                "with auto-MGI enabled — the injected impurities raise the "
                "density and collisionally kill the avalanche.  After a "
                "wall strike, check the first-wall temperature before "
                "re-initiating.",
                AlarmSeverity::Critical);

            // Turbine trips (per unit)
            for(int i=0;i<4;i++){
                char lbl[32];
                snprintf(lbl,32,"Turbine %d trip",i+1);
                snprintf(cause_buf, sizeof(cause_buf),
                    "Turbine %d has tripped.  State: %s.  Rotor speed %.0f RPM, "
                    "SG pressure %.2f MPa, generator breaker %s.  See the "
                    "turbine tab for trip cause (overspeed / low vacuum / "
                    "low SG level / differential current).",
                    i+1, turbStateStr(turbines.unit(i).s.state),
                    turbines.unit(i).s.rpm,
                    turbines.unit(i).s.sg_pressure_MPa,
                    turbines.unit(i).s.generator.breaker_closed ? "CLOSED" : "OPEN");
                alarms.trip(lbl, turbines.unit(i).s.alarm_turb_trip, sim.total_s,
                    cause_buf,
                    "Reset the turbine on its tab (RESET button).  Verify "
                    "the trip cause has cleared before restarting.  Restart "
                    "sequence: START UNIT → wait for synch → CLOSE BREAKER.  "
                    "If the trip was overspeed, inspect the governor.",
                    AlarmSeverity::Critical);
            }
        }

        // ── Render ─────────────────────────────────────────────────────────────
        ImGui_ImplOpenGL3_NewFrame();ImGui_ImplSDL2_NewFrame();ImGui::NewFrame();

        RenderStatusBar(state,sim,paused,alarms,speed);

        float lw=io.DisplaySize.x*LEFT_W_FRAC;
        float rw=io.DisplaySize.x-lw;
        float ch=io.DisplaySize.y-STATUS_H;
        float ry=STATUS_H;

        // Left panel
        RenderLeftPanel(state,sim,paused,req_scram,speed,
                        fuel,turbines,salt,heSystem,alarms);

        // Right area: tab bar header
        BeginTiled("##main_tabs",lw,ry,rw,ch,ImGuiWindowFlags_NoTitleBar);

        // ── Round 4: expanded tab list ─────────────────────────────────────────
        //  Original 11 tabs + 5 new tabs:
        //    11 = H&CD         — NBI/ICRH/ECRH/LHCD actuators
        //    12 = VACUUM & DM  — vessel pumping + disruption mitigation
        //    13 = TRITIUM      — TES / detritiation / accountancy
        //    14 = OPERATIONS   — startup/shutdown procedures with checklists
        //    15 = ALARMS       — full alarm list with cause + remediation
        //    16 = OP SPACE     — live POPCON operating-space diagram
        static const char* tab_labels[]={
            "OVERVIEW","TURBINE 1","TURBINE 2","TURBINE 3","TURBINE 4",
            "ELEC GRID","MOLTEN SALT","HELIUM",
            "TEMP MAP","FUSION HIST","PARTICLES",
            "H&CD","VACUUM&DM","TRITIUM","OPERATIONS","ALARMS","OP SPACE"
        };
        constexpr int N_TABS = 17;
        if(ImGui::BeginTabBar("##main")){
            for(int i=0;i<N_TABS;i++){
                // ── Per-tab style-color stack management ─────────────────────────
                //  Each branch below may push 1 or 2 colors onto ImGui's style
                //  stack.  We track the count in `col_push_count` and pop
                //  exactly that many at the bottom of the loop — using a
                //  bool + PopStyleColor(2) was the source of the
                //  "[imgui-error] Missing PopStyleColor()" report: when i==15
                //  (ALARMS tab) had unacked critical alarms, BOTH the i>=11
                //  branch (pushes 2) and the i==15 branch (pushes 1) fired,
                //  pushing 3 colors total, but the bool-based pop only
                //  removed 2 — leaving 1 orphaned color on the stack.
                int col_push_count = 0;
                // Colour turbine tabs by state
                if(i>=1&&i<=4){
                    TurbineState st=turbines.unit(i-1).s.state;
                    if(st==TurbineState::Online){
                        ImGui::PushStyleColor(ImGuiCol_Tab,Col::GREEN_DARK);
                        ImGui::PushStyleColor(ImGuiCol_TabActive,{0.08f,0.30f,0.10f,1.f});
                        col_push_count += 2;
                    } else if(st==TurbineState::Tripped||st==TurbineState::Tripping){
                        ImGui::PushStyleColor(ImGuiCol_Tab,{0.30f,0.04f,0.04f,1.f});
                        ImGui::PushStyleColor(ImGuiCol_TabActive,{0.45f,0.06f,0.06f,1.f});
                        col_push_count += 2;
                    }
                }
                // Highlight the three plasma-viz tabs (and OP SPACE) in cyan
                if((i>=8&&i<=10)||i==16){
                    ImGui::PushStyleColor(ImGuiCol_TabHovered,{0.10f,0.35f,0.40f,1.f});
                    ImGui::PushStyleColor(ImGuiCol_TabActive,{0.06f,0.25f,0.30f,1.f});
                    col_push_count += 2;
                }
                // Round 4 new tabs in amber (OP SPACE keeps the cyan above)
                if(i>=11&&i<=15){
                    ImGui::PushStyleColor(ImGuiCol_TabHovered,{0.40f,0.30f,0.10f,1.f});
                    ImGui::PushStyleColor(ImGuiCol_TabActive,{0.30f,0.22f,0.08f,1.f});
                    col_push_count += 2;
                }
                // ALARMS tab blinks red when there are unacked critical alarms.
                // This adds ONE more color on top of the amber from the i>=11
                // branch — that's intentional (the red overrides the amber
                // base colour).  We just need to make sure we pop it too.
                if(i==15 && alarms.countCritical()>0){
                    float bl=fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f?1.f:0.f;
                    ImGui::PushStyleColor(ImGuiCol_Tab,{0.40f*bl,0.04f,0.04f,1.f});
                    col_push_count += 1;
                }
                if(ImGui::BeginTabItem(tab_labels[i])){
                    active_tab=i;
                    ImGui::EndTabItem();
                }
                if(col_push_count>0){
                    ImGui::PopStyleColor(col_push_count);
                }
            }
            ImGui::EndTabBar();
        }

        // Tab content
        switch(active_tab){
            case 0: TabOverview(state,h_pfus,h_te,h_ti,h_ne,h_q,alarms,sim); break;
            case 1: case 2: case 3: case 4:
                TabTurbine(turbines.unit(active_tab-1),egrid,active_tab-1); break;
            case 5: TabGrid(egrid,turbines); break;
            case 6: TabMoltenSalt(salt,turbines); break;
            case 7: TabHelium(heSystem);     break;
            // ── Three plasma visualization tabs ──────────────────────────────
            case 8:  TabTemperature(state, plasmaviz); break;
            case 9:  TabFusionReactions(state, plasmaviz); break;
            case 10: TabParticleDistribution(state, plasmaviz,
                                              viz_species_mask, viz_projection); break;
            // ── Round 4 new tabs ─────────────────────────────────────────────
            case 11: TabHCD(state, hcd); break;
            case 12: TabVacuumDM(state, vacuum, dm); break;
            case 13: TabTritiumPlant(state, tritiumplant); break;
            case 14: TabOperations(state, sim); break;
            case 15: TabAlarms(alarms); break;
            // ── OP SPACE: live POPCON operating-space diagram ────────────────
            case 16: TabOpSpace(state); break;
        }
        ImGui::End();

        // SCRAM modal — informational, not blocking.
        //
        //  The simulation continues to run while this popup is up (see the
        //  main-loop comment above) so the operator can watch the plasma
        //  ramp down in real time on the OVERVIEW tab.  The two buttons
        //  below let the operator either do a full cold restart (RESET)
        //  or attempt to recover the discharge (ACKNOWLEDGE).
        if(state.cmd_scram){
            ImGui::SetNextWindowSize({450.f,220.f},ImGuiCond_Always);
            ImGui::SetNextWindowPos(io.DisplaySize,ImGuiCond_Always,{1.f,1.f}); // bottom-right
            ImGui::PushStyleColor(ImGuiCol_WindowBg,{0.14f,0.02f,0.02f,0.97f});
            ImGui::PushStyleColor(ImGuiCol_Border,   Col::RED);
            ImGui::Begin("##scram",nullptr,ImGuiWindowFlags_NoTitleBar|
                ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar);
            float bl=fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f?1.f:0.f;
            ImGui::TextColored({1.f,bl,bl,1.f},"  ⚡  SCRAM ACTIVATED  ⚡");
            ImGui::TextColored(Col::AMBER,"  Plasma ramp-down in progress.");
            ImGui::TextColored(Col::AMBER,"  All turbines receiving trip signal.");
            // Live status so the operator can see the ramp-down happening
            // rather than wondering whether the sim is frozen.
            ImGui::TextColored(Col::WHITE,   "  I_p:    %.2f MA",  state.plasma_current_MA);
            ImGui::TextColored(Col::WHITE,   "  P_fus:  %.1f MW",  state.fusion_power_MW);
            ImGui::TextColored(Col::WHITE,   "  T_e:    %.2f keV", state.electron_temp_keV);
            ImGui::Separator();
            if(ImGui::Button("RESET — COLD RESTART",{434.f,28.f})){
                state={};
                state.sp_plasma_current_MA=15.f;state.sp_electron_temp_keV=20.f;
                state.sp_density_m3=1e20f;state.sp_B_toroidal_T=5.3f;
                state.sp_fuel_rate=0.5f;state.sp_coolant_flow=0.8f;state.D_T_ratio=1.f;
                state.mode=ReactorMode::Startup;state.plasma_status=PlasmaStatus::Cold;
                // Initialize plasma-shape setpoints to ITER nominal
                state.sp_kappa = 1.7f; state.sp_delta = 0.33f;
                sim={};sim.dt_s=SIM_DT_S;alarms.log.clear();alarms.any_unacked=false;
                // ── CRITICAL: clear the latched SCRAM in the control system ──────
                // Without this, the next control.update() call sees the stale
                // scram_latched_ flag and immediately re-sets state.cmd_scram,
                // trapping the simulation in a permanent SCRAM loop.
                control.resetScramLatch();
                // ── CRITICAL: reset EVERY physics module's internal state ──────
                // The ReactorState wipe above only clears the shared summary
                // fields — each module (magnets, thermal hydraulics, helium,
                // molten salt, turbines, etc.) carries its own private state
                // (coil temps, blanket temps, divertor tile temp, tank levels,
                // PIC plasma energy, dump latches, ...) that survives across
                // the ReactorState wipe.  Without resetting those too, the
                // next tick's update() calls write the stale hot temps back
                // into the freshly-zeroed ReactorState and the control
                // system re-trips SCRAM within a handful of ticks — the
                // classic "RESET only advances 16 ticks then SCRAMs again"
                // soft-lock the user reported.
                magnets          .reset();
                fuel             .reset();
                plasmacore       .reset();
                heliumAsh        .reset();
                thermalhydraulics.reset();
                turbines         .reset();
                egrid            .reset();
                salt             .reset();
                heSystem         .reset();
                // ── Round 4 new modules: also reset on cold restart ──────────
                hcd              .reset();
                vacuum           .reset();
                dm               .reset();
                tritiumplant     .reset();
            }

            ImGui::Spacing();
            ImGui::TextColored(Col::GREY,"  — or —");
            ImGui::Spacing();
            // Soft resume: clear the SCRAM latch but keep the current plasma
            // state (so the operator can try to recover the discharge rather
            // than starting from cold).  The control system will re-evaluate
            // the SCRAM triggers on the next tick — if the underlying cause
            // (overheating, quench, etc.) is still present, SCRAM will
            // re-trip immediately, which is the correct safety behaviour.
            if(ImGui::Button("ACKNOWLEDGE — RESUME WITHOUT RESET",{434.f,28.f})){
                state.cmd_scram = false;
                state.mode = ReactorMode::SteadyState;
                control.resetScramLatch();
                // ── Clear the magnet quench latch ──────────────────────────
                //  When SCRAM fired, MagnetSystem::triggerQuenchDump() set
                //  dump_triggered_=true, which makes state.quench_detected
                //  =true every tick.  That quench_detected flag is a SCRAM
                //  trigger, so without clearing it here, the next tick's
                //  runScramLogic immediately re-latches the SCRAM — the
                //  operator presses ACKNOWLEDGE and nothing happens.
                //  clearQuenchLatch() clears the dump latch and per-coil
                //  quenched flags WITHOUT resetting coil currents/temps,
                //  so the operator can resume from the current state.
                magnets.clearQuenchLatch();
                state.quench_detected = false;
                state.alarm_quench = false;
                // ── Restore default setpoints ──────────────────────────────
                //  The SCRAM zeroed ALL setpoints (sp_plasma_current_MA,
                //  sp_electron_temp_keV, sp_density_m3, sp_fuel_rate) via
                //  ControlSystem::update when scram_latched_ was true.
                //  Without restoring them here, the plasma continues to
                //  ramp toward 0 even after the operator acknowledges the
                //  SCRAM — the mode says "STEADY STATE" but the plasma is
                //  still shutting down.  The operator perceives this as
                //  "Ramp Down doesn't reset when SCRAM is reset."
                //
                //  We restore to the default ITER-class setpoints (same as
                //  what RESET uses).  The plasma state itself (I_p, T_e,
                //  ne) is NOT reset — only the setpoints are restored, so
                //  the control system can start driving the plasma back up
                //  from wherever it currently is.
                state.sp_plasma_current_MA = 15.f;
                state.sp_electron_temp_keV = 20.f;
                state.sp_density_m3        = 1e20f;
                state.sp_fuel_rate         = 0.5f;
                state.sp_coolant_flow      = 0.8f;
                // ── Round 4: restore plasma-shape setpoints ──────────────────
                state.sp_kappa = 1.7f; state.sp_delta = 0.33f;
                // If the plasma has fully quenched, go back to Startup so
                // the operator can re-initiate via the INITIATE PLASMA
                // button (which requires plasma_status == Cold).
                if (state.plasma_status == PlasmaStatus::Quenched ||
                    state.plasma_status == PlasmaStatus::Disrupting) {
                    state.plasma_status = PlasmaStatus::Cold;
                    state.mode = ReactorMode::Startup;
                }
                alarms.ackAll();
            }
            ImGui::End();ImGui::PopStyleColor(2);
        }

        ImGui::Render();
        int dw,dh;SDL_GetWindowSize(win,&dw,&dh);
        glViewport(0,0,dw,dh);
        glClearColor(Col::BG.x,Col::BG.y,Col::BG.z,1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(win);
    }

    ImGui_ImplOpenGL3_Shutdown();ImGui_ImplSDL2_Shutdown();ImGui::DestroyContext();
    SDL_GL_DeleteContext(gl);SDL_DestroyWindow(win);SDL_Quit();
    return 0;
}