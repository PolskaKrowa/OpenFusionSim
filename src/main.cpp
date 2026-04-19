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

struct AlarmEntry{std::string msg;double t;bool active,acked;};
struct AlarmSystem{
    std::vector<AlarmEntry> log; bool any_unacked=false;
    void trip(const std::string& m,bool active,double t){
        for(auto& a:log){if(a.msg==m){if(active&&!a.active){a.acked=false;any_unacked=true;}a.active=active;return;}}
        if(active){log.push_back({m,t,true,false});any_unacked=true;}
    }
    void ackAll(){for(auto& a:log)a.acked=true;any_unacked=false;}
    void clearOld(){log.erase(std::remove_if(log.begin(),log.end(),[](const AlarmEntry&a){return!a.active&&a.acked;}),log.end());}
};

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
    ImGui::TextColored(Col::GREEN,  "P_fus %.0f MW",  s.fusion_power_MW);   ImGui::SameLine(490.f);
    ImGui::TextColored(Col::YELLOW, "P_net %.0f MW",  s.net_electric_MW);   ImGui::SameLine(630.f);
    ImGui::TextColored(Col::CYAN,   "Q=%.2f",         s.Q_scientific);      ImGui::SameLine(730.f);
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
    if(ImGui::Button(paused?"▶ RESUME":"■■ PAUSE",{pw,0.f}))paused=!paused;
    ImGui::TextColored(Col::GREY," SPACE=pause  F1=scram  ESC=quit");

    Hdr("REACTOR");
    bool cold=(s.plasma_status==PlasmaStatus::Cold);
    if(cold){
        if(GreenBtn("▶ INITIATE PLASMA",{pw,26.f})){
            s.plasma_status=PlasmaStatus::Initiating;
            s.plasma_current_MA=0.5f;s.electron_temp_keV=1.f;
            s.plasma_density_m3=5e18f;s.mode=ReactorMode::Startup;
        }
    } else {
        if(RedBtn("■ CONTROLLED SHUTDOWN",{pw,24.f})){
            s.mode=ReactorMode::Rampdown;
            s.sp_plasma_current_MA=0.f;s.sp_fuel_rate=0.f;
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

    Hdr("FUEL");
    float pct = s.sp_fuel_rate * 100.f;
    if (ImGui::SliderFloat("##gp", &pct, 0.f, 100.f, "Gas: %.0f%%"))
        s.sp_fuel_rate = pct / 100.f;
    ImGui::SetNextItemWidth(pw);ImGui::SliderFloat("##pf",&s.pellet_frequency_Hz,0.f,10.f,"Pellets: %.1f Hz");
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
                         ScrollBuf& h_pfus,ScrollBuf& h_te,ScrollBuf& h_ne,ScrollBuf& h_q,
                         AlarmSystem& alm,const SimTime& t)
{
    float pw=ImGui::GetContentRegionAvail().x;
    if(ImGui::BeginChild("##ov_left",{pw*0.55f,0.f},false)){
        Hdr("CORE PLASMA");
        ImVec4 qc=(s.q_safety<2.f)?Col::RED:(s.q_safety<2.5f)?Col::AMBER:Col::GREEN;
        Row("Plasma Current",   s.plasma_current_MA,   "%.2f","MA");
        Row("Electron Temp",    s.electron_temp_keV,   "%.1f","keV",Col::AMBER);
        Row("Ion Temp",         s.plasma_temp_keV,     "%.1f","keV",Col::AMBER);
        Row("Density",          s.plasma_density_m3,   "%.2e","m-3",Col::CYAN);
        Row("Safety Factor q95",s.q_safety,            "%.3f","",  qc);
        Row("Beta",             s.beta*100.f,          "%.3f","%");
        Row("He Ash",           s.helium_fraction*100.f,"%.2f","%");

        Hdr("POWER BALANCE");
        Row("Fusion Power",   s.fusion_power_MW,  "%.1f","MW");
        Row("Alpha Heating",  s.alpha_power_MW,   "%.1f","MW",Col::AMBER);
        Row("Radiated Power", s.radiated_power_MW,"%.1f","MW",Col::GREY);
        Row("Q scientific",   s.Q_scientific,     "%.3f","");
        Row("Blanket Heat",   s.blanket_heat_MW,  "%.1f","MW",Col::AMBER);
        Row("Net Electric",   s.net_electric_MW,  "%.1f","MW",Col::YELLOW);
        Row("Grid Frequency", s.grid_frequency_Hz,"%.3f","Hz",
            s.grid_frequency_Hz<49.f?Col::AMBER:Col::GREEN);

        Hdr("MAGNETS");
        Row("B_T",          s.B_toroidal_T,    "%.3f","T");
        Row("Coil Current", s.coil_current_kA, "%.1f","kA");
        Row("Magnet Temp",  s.magnet_temp_K,   "%.2f","K",
            s.magnet_temp_K>6.f?Col::AMBER:Col::GREEN);
        Lamp(" QUENCH ",s.quench_detected);

        Hdr("PLASMA TRENDS");
        ImGui::TextColored(Col::GREY," Fusion Power [MW]");Plot("##op",h_pfus,pw*.55f,40.f,Col::GREEN);
        ImGui::TextColored(Col::GREY," Te [keV]");         Plot("##ote",h_te,pw*.55f,40.f,Col::AMBER);
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
            ImVec4 ac=a.acked?Col::AMBER:ImVec4{1.f,bl,bl,1.f};
            ImGui::TextColored(ac," ▲ [t=%.1fs] %s%s",a.t,a.msg.c_str(),a.acked?" (ack)":"");
        }
        if(!any)ImGui::TextColored(Col::GREEN,"  No active alarms");
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
static void TabMoltenSalt(MoltenSaltSystem& salt)
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
        for(int i=0;i<2;i++){
            auto& p=s.blanket_circ[i];
            ImGui::PushID(i+2000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," BCP-%d: %s  %.0f kg/s",
                               i+1,p.running?"RUN":"STOP",p.flow_kg_s);
            float bw=(cw-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f}))p.running=false;
            ImGui::SetNextItemWidth(cw-16.f);
            ImGui::SliderFloat("##bcps",&p.speed_frac,0.f,1.f,"Speed: %.0f%%");
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
            ImGui::TextColored(Col::ORANGE," SG%d: %.0f→%.0f K  Q=%.0f MW",
                i+1,s.sg_salt_inlet_K[i],s.sg_salt_outlet_K[i],s.sg_heat_MW[i]);
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
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f}))p.running=false;
            ImGui::SetNextItemWidth(pw*.34f-16.f);
            ImGui::SliderFloat("##hls",&p.speed_frac,0.f,1.f,"%.0f%%");
            ImGui::PopID();
        }

        Hdr("COLDLEG PUMPS");
        for(int i=0;i<4;i++){
            auto& p=s.coldleg[i];
            ImGui::PushID(i+4000);
            ImGui::TextColored(p.running?Col::GREEN:Col::GREY," CL-%d (SG%d return)",i+1,i+1);
            ImGui::TextColored(Col::GREY,"   %.0f kg/s",p.flow_kg_s);
            float bw=(pw*.34f-16.f)*.5f;
            if(GreenBtn("Start",{bw,0.f})){p.running=true;p.speed_frac=1.f;}ImGui::SameLine();
            if(RedBtn("Stop",{bw,0.f}))p.running=false;
            ImGui::SetNextItemWidth(pw*.34f-16.f);
            ImGui::SliderFloat("##cls",&p.speed_frac,0.f,1.f,"%.0f%%");
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
            ImGui::SliderFloat("##rps",&p.speed_frac,0.f,1.f,"Speed: %.0f%%");
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

    ControlConfig cc; ControlSystem  control(cc);
    MagnetConfig  mc; MagnetSystem   magnets(mc);
    FuelConfig    fc; FuelSystem     fuel(fc);
    PlasmaConfig  pc; pc.Nx=64;pc.Ny=64;pc.Nz=64;pc.pic_dt=1e-12f;
    PlasmaCoreBridge plasmacore(pc);
    HeliumConfig             hc;  HeliumSystem         heliumAsh(hc);
    ThermalHydraulicsConfig  tc2; tc2.coolant=CoolantType::FLiBe;
    ThermalHydraulics        thermalhydraulics(tc2);

    TurbineSystem       turbines;
    ElectricalGridSystem egrid;
    MoltenSaltSystem    salt;
    HeliumCoolingSystem heSystem;

    // Initialise molten salt pumps
    for(auto& p:salt.saltState().hotleg)    {p.running=true;p.speed_frac=1.f;}
    for(auto& p:salt.saltState().coldleg)   {p.running=true;p.speed_frac=1.f;}
    for(auto& p:salt.saltState().blanket_circ){p.running=true;p.speed_frac=1.f;}

    // History buffers
    ScrollBuf h_pfus,h_te,h_ne,h_q,h_elec,h_freq;

    // UI state
    AlarmSystem alarms;
    bool paused=false,req_scram=false;
    float speed=1.f;
    int active_tab=0; // 0=overview, 1-4=turbines, 5=grid, 6=salt, 7=helium
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
        if(!paused&&!state.cmd_scram){
            auto   now=std::chrono::steady_clock::now();
            double wdt=std::chrono::duration<double>(now-last_wall).count();
            last_wall=now;
            int ticks=std::clamp((int)(wdt*speed/SIM_DT_S),1,500);

            for(int i=0;i<ticks;i++){
                if(req_scram){state.cmd_scram=true;state.mode=ReactorMode::Emergency;req_scram=false;}
                control          .update(state,sim);
                magnets          .update(state,sim);
                fuel             .update(state,sim);
                plasmacore       .update(state,SIM_DT_S);
                heliumAsh        .update(state,sim);
                thermalhydraulics.update(state,sim);

                // New module chain
                salt    .update(state,sim,state.blanket_heat_MW);
                turbines.update(state,sim,state.grid_frequency_Hz,salt.sgHeatMW());
                egrid   .update(state,sim,turbines);
                heSystem.update(state,sim);

                sim.advance();
            }
            h_pfus.push(state.fusion_power_MW);
            h_te  .push(state.electron_temp_keV);
            h_ne  .push((float)(state.plasma_density_m3*1e-20));
            h_q   .push(state.q_safety);
            h_elec.push(state.net_electric_MW);
            h_freq.push(state.grid_frequency_Hz);
        } else { last_wall=std::chrono::steady_clock::now(); }

        // Alarms
        alarms.trip("Disruption risk",    state.alarm_disruption,    sim.total_s);
        alarms.trip("Magnet quench",      state.alarm_quench,        sim.total_s);
        alarms.trip("Overtemperature",    state.alarm_overtemp,      sim.total_s);
        alarms.trip("Low tritium",        state.alarm_low_tritium,   sim.total_s);
        alarms.trip("Grid underfrequency",egrid.grid().underfrequency_alarm, sim.total_s);
        alarms.trip("Cryostat vacuum loss",!heSystem.heState().cryostat.vacuum_ok,sim.total_s);
        alarms.trip("Magnet He hi temp",  heSystem.heState().magnet_circuit.hi_temp_alarm,sim.total_s);
        alarms.trip("Hot tank lo level",  salt.saltState().hot_tank.lo_level_alarm,sim.total_s);
        alarms.trip("Hot tank freeze",    salt.saltState().hot_tank.lo_temp_alarm, sim.total_s);
        for(int i=0;i<4;i++){
            char lbl[32];
            snprintf(lbl,32,"Turbine %d trip",i+1);
            alarms.trip(lbl,turbines.unit(i).s.alarm_turb_trip,sim.total_s);
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

        static const char* tab_labels[]={"OVERVIEW","TURBINE 1","TURBINE 2","TURBINE 3","TURBINE 4","ELEC GRID","MOLTEN SALT","HELIUM"};
        if(ImGui::BeginTabBar("##main")){
            for(int i=0;i<8;i++){
                // Colour turbine tabs by state
                bool col_pushed=false;
                if(i>=1&&i<=4){
                    TurbineState st=turbines.unit(i-1).s.state;
                    if(st==TurbineState::Online){
                        ImGui::PushStyleColor(ImGuiCol_Tab,Col::GREEN_DARK);
                        ImGui::PushStyleColor(ImGuiCol_TabActive,{0.08f,0.30f,0.10f,1.f});
                        col_pushed=true;
                    } else if(st==TurbineState::Tripped||st==TurbineState::Tripping){
                        ImGui::PushStyleColor(ImGuiCol_Tab,{0.30f,0.04f,0.04f,1.f});
                        ImGui::PushStyleColor(ImGuiCol_TabActive,{0.45f,0.06f,0.06f,1.f});
                        col_pushed=true;
                    }
                }
                if(ImGui::BeginTabItem(tab_labels[i])){
                    active_tab=i;
                    ImGui::EndTabItem();
                }
                if(col_pushed)ImGui::PopStyleColor(2);
            }
            ImGui::EndTabBar();
        }

        // Tab content
        switch(active_tab){
            case 0: TabOverview(state,h_pfus,h_te,h_ne,h_q,alarms,sim); break;
            case 1: case 2: case 3: case 4:
                TabTurbine(turbines.unit(active_tab-1),egrid,active_tab-1); break;
            case 5: TabGrid(egrid,turbines); break;
            case 6: TabMoltenSalt(salt);     break;
            case 7: TabHelium(heSystem);     break;
        }
        ImGui::End();

        // SCRAM modal
        if(state.cmd_scram){
            ImGui::SetNextWindowSize({450.f,200.f},ImGuiCond_Always);
            ImGui::SetNextWindowPos(io.DisplaySize,ImGuiCond_Always,{1.f,1.f}); // bottom-right
            ImGui::PushStyleColor(ImGuiCol_WindowBg,{0.14f,0.02f,0.02f,0.97f});
            ImGui::PushStyleColor(ImGuiCol_Border,   Col::RED);
            ImGui::Begin("##scram",nullptr,ImGuiWindowFlags_NoTitleBar|
                ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar);
            float bl=fmodf((float)ImGui::GetTime()*2.f,1.f)>0.5f?1.f:0.f;
            ImGui::TextColored({1.f,bl,bl,1.f},"  ⚡  SCRAM ACTIVATED  ⚡");
            ImGui::TextColored(Col::AMBER,"  Plasma ramp-down in progress.");
            ImGui::TextColored(Col::AMBER,"  All turbines receiving trip signal.");
            ImGui::Separator();
            if(ImGui::Button("RESET — COLD RESTART",{434.f,28.f})){
                state={};
                state.sp_plasma_current_MA=15.f;state.sp_electron_temp_keV=20.f;
                state.sp_density_m3=1e20f;state.sp_B_toroidal_T=5.3f;
                state.sp_fuel_rate=0.5f;state.sp_coolant_flow=0.8f;state.D_T_ratio=1.f;
                state.mode=ReactorMode::Startup;state.plasma_status=PlasmaStatus::Cold;
                sim={};sim.dt_s=SIM_DT_S;alarms.log.clear();alarms.any_unacked=false;
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