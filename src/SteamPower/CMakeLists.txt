# src/SteamPower/CMakeLists.txt
# Rankine cycle: steam generator, turbine, condenser, generator.

add_library(steampower STATIC
    SteamPower.cpp
    SteamCyclePhysics.cpp
    iapws_if97_coefficients.cpp    # NEW: full IAPWS-IF97 coefficient tables
)

target_include_directories(steampower
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(steampower PRIVATE -Wall -Wextra)