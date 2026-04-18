# src/SteamPower/CMakeLists.txt
# Rankine cycle: steam generator, turbine, condenser, generator.

add_library(steampower STATIC
    SteamPower.cpp
)

target_include_directories(steampower
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(steampower PRIVATE -Wall -Wextra)