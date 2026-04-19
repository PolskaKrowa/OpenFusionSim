# src/ThermalHydraulics/CMakeLists.txt
# Breeding blanket, primary coolant loop, and first-wall thermal model.

add_library(thermalhydraulics STATIC
    ThermalHydraulics.cpp
    MoltenSaltPhysics.cpp
)

target_include_directories(thermalhydraulics
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(thermalhydraulics PRIVATE -Wall -Wextra)