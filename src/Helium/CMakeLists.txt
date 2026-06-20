# src/Helium/CMakeLists.txt
# Helium ash management: divertor, exhaust pumping, impurity control.

add_library(helium STATIC
    Helium.cpp
    HeliumCoolingPhysics.cpp
)

target_include_directories(helium
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(helium PRIVATE -Wall -Wextra)