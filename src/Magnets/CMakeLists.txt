# src/Magnets/CMakeLists.txt
# Superconducting magnet systems: TF coils, PF coils, CS solenoid.
# Pure C++ — no CUDA.

add_library(magnets STATIC
    Magnets.cpp
    MagnetPhysics.cpp
)

target_include_directories(magnets
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(magnets PRIVATE -Wall -Wextra)