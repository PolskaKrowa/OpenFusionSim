# src/PlasmaCore/CMakeLists.txt
# CUDA PIC simulation kernel — GPU-side only.
# Exposes a thin C++ host wrapper (PlasmaCoreBridge) to the rest of the game.

cmake_minimum_required(VERSION 3.20)

# ─── CUDA sources ─────────────────────────────────────────────────────────────
set(PLASMACORE_CUDA_SOURCES
    boris_push.cu
    field_interpolation.cu
    charge_deposit.cu
    field_solve.cu
    collisions.cu
    fusion_reactions.cu
    neutron_transport.cu
    sorting.cu
    plasmacore.cu
)

# ─── Host-side bridge (plain C++) ────────────────────────────────────────────
set(PLASMACORE_CXX_SOURCES
    PlasmaCoreBridge.cpp
)

add_library(plasmacore STATIC
    ${PLASMACORE_CUDA_SOURCES}
    ${PLASMACORE_CXX_SOURCES}
)

target_include_directories(plasmacore
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}    # shared include/
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} # types.cuh etc. local to PlasmaCore
)

target_link_libraries(plasmacore
    PUBLIC
        CUDA::cudart
        CUDA::cufft
        CUDA::curand
        CUDA::cusparse
        CUDA::cusolver
)

set_target_properties(plasmacore PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_compile_options(plasmacore PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --extra-device-vectorization
        -lineinfo
        --ptxas-options=-v
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -Wall -Wextra
    >
)