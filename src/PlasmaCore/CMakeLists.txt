# src/PlasmaCore/CMakeLists.txt
# CUDA PIC simulation kernel.
# constants.cu is the single definition point for all shared __constant__ vars.

cmake_minimum_required(VERSION 3.20)

set(PLASMACORE_CUDA_SOURCES
    constants.cu              # <-- MUST be first: defines all __constant__ symbols
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

set(PLASMACORE_CXX_SOURCES
    PlasmaCoreBridge.cpp
)

add_library(plasmacore STATIC
    ${PLASMACORE_CUDA_SOURCES}
    ${PLASMACORE_CXX_SOURCES}
)

target_include_directories(plasmacore
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(plasmacore
    PUBLIC
        CUDA::cudart
        CUDA::cufft
        CUDA::curand
        CUDA::cusparse
        CUDA::cusolver
)

# Separable compilation is required so that __constant__ symbols defined in
# constants.cu are visible to kernels in other .cu files at device-link time.
# Without CUDA_SEPARABLE_COMPILATION=ON, each .cu is compiled as a standalone
# device binary and the extern __constant__ references would go unresolved.
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