cmake_minimum_required(VERSION 3.20)
project(FusionSim VERSION 0.1.0 LANGUAGES CXX CUDA)

# ─── Global standards ──────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD  17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ─── CUDA architecture targets ────────────────────────────────────────────────
# Volta (V100), Ampere (A100/RTX30), Ada (RTX40), Hopper (H100)
set(CMAKE_CUDA_ARCHITECTURES 70 80 86 89 90)

# ─── Find packages ────────────────────────────────────────────────────────────
find_package(CUDAToolkit REQUIRED)

# Optional: SDL2 / OpenGL for the game frontend
# find_package(SDL2 QUIET)
# find_package(OpenGL QUIET)

# ─── Shared include directory ─────────────────────────────────────────────────
# All modules can see include/ without any extra target_include_directories.
set(FUSIONSIM_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# ─── Modules ──────────────────────────────────────────────────────────────────
add_subdirectory(src/PlasmaCore)
add_subdirectory(src/Magnets)
add_subdirectory(src/Fuel)
add_subdirectory(src/Helium)
add_subdirectory(src/ThermalHydraulics)
add_subdirectory(src/SteamPower)
add_subdirectory(src/Control)

# ─── Main executable ─────────────────────────────────────────────────────────
add_executable(fusionsim src/main.cpp)

target_include_directories(fusionsim PRIVATE ${FUSIONSIM_INCLUDE_DIR})

target_link_libraries(fusionsim PRIVATE
    plasmacore
    magnets
    fuel
    helium
    thermalhydraulics
    steampower
    control
)

# ─── Compiler warnings ────────────────────────────────────────────────────────
target_compile_options(fusionsim PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>
)

# ─── Build type defaults ──────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "FusionSim ${PROJECT_VERSION} — build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")