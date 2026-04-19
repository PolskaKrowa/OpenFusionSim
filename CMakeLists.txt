cmake_minimum_required(VERSION 3.20)
project(FusionSim VERSION 0.1.0 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD  17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 80 86 89 90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ─── Dependencies ─────────────────────────────────────────────────────────────
find_package(CUDAToolkit REQUIRED)
find_package(SDL2        REQUIRED)
find_package(OpenGL      REQUIRED)

# ─── Dear ImGui (vendored in third_party/imgui) ───────────────────────────────
#
#  Setup (one-time):
#    git clone https://github.com/ocornut/imgui third_party/imgui
#
#  No CMakeLists in that directory — we compile the sources directly here.
#
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/third_party/imgui)

add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/backends/imgui_impl_sdl2.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${SDL2_INCLUDE_DIRS}
)
target_link_libraries(imgui PUBLIC ${SDL2_LIBRARIES} OpenGL::GL)

# ─── Shared headers ───────────────────────────────────────────────────────────
set(FUSIONSIM_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# ─── Physics modules ──────────────────────────────────────────────────────────
add_subdirectory(src/PlasmaCore)
add_subdirectory(src/Magnets)
add_subdirectory(src/Fuel)
add_subdirectory(src/Helium)
add_subdirectory(src/ThermalHydraulics)
add_subdirectory(src/SteamPower)
add_subdirectory(src/Control)
add_subdirectory(src/TurbineSystem)
add_subdirectory(src/ElectricalGrid)
add_subdirectory(src/MoltenSalt)
add_subdirectory(src/HeliumSystem)

# ─── Executable ───────────────────────────────────────────────────────────────
add_executable(fusionsim src/main.cpp)

target_include_directories(fusionsim PRIVATE
    ${FUSIONSIM_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src      # #include "Control/Control.h" etc.
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(fusionsim PRIVATE
    plasmacore thermalhydraulics steampower
    magnets fuel helium control
    turbinesystem electricalgrid moltensalt heliumsystem
    imgui ${SDL2_LIBRARIES} OpenGL::GL CUDA::cudart
)

target_compile_options(fusionsim PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>
)

# Copy assets/ (fonts etc.) next to the binary after build
add_custom_command(TARGET fusionsim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/assets $<TARGET_FILE_DIR:fusionsim>/assets
    COMMENT "Copying assets to binary directory")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

message(STATUS "FusionSim ${PROJECT_VERSION}  |  ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA arch: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "ImGui:     ${IMGUI_DIR}")
message(STATUS "SDL2:      ${SDL2_LIBRARIES}")