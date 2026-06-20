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

# ═══════════════════════════════════════════════════════════════════════════════
# BUILD TYPE SELECTION
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Pick a build type by passing -DCMAKE_BUILD_TYPE=<type> on the configure line.
#  Supported types:
#
#    Debug           — no optimisation, full debug info, asserts ON.
#                      Use this with gdb / cuda-gdb to step through code.
#                      PIC kernels will be slow but stack traces are clean.
#
#    Release         — -O3, no debug info, asserts OFF, CUDA --use_fast_math ON.
#                      Use this for benchmarking / final runs.
#
#    RelWithDebInfo  — -O2 -g, asserts ON, CUDA --use_fast_math ON.
#                      THE SWEET SPOT for development: fast + debuggable.
#                      Recommended default for working on the simulator.
#
#    FastDebug       — -O1 -g, asserts ON, no fast-math.  A custom type for
#                      when Debug is too slow but RelWithDebInfo's -O2 hides
#                      too much under gdb.  Useful for hunting memory bugs.
#
#  Default if user doesn't specify: RelWithDebInfo (matches typical user
#  expectation — optimised but still debuggable).
#
#  For multi-config generators (Visual Studio, Ninja Multi-Config), the build
#  type is chosen at BUILD time via --config <type> rather than at configure
#  time, so we set CMAKE_CONFIGURATION_TYPES instead.
#
if(NOT CMAKE_CONFIGURATION_TYPES)
    # Single-config generator (Make, Ninja)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
            "Build type (Debug|Release|RelWithDebInfo|FastDebug)" FORCE)
    endif()
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo FastDebug)
else()
    # Multi-config generator (VS, Ninja Multi-Config)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;FastDebug"
        CACHE STRING "Available configurations" FORCE)
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# OPTIONAL FEATURE FLAGS
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Pass -DFUSIONSIM_ENABLE_ASAN=ON to enable the AddressSanitizer.
#  Catches buffer overruns, use-after-free, stack-overflow, etc. at the cost
#  of ~2× runtime slowdown and ~3× memory.  Incompatible with CUDA device
#  code (host-only), but invaluable for catching the kind of memory bugs
#  that have been showing up in the bridge code.
#
option(FUSIONSIM_ENABLE_ASAN "Enable AddressSanitizer (host code only)" OFF)
#
#  Pass -DFUSIONSIM_CUDA_DEBUG=ON to compile CUDA kernels with full debug info
#  and NO optimisations.  Required for cuda-gdb stepping into device code.
#  Auto-enabled in the Debug build type; can be forced on for other types.
#
option(FUSIONSIM_CUDA_DEBUG "Compile CUDA with -G (cuda-gdb friendly)" OFF)
#
#  Pass -DFUSIONSIM_WARNINGS_AS_ERRORS=ON to treat compiler warnings as errors.
#  Useful for CI; too noisy for interactive development.
#
option(FUSIONSIM_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

# ═══════════════════════════════════════════════════════════════════════════════
# PER-BUILD-TYPE FLAGS
# ═══════════════════════════════════════════════════════════════════════════════
#
#  CMake's built-in defaults (CMAKE_CXX_FLAGS_DEBUG etc.) are reasonable but
#  not aggressive enough.  We override them with project-appropriate values.
#
#  Note: CUDA flags are separate from C++ flags — nvcc has its own dialect
#  (--use_fast_math, --opt-level, -G, etc.) that doesn't map cleanly to g++
#  flags.  We use generator expressions to apply them per-language.
#

# ─── Debug ────────────────────────────────────────────────────────────────────
#  No optimisation, full debug info, all warnings, asserts ON (NDEBUG not set).
set(CMAKE_CXX_FLAGS_DEBUG
    "-O0 -g3 -ggdb3 -Wall -Wextra -Wpedantic -D_GLIBCXX_DEBUG -D_GLIBCXX_ASSERTIONS"
    CACHE STRING "C++ Debug flags" FORCE)
set(CMAKE_CUDA_FLAGS_DEBUG
    "-O0 -g -G --ptxas-options=-v -lineinfo -Xcompiler -Wall -Xcompiler -Wextra"
    CACHE STRING "CUDA Debug flags" FORCE)
# Note: -G is nvcc's device-debug flag (disables all device-side optimisations
# and emits full line tables).  Essential for cuda-gdb but makes kernels ~10×
# slower.  Don't use for performance work.

# ─── Release ──────────────────────────────────────────────────────────────────
#  Full optimisation, NDEBUG set (asserts compiled out), CUDA fast-math ON.
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG -fno-stack-protector"
    CACHE STRING "C++ Release flags" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE
    "-O3 --use_fast_math --extra-device-vectorization -DNDEBUG"
    CACHE STRING "CUDA Release flags" FORCE)

# ─── RelWithDebInfo ───────────────────────────────────────────────────────────
#  Optimised + debuggable.  -O2 (not -O3) so the debugger view isn't too
#  scrambled by inlining.  CUDA gets --use_fast_math but NOT -G (device debug
#  info is still emitted via -lineinfo, which is much cheaper than -G).
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "-O2 -g3 -ggdb3 -DNDEBUG -D_GLIBCXX_ASSERTIONS"
    CACHE STRING "C++ RelWithDebInfo flags" FORCE)
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO
    "-O2 -lineinfo --use_fast_math -DNDEBUG"
    CACHE STRING "CUDA RelWithDebInfo flags" FORCE)

# ─── FastDebug (custom) ──────────────────────────────────────────────────────
#  -O1 keeps the debugger view mostly intact while still applying basic
#  optimisations (instruction combining, etc.).  Asserts ON.  CUDA gets -lineinfo
#  but no -G and no --use_fast_math (so numerical behaviour matches Debug).
set(CMAKE_CXX_FLAGS_FASTDEBUG
    "-O1 -g3 -ggdb3 -Wall -Wextra -D_GLIBCXX_ASSERTIONS"
    CACHE STRING "C++ FastDebug flags" FORCE)
set(CMAKE_CUDA_FLAGS_FASTDEBUG
    "-O1 -lineinfo -Xcompiler -Wall -Xcompiler -Wextra"
    CACHE STRING "CUDA FastDebug flags" FORCE)
# Tell CMake this is a valid configuration (it doesn't know about FastDebug by default)
if(NOT CMAKE_CONFIGURATION_TYPES)
    list(APPEND CMAKE_CXX_FLAGS_FASTDEBUG_VALID Build FastDebug)
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# ASAN / CUDA-DEBUG / WARNINGS-AS-ERRORS HANDLING
# ═══════════════════════════════════════════════════════════════════════════════
#
#  These are applied as INTERFACE options on a custom target so they propagate
#  to every library and executable in the project without having to touch each
#  subdirectory's CMakeLists.txt.  (Pattern: an "interface library" is CMake's
#  recommended way to share compiler flags across a project.)
#
add_library(fusionsim_compile_options INTERFACE)

if(FUSIONSIM_ENABLE_ASAN)
    message(STATUS "AddressSanitizer: ENABLED (host code only)")
    target_compile_options(fusionsim_compile_options INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize-recover=all>
    )
    target_link_options(fusionsim_compile_options INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address -fsanitize=undefined>
    )
    # ASAN + CUDA: nvcc doesn't understand -fsanitize, so we only apply it to
    # host C++.  Device code is unaffected.  Note that ASAN will not catch
    # device-side memory bugs (use cuda-memcheck for those).
endif()

if(FUSIONSIM_CUDA_DEBUG)
    message(STATUS "CUDA device debug: ENABLED (-G, no fast-math)")
    target_compile_options(fusionsim_compile_options INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:-G -lineinfo>
    )
    # Disable fast-math if it was enabled by the build type
    target_compile_options(fusionsim_compile_options INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:--no-use_fast_math>
    )
endif()

if(FUSIONSIM_WARNINGS_AS_ERRORS)
    message(STATUS "Warnings as errors: ENABLED")
    target_compile_options(fusionsim_compile_options INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Werror>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -Werror>
    )
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# DEPENDENCIES
# ═══════════════════════════════════════════════════════════════════════════════
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
target_link_libraries(imgui PUBLIC ${SDL2_LIBRARIES} OpenGL::GL fusionsim_compile_options)

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
# ── New modules (Round 4) ─────────────────────────────────────────────────────
#  HCD                  — NBI/ICRH/ECRH/LHCD actuators (replaces hardcoded
#                        P_aux = 50 MW in PlasmaCoreBridge)
#  VacuumVessel         — roughing + turbo pumps, bakeout, boronization;
#                        interlocks plasma initiation on base vacuum
#  DisruptionMitigation — MGI/SPI systems for operator-triggered disruption
#                        mitigation
#  TritiumPlant         — TES extraction, isotope separation, detritiation
add_subdirectory(src/HCD)
add_subdirectory(src/VacuumVessel)
add_subdirectory(src/DisruptionMitigation)
add_subdirectory(src/TritiumPlant)

# ─── Executable ───────────────────────────────────────────────────────────────
add_executable(fusionsim
    src/main.cpp
    src/PlasmaVizTabs.cpp      # 3 plasma visualization tabs (temp/fusion/particles)
)

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
    hcd vacuumvessel disruptionmitigation tritiumplant
    imgui ${SDL2_LIBRARIES} OpenGL::GL CUDA::cudart
    fusionsim_compile_options     # propagates ASAN / CUDA_DEBUG / Werror flags
)

target_compile_options(fusionsim PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>
)

# Copy assets/ (fonts etc.) next to the binary after build
add_custom_command(TARGET fusionsim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/assets $<TARGET_FILE_DIR:fusionsim>/assets
    COMMENT "Copying assets to binary directory")

# ═══════════════════════════════════════════════════════════════════════════════
# BUILD SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
if(CMAKE_CONFIGURATION_TYPES)
    message(STATUS "FusionSim ${PROJECT_VERSION}  |  multi-config: ${CMAKE_CONFIGURATION_TYPES}")
else()
    message(STATUS "FusionSim ${PROJECT_VERSION}  |  build type: ${CMAKE_BUILD_TYPE}")
endif()
message(STATUS "CUDA arch:    ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "ASAN:         ${FUSIONSIM_ENABLE_ASAN}")
message(STATUS "CUDA debug:   ${FUSIONSIM_CUDA_DEBUG}")
message(STATUS "Werror:       ${FUSIONSIM_WARNINGS_AS_ERRORS}")
message(STATUS "ImGui:        ${IMGUI_DIR}")
message(STATUS "SDL2:         ${SDL2_LIBRARIES}")

# ═══════════════════════════════════════════════════════════════════════════════
# USAGE
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Standard debug build:
#    cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug
#    cmake --build build-debug -j
#
#  Optimised release build (for benchmarking):
#    cmake -B build-release -DCMAKE_BUILD_TYPE=Release
#    cmake --build build-release -j
#
#  RelWithDebInfo (default — recommended for development):
#    cmake -B build
#    cmake --build build -j
#
#  ASAN build (catch host memory bugs):
#    cmake -B build-asan -DCMAKE_BUILD_TYPE=FastDebug -DFUSIONSIM_ENABLE_ASAN=ON
#    cmake --build build-asan -j
#    ASAN_OPTIONS=detect_leaks=1:halt_on_error=0 ./build-asan/bin/fusionsim
#
#  cuda-gdb build (step into device code):
#    cmake -B build-cudbg -DCMAKE_BUILD_TYPE=Debug -DFUSIONSIM_CUDA_DEBUG=ON
#    cmake --build build-cudbg -j
#    cuda-gdb ./build-cudbg/bin/fusionsim
#
#  Warnings-as-errors (for CI):
#    cmake -B build-ci -DCMAKE_BUILD_TYPE=RelWithDebInfo -DFUSIONSIM_WARNINGS_AS_ERRORS=ON
#
#  Switch build type without reconfiguring from scratch:
#    cmake -B build -DCMAKE_BUILD_TYPE=Release
#    cmake --build build -j   # only the changed files recompile
#
