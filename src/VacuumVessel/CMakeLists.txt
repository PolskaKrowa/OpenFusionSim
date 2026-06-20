# src/VacuumVessel/CMakeLists.txt
# Vacuum vessel pumping: roughing + turbo + bakeout + boronization.
# Interlocks plasma initiation on base vacuum.

add_library(vacuumvessel STATIC
    VacuumVessel.cpp
)

target_include_directories(vacuumvessel
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(vacuumvessel PRIVATE -Wall -Wextra)
