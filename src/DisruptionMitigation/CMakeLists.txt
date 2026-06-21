# src/DisruptionMitigation/CMakeLists.txt
# Disruption mitigation: MGI (Massive Gas Injection) + SPI (Shattered Pellet
# Injection).  Operator-triggered, reduces halo current forces.

add_library(disruptionmitigation STATIC
    DisruptionMitigation.cpp
)

target_include_directories(disruptionmitigation
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(disruptionmitigation PRIVATE -Wall -Wextra)
