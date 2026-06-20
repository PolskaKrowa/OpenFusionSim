# src/HCD/CMakeLists.txt
# Heating and Current Drive actuators: NBI, ICRH, ECRH, LHCD.
# Replaces the hardcoded `P_aux = 50 MW` in PlasmaCoreBridge.

add_library(hcd STATIC
    HCD.cpp
)

target_include_directories(hcd
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(hcd PRIVATE -Wall -Wextra)
