# src/Control/CMakeLists.txt
# Plasma control: PID loops, disruption prediction, SCRAM, runaway electrons.
# The ONLY module that writes setpoints to ReactorState.

add_library(control STATIC
    Control.cpp
    ControlPhysics.cpp
)

target_include_directories(control
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(control PRIVATE -Wall -Wextra)