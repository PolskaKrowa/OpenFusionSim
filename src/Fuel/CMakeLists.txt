# src/Fuel/CMakeLists.txt
# D-T fuel management: gas puffing, pellet injection, inventory tracking.

add_library(fuel STATIC
    Fuel.cpp
)

target_include_directories(fuel
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(fuel PRIVATE -Wall -Wextra)