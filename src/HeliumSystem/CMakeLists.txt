# src/HeliumSystem/CMakeLists.txt
add_library(heliumsystem STATIC HeliumCoolingSystem.cpp)
target_include_directories(heliumsystem
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(heliumsystem PRIVATE -Wall -Wextra)