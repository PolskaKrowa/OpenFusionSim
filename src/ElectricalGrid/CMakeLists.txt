# src/ElectricalGrid/CMakeLists.txt
add_library(electricalgrid STATIC ElectricalGrid.cpp)
target_include_directories(electricalgrid
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(electricalgrid PRIVATE -Wall -Wextra)