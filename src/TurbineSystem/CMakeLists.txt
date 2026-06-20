# src/TurbineSystem/CMakeLists.txt
add_library(turbinesystem STATIC
    TurbineUnit.cpp
    TurbineSystem.cpp
)
target_include_directories(turbinesystem
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(turbinesystem PRIVATE -Wall -Wextra)