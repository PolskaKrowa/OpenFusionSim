# src/MoltenSalt/CMakeLists.txt
add_library(moltensalt STATIC MoltenSaltSystem.cpp)
target_include_directories(moltensalt
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_options(moltensalt PRIVATE -Wall -Wextra)