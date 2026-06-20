# src/TritiumPlant/CMakeLists.txt
# Tritium plant: TES extraction, isotope separation, detritiation.

add_library(tritiumplant STATIC
    TritiumPlant.cpp
)

target_include_directories(tritiumplant
    PUBLIC  ${FUSIONSIM_INCLUDE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(tritiumplant PRIVATE -Wall -Wextra)
