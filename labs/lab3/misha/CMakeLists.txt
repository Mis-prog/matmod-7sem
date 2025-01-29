cmake_minimum_required(VERSION 3.15)


add_subdirectory(src)

add_executable(lab3_misha main.cpp)
target_include_directories(
        lab3_misha
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(lab3_misha lab3_misha_lib)

target_compile_options(
        lab3_misha
        PRIVATE
        -fno-fast-math
        -fno-finite-math-only
)