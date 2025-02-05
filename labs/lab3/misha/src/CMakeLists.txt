cmake_minimum_required(VERSION 3.15)

add_executable(lab3_misha main_sol.c)
target_compile_options(lab3_misha PRIVATE -O3)
add_library(lab3_misha_lib SHARED main_lib.c)
target_compile_options(lab3_misha_lib PRIVATE -O3)
