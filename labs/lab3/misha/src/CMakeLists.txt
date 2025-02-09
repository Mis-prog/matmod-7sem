cmake_minimum_required(VERSION 3.15)

add_executable(lab3_misha_c main_sol.c)
target_compile_options(lab3_misha_c PRIVATE -O3)

add_executable(lab3_misha_cpp main_sol.cpp)
target_compile_options(lab3_misha_cpp PRIVATE -O3)