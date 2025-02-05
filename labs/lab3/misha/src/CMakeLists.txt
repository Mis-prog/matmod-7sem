cmake_minimum_required(VERSION 3.15)

add_executable(lab3_misha main_sol.c)
target_compile_options(lab3_misha PRIVATE -O3)

add_library(lab3_misha_lib SHARED main_lib.c)
set_target_properties(lab3_misha_lib PROPERTIES COMPILE_FLAGS "-fPIC")
set_target_properties(lab3_misha_lib PROPERTIES OUTPUT_NAME "lab3")
