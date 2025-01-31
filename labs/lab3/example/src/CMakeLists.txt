cmake_minimum_required(VERSION 3.15) 

add_executable(lab3_my_example main_my.cpp)

target_compile_options(lab3_my_example
            PRIVATE
            -O3
            -march=native
            -funroll-loops
            -flto
            -fomit-frame-pointer
            -pipe
        )

add_executable(lab3_ernest_example main_er.cpp)

add_executable(lab3_angelina_example main_angelina.cpp)

add_executable(lab3_gumerov_example main_gumerov.cpp)