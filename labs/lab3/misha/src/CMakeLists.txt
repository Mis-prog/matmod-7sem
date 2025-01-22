cmake_minimum_required(VERSION 3.15) 

add_executable(lab3_misha main.cpp)
target_compile_options(lab3_misha 
            PRIVATE 
            -O3 
            -march=native 
            -funroll-loops 
            -flto 
            -fomit-frame-pointer 
            -pipe
        )