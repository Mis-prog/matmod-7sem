cmake_minimum_required(VERSION 3.15)


set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
    -Ofast \
    -march=native \
    -flto \
    -funroll-loops \
    -ftree-vectorize \
    -ffast-math \
    -pipe \
    -DNDEBUG"
)

find_package(OpenMP REQUIRED)

add_executable(lab1_task1 task1.cpp)

target_include_directories(lab1_task1 PRIVATE "C:/cygwin64/usr/include/gsl")


target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)

