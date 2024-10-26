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

add_executable(lab1_task1_sun task1.cpp)

target_include_directories(lab1_task1_sun PRIVATE "C:/cygwin64/usr/include")

target_link_libraries(lab1_task1_sun PRIVATE gsl gslcblas)
