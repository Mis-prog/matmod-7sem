cmake_minimum_required(VERSION 3.15)


set(PATH_LIB "C:/cygwin64/usr/include/gsl")
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

add_executable(lab1_task1 task1.cpp)
add_executable(lab1_task2 task2.cpp)

target_include_directories(lab1_task1 PRIVATE ${PATH_LIB})
target_include_directories(lab1_task2 PRIVATE ${PATH_LIB})

target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)
target_link_libraries(lab1_task2 PRIVATE gsl gslcblas)

