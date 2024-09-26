cmake_minimum_required(VERSION 3.15)


add_executable(lab1_task1 task1.cpp)
add_executable(lab1_task2 task2.cpp)

target_include_directories(lab1_task1 PRIVATE "C:/Cygwin/usr/include")

target_link_directories(lab1_task1 PRIVATE "C:/Cygwin/usr/lib")

target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)
