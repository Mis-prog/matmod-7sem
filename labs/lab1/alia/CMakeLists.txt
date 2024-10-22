cmake_minimum_required(VERSION 3.15)


add_executable(lab1_task1_sun task1.cpp)

target_include_directories(lab1_task1_sun PRIVATE "C:/cygwin64/usr/include")

target_link_directories(lab1_task1_sun PRIVATE "C:/cygwin64/usr/share")

target_link_libraries(lab1_task1_sun PRIVATE gsl gslcblas)
