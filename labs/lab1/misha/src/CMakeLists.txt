cmake_minimum_required(VERSION 3.15)

set(PATH_BOOST_INCLUDE "C:/cygwin64/usr/include/boost")

find_package(OpenMP REQUIRED)

add_executable(lab1_task1 task1.cpp)
target_include_directories(lab1_task1 PRIVATE ${PATH_BOOST_INCLUDE})

add_executable(lab1_task2 task2.cpp)
target_include_directories(lab1_task2 PRIVATE ${PATH_BOOST_INCLUDE})
target_link_libraries(lab1_task2 PRIVATE OpenMP::OpenMP_CXX)