cmake_minimum_required(VERSION 3.15)

set(PATH_BOOST_INCLUDE "C:/cygwin64/usr/include/boost")
find_package(OpenMP REQUIRED)

add_executable(lab1_task1_sun task1.cpp)
target_include_directories(lab1_task1_sun PRIVATE ${PATH_BOOST_INCLUDE})

add_executable(lab1_task2_sun task2.cpp)
target_include_directories(lab1_task2_sun PRIVATE ${PATH_BOOST_INCLUDE})
target_link_libraries(lab1_task2_sun PRIVATE OpenMP::OpenMP_CXX)