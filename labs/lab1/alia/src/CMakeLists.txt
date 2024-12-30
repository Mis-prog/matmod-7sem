cmake_minimum_required(VERSION 3.15)

set(PATH_BOOST_INCLUDE "C:/cygwin64/usr/include/boost")

add_executable(lab1_task1_sun task1.cpp)
target_include_directories(lab1_task1_sun PRIVATE ${PATH_BOOST_INCLUDE})