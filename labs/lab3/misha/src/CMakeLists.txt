cmake_minimum_required(VERSION 3.15)

file(GLOB SRCS *.cpp)

add_library(
        lab3_misha_lib
        STATIC
        ${SRCS}
)