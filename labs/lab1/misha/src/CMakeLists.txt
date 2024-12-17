cmake_minimum_required(VERSION 3.15)


set(PATH_GSL_INCLUDE "C:/Cygwin/usr/include/gsl")
set(PATH_BOOST_INCLUDE "C:/Cygwin/usr/include/boost")

#set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_RELEASE} \
#    -Ofast \
#    -march=native \
#    -flto \
#    -funroll-loops \
#    -ftree-vectorize \
#    -ffast-math\
#"
#)

#set (CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set (CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#add_executable(lab1_task1 task1.cpp)
#add_executable(lab1_task2 task2.cpp)
#
#target_include_directories(lab1_task1 PRIVATE ${PATH_GSL_INCLUDE})
#target_include_directories(lab1_task2 PRIVATE ${PATH_GSL_INCLUDE})
#
#target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)
#target_link_libraries(lab1_task2 PRIVATE gsl gslcblas)

find_package(OpenMP REQUIRED)

#find_package(Boost REQUIRED)
add_executable(lab1_task1 task1.cpp)
target_include_directories(lab1_task1 PRIVATE ${PATH_BOOST_INCLUDE})

add_executable(lab1_task2 task2.cpp)
target_include_directories(lab1_task2 PRIVATE ${PATH_BOOST_INCLUDE})
target_link_libraries(lab1_task2 PRIVATE OpenMP::OpenMP_CXX)