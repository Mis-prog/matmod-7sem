cmake_minimum_required(VERSION 3.15)


set(PATH_GSL_INCLUDE "C:/Cygwin/usr/include/gsl")

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

target_include_directories(lab1_task1 PRIVATE ${PATH_GSL_INCLUDE})
target_include_directories(lab1_task2 PRIVATE ${PATH_GSL_INCLUDE})

target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)
target_link_libraries(lab1_task2 PRIVATE gsl gslcblas)


find_package(OpenMP)

set(SUNDIALS_INCLUDE_DIR "/usr/include" CACHE PATH "Location of SUNDIALS header files")
set(SUNDIALS_LIBRARY_DIR "/usr/lib" CACHE PATH "Location of SUNDIALS libraries")

find_library(SUNDIALS_SOLVER_LIB
        sundials_cvode
        PATHS ${SUNDIALS_LIBRARY_DIR}
        DOC "CVODE library")
find_library(SUNDIALS_NVEC_LIB
        sundials_nvecserial
        PATHS ${SUNDIALS_LIBRARY_DIR}
        DOC "NVECTOR_SERIAL library")

if (NOT SUNDIALS_SOLVER_LIB OR NOT SUNDIALS_NVEC_LIB)
    message(FATAL_ERROR "SUNDIALS libraries not found! Check SUNDIALS_LIBRARY_DIR.")
endif()

add_executable(lab1_sundials task1-sundials.cpp)

target_include_directories(lab1_sundials PRIVATE ${SUNDIALS_INCLUDE_DIR})

target_link_libraries(lab1_sundials PRIVATE ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB})
target_link_libraries(lab1_sundials PUBLIC OpenMP::OpenMP_CXX)

