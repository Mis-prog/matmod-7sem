cmake_minimum_required(VERSION 3.15)

# Укажите исполняемый файл
add_executable(lab1_task1 task1.cpp)

target_include_directories(lab1_task1 PRIVATE "C:/Cygwin/usr/include")

# Явно добавьте пути к библиотекам
target_link_directories(lab1_task1 PRIVATE "C:/Cygwin/usr/lib")

# Явно укажите библиотеки, используя их имена
target_link_libraries(lab1_task1 PRIVATE gsl gslcblas)
