cmake_minimum_required(VERSION 3.15)


add_subdirectory(src)

add_executable(lab3_misha main.cpp)
target_include_directories(
        lab3_misha
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(lab3_misha lab3_misha_lib)

#target_compile_options(
#        lab3_misha
#        PRIVATE
#        -O3                      # Оптимизация для скорости
#        -funroll-loops            # Развёртывание циклов
#        -march=native             # Использовать инструкцию для текущей архитектуры процессора
#        -flto                     # Включение Link Time Optimization (LTO) для улучшения оптимизации на этапе линковки
#        -fno-rtti                 # Отключение RTTI (если не используется, можно ускорить выполнение)
#        -fomit-frame-pointer
#)