import os
import cffi

ffi = cffi.FFI()

# Путь к вашему заголовочному файлу
PATH = os.getcwd()

# Открываем и читаем заголовочный файл
with open(os.path.join(PATH, "main_lib.h")) as f:
    ffi.cdef(f.read())

# Устанавливаем путь к файлам заголовков
ffi.set_source("lab3",
               '''
               #include "main_lib.h"  # Убедитесь, что путь к заголовочному файлу правильный
               ''',
               sources=["main_lib.c"]
               )

# Компиляция
ffi.compile(tmpdir='.')  # Генерация C расширения для Python
