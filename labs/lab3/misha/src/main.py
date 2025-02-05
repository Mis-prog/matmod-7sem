from lab3 import lib
import numpy as np
import cffi
import matplotlib.pyplot as plt

ffi = cffi.FFI()

# Инициализация данных
N = 1000
count = 1000000
tau = 0.01
alpha = 0
beta = 100
q = np.zeros(N, dtype=np.float64)
v = np.zeros(N, dtype=np.float64)
a = np.zeros(N, dtype=np.float64)
q[N // 2] = 0.5
q[N // 2 - 1] = -0.5

q_ptr = ffi.cast("double*", ffi.from_buffer(q))
v_ptr = ffi.cast("double*", ffi.from_buffer(v))
a_ptr = ffi.cast("double*", ffi.from_buffer(a))


# lib.SimplexVerle(q_ptr, v_ptr, N, 0.0, 100.0, 0.01, 1000000, 1.0)

for i in range(count):
    lib.SimplexVerleNew(q_ptr, v_ptr, a_ptr, N, alpha, beta, tau, 1)

