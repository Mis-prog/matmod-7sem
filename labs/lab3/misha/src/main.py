from cffi import FFI
import numpy as np

ffi = FFI()

lib = ffi.dlopen(r"cyglab3_misha_lib.dll")

ffi.cdef("""
   void SimplexVerle(double *q, double *v, int size, double alpha, double beta, double tau, int N, double m, double *result_v);
""")

lib.SimplexVerle()

alpha = 0
beta = 100
tau = 0.01
N = 1000

q = np.zeros(N, dtype=np.float64)
v = np.zeros(N, dtype=np.float64)
result_v = np.zeros(N, dtype=np.float64)

q_c = ffi.cast("double*", q.ctypes.data)
v_c = ffi.cast("double*", v.ctypes.data)
result_v_c = ffi.cast("double*", result_v.ctypes.data)

print(q_c)
# lib.SimplexVerle(q_c, v_c, N, alpha, beta, tau, 1000000, 1, result_v_c)

# print(result_v[:10])
