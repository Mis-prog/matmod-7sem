from cffi import FFI

ffi = FFI()

lib = ffi.dlopen(r"cyglab3_misha_lib.dll")

ffi.cdef("""
   double F(double q_last, double q, double q_next, double m, double alpha, double beta);
""")

