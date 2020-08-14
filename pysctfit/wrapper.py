import ctypes
lib = ctypes.cdll.LoadLibrary('../build/libfit_lib.so')

import numpy as np

def wrap_function(lib, funcname, restype, argtypes):
    """Simplify wrapping ctypes functions"""
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func

# Wrappers 

particle_from_momentum = wrap_function(lib, 'particle_from_momentum', ctypes.c_void_p, 
    [ctypes.c_int, np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')])

particle_from_daughters = wrap_function(lib, 'particle_from_daughters', ctypes.c_void_p, 
    [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p)])

get_momentum = wrap_function(lib, 'get_momentum', np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(4,)),
    [ctypes.c_void_p])

create_manager = wrap_function(lib, 'create_manager', ctypes.c_void_p, [ctypes.c_void_p])
    
manager_fit = wrap_function(lib, 'manager_fit', None, [ctypes.c_void_p])

manager_momentum = wrap_function(lib, 'manager_momentum', np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(4,)),
    [ctypes.c_void_p, ctypes.c_void_p])

manager_chi = wrap_function(lib, 'manager_chi', ctypes.c_double, [ctypes.c_void_p])
