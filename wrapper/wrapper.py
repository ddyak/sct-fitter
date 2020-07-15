#!/usr/bin/python3

import ctypes
lib = ctypes.cdll.LoadLibrary('../build/libfit_lib.so')

def wrap_function(lib, funcname, restype, argtypes):
    """Simplify wrapping ctypes functions"""
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func

class Particle(ctypes.Structure):
    # This implementation is due to recursive definition (tree).
    pass

Particle._fields_ = [
    ('momentum', ctypes.c_double * 3), 
    ('pdg', ctypes.c_int), 
    ('daughters_size', ctypes.c_int), 
    ('daughters', ctypes.POINTER(ctypes.POINTER(Particle)))
]

def constr(pdg, momentum=None, daughters=None):
    if daughters:
        get_node = wrap_function(lib, 'Particle_from_daughters', ctypes.POINTER(Particle), [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(Particle))])
        array_type = ctypes.POINTER(Particle) * len(daughters)
        return get_node(pdg, len(daughters), array_type(*daughters))
    
    if momentum:
        get_node = wrap_function(lib, 'Particle_from_momentum', ctypes.POINTER(Particle), [ctypes.c_int, ctypes.POINTER(ctypes.c_double)])
        array_type = ctypes.c_double * len(momentum)
        return get_node(pdg, array_type(*momentum))

a = constr(1, momentum=[30, 30, 30])
b = constr(2, momentum=[100, 1, 20])
c = constr(3, daughters=[a, b])

fit = wrap_function(lib, 'fit', None, [ctypes.POINTER(Particle)])
fit(c)

print(a.contents.pdg)
print(b.contents.pdg)
print(c.contents.pdg)
