#!/usr/bin/python3

import numpy as np

import ctypes
lib = ctypes.cdll.LoadLibrary('../build/libfit_lib.so')

class ParticlePOD(ctypes.Structure):
    # This implementation is due to recursive definition (tree).
    def hw(self):
        print("Hi!")

ParticlePOD._fields_ = [
    ('pdg', ctypes.c_int), 
    ('daughters_size', ctypes.c_int), 
    ('daughters', ctypes.POINTER(ctypes.POINTER(ParticlePOD))),
    ('momentum', ctypes.c_double * 3),
    ('covariance', (ctypes.c_double * 3) * 3)
]

def wrap_function(lib, funcname, restype, argtypes):
    """Simplify wrapping ctypes functions"""
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func

# Wrappers 
particle_from_daughters = wrap_function(lib, 'Particle_from_daughters', ctypes.POINTER(ParticlePOD), 
    [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(ParticlePOD))])

particle_from_momentum = wrap_function(lib, 'Particle_from_momentum', ctypes.POINTER(ParticlePOD), 
    [ctypes.c_int, ctypes.POINTER(ctypes.c_double), np.ctypeslib.ndpointer(dtype=np.float64,
            ndim=2,
            flags='C_CONTIGUOUS'
            )]
        )
    # ctypes.POINTER(ctypes.c_double)])

fit = wrap_function(lib, 'fit', ctypes.c_double, [ctypes.POINTER(ParticlePOD)])

# Python API for sct::Particle
class Particle:
    def __init__(self, pdg, momentum=None, cov=None, daughters=None):
        if daughters:
            daughters = [daughter.particle for daughter in daughters]
            array_type = ctypes.POINTER(ParticlePOD) * len(daughters)
            self.particle = particle_from_daughters(pdg, len(daughters), array_type(*daughters))
        elif momentum is not None and cov is not None:
            mom = (ctypes.c_double * len(momentum)) (*momentum)
            cov = np.array(cov, dtype=np.float64)
            self.particle = particle_from_momentum(pdg, mom, cov)
        else:
            raise RuntimeError("Particle created from momentum/covariance or daughters: define it.") 

    def fit(self):
        return fit(self.particle)

    def momentum(self):
        return np.array(self.particle.contents.momentum)


if __name__ == "__main__":
    a = Particle(211, momentum=[0.15730242, 0.0807638,  0.0967367], cov=[[9e-6, 0, 0], [0, 9e-6, 0], [0, 0, 9e-6]])
    b = Particle(-211, momentum=[-0.15941372, -0.08322166, -0.10307127], cov=[[9e-6, 0, 0], [0, 9e-6, 0], [0, 0, 9e-6]])
    c = Particle(311, daughters=[a, b])
    c.fit()    
