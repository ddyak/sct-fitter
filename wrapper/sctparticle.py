from wrapper import *
import numpy as np


class Particle:
    def __init__(self, pdg, momentum=None, cov=None, daughters=None, track=None, cluster=None):
        if momentum is not None and cov is not None:
            momentum = np.array(momentum, dtype=np.float64)
            cov = np.array(cov, dtype=np.float64)
            self.particle = particle_from_momentum(pdg, momentum, cov)
        elif daughters is not None:
            daughters = [daughter.particle for daughter in daughters]
            daughters = np.array(daughters)
            self.particle = particle_from_daughters(pdg, len(daughters), daughters.ctypes.data_as(ctypes.POINTER(ctypes.c_void_p)))
        elif track is not None:
            pass
        elif cluster is not None:
            pass
        else:
            raise Exception("Invalid arguments.")

    def momentum(self):
        return get_momentum(self.particle)
    
    def momentum3(self):
        pass

    def energy(self):
        pass

if __name__ == "__main__":
    a = Particle(211, momentum=np.random.rand(3), cov=np.random.rand(3, 3))
    b = Particle(-211, momentum=np.random.rand(3), cov=np.random.rand(3, 3))
    c = Particle(311, daughters=[a, b])
    for part in [a, b, c]:
        print(part.momentum())