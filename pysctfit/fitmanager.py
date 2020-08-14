from sctparticle import Particle
from wrapper import *

class FitManager:
    def __init__(self, particle: Particle):
        self.fitmanager = create_manager(particle.particle)

    def fit(self):
        manager_fit(self.fitmanager)

    def momentum(self, part: Particle):
        return manager_momentum(self.fitmanager, part.particle)

    def cov_momentum(self, part: Particle):
        pass

    def momentum3(self, part: Particle):
        pass

    def cov_momentum3(self, part: Particle):
        pass

    def chi2(self):
        return manager_chi(self.fitmanager)


if __name__ == "__main__":
    a = Particle(211, momentum=np.random.rand(3), cov=np.random.rand(3, 3))
    b = Particle(-211, momentum=np.random.rand(3), cov=np.random.rand(3, 3))
    c = Particle(311, daughters=[a, b])

    fm = FitManager(c)
    fm.fit()
    print(fm.momentum(c))