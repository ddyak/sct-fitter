import numpy as np
import sys

sys.path.append("../pysctfit")

from sctparticle import Particle
from fitmanager  import FitManager

pip = Particle(211, momentum=np.random.rand(3), cov=np.diag([1e-3, 1e-3, 1e-3]))
pim = Particle(-211, momentum=np.random.rand(3), cov=np.diag([1e-3, 1e-3, 1e-3]))
kaon = Particle(311, daughters=[pip, pim])

fm = FitManager(kaon)
fm.fit()
chi2 = fm.chi2()