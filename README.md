[![Build Status](https://travis-ci.org/ddyak/sct-fitter.svg?branch=master)](https://travis-ci.org/ddyak/sct-fitter)

# TreeFitter

## About

Kinematic fitter based on Belle II paper ([arxiv](https://arxiv.org/pdf/1901.11198.pdf)). This project develops as stand-alone, but the results will be used in the Super C-Tau factory framework Aurora.

## Install

#### clone 
```
git clone --recursive https://github.com/ddyak/sct-fitter
```

or after clone

```
git submodule update --init
```

#### build

```
mkdir build
cd build
cmake ..
make -j
```

## Example

All examples in `examples` directory.

#### Python
```python
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
```

#### C++
``` C++
todo
```