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
import wrapper as fitter

pip = fitter.Particle(211, momentum=[0.1, 0.2, 0.3])
pim = fitter.Particle(-211, momentum=[-0.1, -0.2, -0.3])
kaon = fitter.Particle(311, daughters=[pip, pim])
chi = kaon.fit()
```

#### C++
``` C++
todo
```