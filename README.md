[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![codecov](https://codecov.io/gh/dbbs-lab/cerebellum/branch/master/graph/badge.svg)](https://codecov.io/gh/dbbs-lab/cerebellum)

# Cerebellum: A full scaffold model of the cerebellum, using the BSB.
This repository provides the code, configuration and morphology data to reconstruct and simulate 
cerebellar cortex circuits using the [Brain Scaffold Builder](https://github.com/dbbs-lab/bsb) 
framework. 

## Installation
The `cerebellum` package requires Python 3.9+.

Developers best use pip's *editable* install. This creates a live link between the
installed package and the local git repository:

```
 git clone git@github.com:dbbs-lab/cerebellum
 cd bsb
 pip install -e .[dev]
 pre-commit install
```

## Usage

### Reconstruction of the cerebellar cortex
Reconstruction 

# Acknowledgements
This research has received funding from the European Union’s Horizon 2020 Framework
Program for Research and Innovation under the Specific Grant Agreement No. 945539
(Human Brain Project SGA3) and Specific Grant Agreement No. 785907 (Human Brain
Project SGA2) and from Centro Fermi project “Local Neuronal Microcircuits” to ED. We
acknowledge the use of EBRAINS platform and Fenix Infrastructure resources, which are
partially funded from the European Union’s Horizon 2020 research and innovation
programme through the ICEI project under the grant agreement No. 800858
