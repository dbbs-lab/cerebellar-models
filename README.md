[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![codecov](https://codecov.io/gh/dbbs-lab/cerebellum/branch/master/graph/badge.svg)](https://codecov.io/gh/dbbs-lab/cerebellum)

# Cerebellum: A full scaffold model of the cerebellum, using the BSB.
This repository provides the code, configuration and morphology data to reconstruct and simulate 
cerebellar cortex circuits using the [Brain Scaffold Builder](https://github.com/dbbs-lab/bsb) 
framework. These models are based on the iterative work of the [Department of Brain and Behavioral 
Sciences](https://dbbs.dip.unipv.it/en) of the university of Pavia. 

## Installation
The `cerebellum` package requires Python 3.9+.

Developers best use pip's *editable* install. This creates a live link between the
installed package and the local git repository:

```bash
 git clone git@github.com:dbbs-lab/cerebellum
 cd cerebellum
 pip install -e .[dev]
 pre-commit install
```

## Content

### Code
`cerebellum` is a package that contains implementation of connectivity or placement rules for BSB. 

You can produce the documentation for this code using the following commands, assuming you are in 
the `cerebellum` folder:
```bash
pip install -e .[docs]
cd docs
sphinx-build -n -b html . _build/html
```

You can also run the unit tests and functional tests using the following commands, assuming you are 
in the `cerebellum` folder:
```bash
pip install -e .[tests]
pytest test/
```

### Morphologies
Cerebellar cortex neuron morphology reconstructions used in our microcircuits are stored in the 
[morphologies](https://github.com/dbbs-lab/cerebellum/configurations) folder. The folder also 
contains a README with more information about each file.

### BSB configuraion files for the cerebellar cortex
In this repository, the BSB configurations are stored in the 
[configurations](https://github.com/dbbs-lab/cerebellum/configurations) folder. Sub-folders within 
`configurations` corresponds to different species reconstructions. Each specie should have its 
default configuration to reconstruct the models as well as `extensions` that can be added to the 
configuration. This includes the different simulation' paradigms.

## Building a circuit

Depending on the circuit you wish to obtain and/or simulate, the process to obtain the results will 
vary.

As an example, we present the process to reconstruct the cerebellar cortex of the rodent brain, 
based on the reconstruction of [De Schepper et al. 2022](https://doi.org/10.1038/s42003-022-04213-y) 
with BSB.

Assuming you are in the `cerebellum` folder, run the following command in your terminal:
```bash
bsb compile configurations/mouse/mouse_cerebellar_cortex.yaml -v4 --clear
```
This command will produce a microcircuit of the mouse cerebellar cortex and store it in the
`cerebellum.hdf5` file. This process might take a while depending on your machine.

# Acknowledgements
This research has received funding from the European Union’s Horizon 2020 Framework
Program for Research and Innovation under the Specific Grant Agreement No. 945539
(Human Brain Project SGA3) and Specific Grant Agreement No. 785907 (Human Brain
Project SGA2) and from Centro Fermi project “Local Neuronal Microcircuits” to ED. We
acknowledge the use of EBRAINS platform and Fenix Infrastructure resources, which are
partially funded from the European Union’s Horizon 2020 research and innovation
programme through the ICEI project under the grant agreement No. 800858
