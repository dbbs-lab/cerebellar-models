[![Build and Test](https://github.com/dbbs-lab/cerebellar-models/actions/workflows/integration.yml/badge.svg)](https://github.com/dbbs-lab/cerebellar-models/actions/workflows/integration.yml)
[![Documentation Status](https://readthedocs.org/projects/cerebellar-models/badge/?version=latest)](https://bsb.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![codecov](https://codecov.io/gh/dbbs-lab/cerebellar-models/graph/badge.svg?token=KBEE3D83YQ)](https://codecov.io/gh/dbbs-lab/cerebellar-models)

# DBBS Cerebellar models: Models of the cerebellum, using the BSB.
This repository provides the code, configuration and morphology data to reconstruct and simulate 
cerebellar cortex circuits using the [Brain Scaffold Builder](https://github.com/dbbs-lab/bsb) 
framework. These models are based on the iterative work of the [Department of Brain and Behavioral 
Sciences](https://dbbs.dip.unipv.it/en) of the university of Pavia. 

<h3>:closed_book: Read the documentation on  https://cerebellar-models.readthedocs.io/en/latest

## Installation
`cerebellar-models` is a package that contains implementation of connectivity or placement rules for BSB.
The `cerebellar-models` package requires Python 3.10+.

### Pre-requirements

BSB parallelizes the network reconstruction using MPI, and translates simulator instructions to 
the simulator backends with it as well, for effortless parallel simulation. 
To use MPI from Python the mpi4py package is required, which in turn needs a working MPI 
implementation installed in your environment.

On your local machine you can install OpenMPI:
```bash
sudo apt-get update && sudo apt-get install -y libopenmpi-dev openmpi-bin
```

Then, depending on the types of simulations, you want to perform you will need to install the 
simulator(s) of your interest. Please follow the instructions:
- For the [NEST](https://nest-simulator.readthedocs.io/en/stable/installation/index.html) simulator
- For the [NEURON](https://nrn.readthedocs.io/en/8.2.4/install/install.html) simulator

### Cerebellar-models installation

#### pip
Cerebellar-models can be installed from PyPI through `pip`:
```bash
pip install cerebellar-models
```

#### Developers
Developers best use pip's *editable* install. This creates a live link between the
installed package and the local git repository:

```bash
 git clone git@github.com:dbbs-lab/cerebellar-models
 cd cerebellar_models
 pip install -e .
```

## Contents

### Morphologies
Cerebellar cortex neuron morphology reconstructions used in our microcircuits are stored in the 
[morphologies](https://github.com/dbbs-lab/cerebellar-models/tree/master/configurations) folder. 
The folder contains also more information about each file.

### BSB configuration files for cerebellar cortex circuits
In this repository, the BSB configurations are stored in the 
[configurations](https://github.com/dbbs-lab/cerebellar-models/tree/master/configurations) folder. 
Sub-folders within `configurations` corresponds to different species reconstructions. Each specie 
have its default configuration to reconstruct the models as well as `extensions` that can be added 
to override or extend the default one. This includes the different simulation' paradigms.

## Building a circuit

Depending on the circuit you wish to obtain and/or simulate, the process will vary. Check the 
documentation for more information. 

As an example, we present the process to reconstruct the cerebellar cortex of the rodent brain, 
based on the reconstruction of [De Schepper et al. 2022](https://doi.org/10.1038/s42003-022-04213-y) 
with BSB.

Assuming you are in the `cerebellar-models` folder, run the following command in your terminal:
```bash
bsb compile configurations/mouse/mouse_cerebellar_cortex.yaml -v4 --clear
```
This command will produce a microcircuit of the mouse cerebellar cortex and store it in the
`cerebellum.hdf5` file. This process might take a while depending on your machine.

## Running a simulation
As for the previous paragraph the following command might vary depending on your reconstruction and 
simulation. 

Assuming you are in the `cerebellar-models` folder, and you want to run the simulation `simulation_name`, 
run the following command in your terminal:
```bash
bsb simulate cerebellum.hdf5 simulation_name -o output_file_name -v4
```

## Acknowledgements
This research has received funding from the European Union’s Horizon 2020 Framework
Program for Research and Innovation under the Specific Grant Agreement No. 945539
(Human Brain Project SGA3) and Specific Grant Agreement No. 785907 (Human Brain
Project SGA2) and from Centro Fermi project “Local Neuronal Microcircuits” to ED.
The project is also receiving funding from the Virtual Brain Twin Project under the 
European Union's Research and Innovation Program Horizon Europe under grant agreement
No 101137289.

We acknowledge the use of EBRAINS platform and Fenix Infrastructure resources, which are
partially funded from the European Union’s Horizon 2020 research and innovation
programme under the Specific Grant Agreement No. 101147319 (EBRAINS 2.0 Project) and 
through the ICEI project under the grant agreement No. 800858 respectively.
