Building a first circuit
------------------------

In this section, you will learn how to reproduce the different cerebellar cortex circuits built at
the `DBBS <https://dbbs.dip.unipv.it/en)>`_ laboratory of the University of Pavia.

Introduction to BSB
~~~~~~~~~~~~~~~~~~~

The :doc:`Brain Scaffold Builder <bsb:index>` (BSB) is a framework that provides tools and pipelines
to reconstruct and simulate neural networks *in silico*.

BSB consumes a model description or ``configuration`` that describes the process of building a
neural circuit from scratch. This includes the creation of its different sub-regions, the placement
of the different cells that compose it and the connections of these cells to form an network. BSB
implements also simulators backends that allow for the simulation of the produced circuits.
Hence, the configuration can be extended to include simulation paradigms. Model configurations are
stored as dictionaries, and usually loaded from ``yaml`` or ``json`` files.

A modular configuration file to reconstruct cerebellar cortex circuit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configuration files will vary according to the specificities of each model, including the
specie, the sub-regions of interest, or the subject disease.

To help with the reproducibility of the different results from the `DBBS` we are here decomposing
the configuration files used in our models in sub-configurations that can be assembled and/or switch
on based on the model you wish to reproduce. This means that sub-configurations can be used to
override parameters (e.g. modifying sub-region size to match disease conditions) or introduce new
features to the circuit (e.g. adding a new cell type to the circuit).
For each configuration file, we also associate a full description of its parameters, including their
provenance and rationale behind them.

We have split our configuration based on the species studied. For each specie, we also produced a
default model or ``canonical circuit`` that will be used as a base for their alternatives. The other
sub-configurations (that can be used to override or extend the default one) are called
``extensions``.

Learn more about our configuration in the :doc:`configurations section <../configurations/configurations>`

Reconstruction of the rodent canonical circuit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, we present here the process to reconstruct the cerebellar cortex of the rodent brain,
based on the reconstruction of `De Schepper et al. 2022 <https://doi.org/10.1038/s42003-022-04213-y>`_
with BSB. This corresponds to our ``canonical circuit`` for the mouse so there is no need to assemble
multiple configuration files.

Assuming you are in the ``cerebellum`` folder, run the following command in your terminal:

.. code-block:: bash

    bsb compile configurations/mouse/mouse_cerebellar_cortex.yaml -v4 --clear

This command will produce a microcircuit of the mouse cerebellar cortex and store it in the
``cerebellum.hdf5`` file. This process might take a while depending on your machine.

Congratulations you have built your first cerebellar cortex model !

Simulation of the reconstructed circuits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A list of simulations can be attached to a BSB reconstruction through its configuration file (read
more :doc:`here <bsb:simulation/simulation>`). These simulations are dependant on the simulator you
choose and their paradigm. Some examples are provided for each simulator sustained by the
`cerebellum` package in the :doc:`configurations section <../configurations/configurations>`.
