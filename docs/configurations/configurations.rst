Circuit configurations
======================

The ``configurations`` folder contains configuration files for cerebellar reconstruction and
simulation with BSB.

These configuration files are split into subdirectories according to their respective species. Each
subdirectory also contains documentation explaining the parameters of the configuration files
(i.e., their rationale and sources) as well as their extensions.

By default, the base configuration covers the reconstruction of the circuit with BSB. Simulation
parameters are distributed by simulator in state-specific subdirectories and are selectable via the
CLI (see the :doc:`CLI section </cli>` for details).

The following is the list of the species for which canonical configurations have been implemented.

.. toctree::
    :maxdepth: 2
    :caption: Configurations by species
    :glob:

    mouse/mouse