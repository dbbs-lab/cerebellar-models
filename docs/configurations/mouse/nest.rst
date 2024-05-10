NEST
~~~~

Installation
^^^^^^^^^^^^
To reproduce the experiments presented below, you should install the NEST simulator following the
`installation instructions <https://nest-simulator.readthedocs.io/en/stable/installation/index.html>`_.
Do not forget to add the `nest_vars.sh` script to your `.bashrc` so that the installation of the
NEST modules work correctly at the end of the installation.

Nest modules are automatically compiled with BSB as ``components`` and deployed as ``cerebmodule``
through the provided configurations:

.. code-block:: yaml

    components:
      - cerebellum.nest_models.build_models
    simulations:
        simulation_name:
            simulator: nest
            modules:
                - cerebmodule

Alternatively, you can manually compile them, running the
`build_models.py <https://github.com/dbbs-lab/cerebellum/blob/master/cerebellum/nest_models/build_models.py>`_
python script.

Neuron models
^^^^^^^^^^^^^

All the default NEST neuron models are available for simulation with BSB. Additionally, BSB allows
for the different populations' parameters to be described as distribution or plain values:

.. code-block:: yaml

    simulations:
        simulation_name:
            simulator: nest
            seed: 1234
            cell_models:
              granule_cell:
                model: eglif_cond_alpha_multisyn
                constants:
                  t_ref: 1.5
                  I_e: -0.888
                  Vinit:
                    distribution: normal
                    mean: -62.0
                    std: 20.0

The ``distribution`` parameter correspond to one of the functions of the NEST
`random module <https://nest-simulator.readthedocs.io/en/stable/nest_behavior/random_numbers.html#the-nest-random-module>`_.
It can be useful to randomize certain parameters such as the initial membrane potential to avoid
synchrony of activity at the start of the simulation.

.. note::
    As in the previous configuration snippet, you can set the NEST random seed with BSB using the
    ``seed`` parameter.

Parrot neurons
##############

`Parrot neurons <https://nest-simulator.readthedocs.io/en/latest/models/parrot_neuron.html>`_ are
the most basic neuron model of NEST. They basically emits one spike for every incoming spike they
receive. We are using them here to represent the mossy fibers and glomerulus population as these
fibers only transmit spikes coming from other regions of the brain.

.. include:: ../../../cerebellum/nest_models/eglif_cond_alpha_multisyn.nestml
    :start-after: """
    :end-before: References

Neuron parameters
+++++++++++++++++

The parameters for the eglif models are extracted from Table 2 and Table 3 in Geminiani et al. (2019) [#geminiani_2019]_.
First, the default LIF parameters are as follows:

.. csv-table::
   :header-rows: 1
   :delim: ;

   Cell name;:math:`C_m\ (pF)`;:math:`\tau_m\ (ms)`;:math:`E_L\ (mV)`;:math:`t_{ref}\ (ms)`;:math:`V_{reset}\ (mV)`;:math:`V_{th}\ (mV)`; References
   Granule Cell; 7 (5.5 :math:`\pm` 0.5); 24.15 (24.15 :math:`\pm` 2); -62 (-62 :math:`\pm` 11); 1.5 (1.5 :math:`\pm` 0.4); -70 (-70); -41 (-41 :math:`\pm` 3);
   Golgi Cell; 145 (145 :math:`\pm` 73); 44 (44 :math:`\pm` 22); -62 (-62); 2 (2 :math:`\pm` 0.4); -75 (-75); -55 (-55 :math:`\pm` 1);
   Purkinje Cell; 334 (334 :math:`\pm` 106); 47 (47 :math:`\pm` 32); -59 (-59 :math:`\pm` 6); 0.5 (0.5 :math:`\pm` 0.1); -69 (-69); -43 (-43 :math:`\pm` 2);
   Basket and Stellate Cells; 14.6 (14.6); 9.125 (9.125); -68 (-68); 1.59 (1.59); -78 (-78); -53 (-53);

Then, the following parameters are optimized according to the method described in Geminiani et al.
(2018) [#geminiani_2018]_:

.. csv-table::
   :header-rows: 1
   :delim: ;

    Cell name;:math:`k_{adap}\ (nS \cdot ms^{-1})`;:math:`k_1\ (ms^{-1})`;:math:`k_2\ (ms^{-1})`;:math:`A_1\ (pA)`;:math:`A_2\ (pA)`;:math:`I_e\ (pA)`
    Granule Cell; 0.022; 0.311; 0.041; 0.01; -0.94; -0.888
    Golgi Cell; 0.217; 0.031; 0.023; 259.988; 178.01; 16.214
    Purkinje Cell; 1.491; 0.195; 0.041; 157.622; 172.622; 742.534
    Basket and Stellate Cells; 2.025; 1.887; 1.096; 5.953; 5.863; 3.711

It is not clear how the spiking parameters (i.e :math:`\lambda_0` and :math:`\tau_V`) are obtained
in the paper so the default values are used here.

Synapse models
^^^^^^^^^^^^^^

Static synapses
###############

Description
+++++++++++

By default, NEST `static synapses <https://nest-simulator.readthedocs.io/en/latest/models/static_synapse.html>`_
are used to connect the different neurons together. This model only transmit spikes as weights to
postsynaptic neurons after a provided delay.


Synapse parameters
++++++++++++++++++

The synaptic parameters used for the `canonical circuit` corresponds to the one listed in Table B of
supplementary document 1 in Geminiani et al. (2024) [#geminiani_2024]_. However, it is currently
unclear how these parameters were optimized, or which features were targeted:

.. csv-table::
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`
    Mf-glom;1;1
    glom-GrC;0.23;1
    glom-GoC;0.24;1
    GoC-GrC;-0.24;2
    GoC-GoC;-0.007;4
    GrC(aa)-GoC ;0.82;2
    GrC(aa)-PC;0.88;2
    GrC(pf)-GoC;0.05;5
    GrC(pf)-PC;0.14;5
    GrC(pf)-SC;0.18;5
    GrC(pf)-BC;0.1;5
    BC-PC;-0.44;4
    SC-PC;-1.64;5
    BC-BC ;-0.006;4
    SC-SC;-0.005;4

Simulation paradigms
^^^^^^^^^^^^^^^^^^^^

Different configuration files are available to reproduce experiments with the cerebellar cortex
circuit. As for the circuit reconstructions, the simulations are based on a basic paradigm (see the
following section) and can be extended with simulation extensions.

Basal activity
##############

The basal activity configuration file
`basal.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/nest/basal.yaml>`_
corresponds to the default activity of the cerebellar cortex circuit. Neurons are represented as
`eglif_cond_alpha_multisyn` and are connected with `static synapses`.

This simulation is set to last ``5000 ms`` (with a ``0.1 ms`` timestep) during which the neurons
are only stimulated with ``background noise`` represented as a ``4 Hz`` ``Poisson spike generator``
on the mossy fibers population:

.. code-block:: yaml

   simulations:
     basal_activity:
       simulator: nest
       resolution: 0.1
       duration: 5000
       modules:
        - cerebmodule
       seed: 1234
       devices:
        device: poisson_generator
        rate: 4
        targetting:
          strategy: cell_model
          cell_models:
            - mossy_fibers
        weight: 1
        delay: 0.1

Each neuron spiking activity is additionally recorded.

This simulation should demonstrate the activity of the network in a stable state. The results of
this simulation serves as a baseline for the following ones. To analyze these results,
for each neuron population, we define:

* the mean firing rate as the mean of each of its neuron's total number of
  spikes over the simulation time, expressed in Hz.
* the mean Inter-Spike intervals (ISI) as the mean of each of its neuron mean duration between each
  of its pair of consecutive spike, expressed in ms.

For this simulation, the mean firing rates and mean ISI obtained for each neuron population are as
follows (expressed in mean :math:`\pm` standard deviation):

.. csv-table::
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`3.9 \pm 0.99`; :math:`260 \pm 73`
    Granule cell; :math:`3.3 \pm 2.4`; :math:`470 \pm 490`
    Golgi cell;:math:`11 \pm 3.7`; :math:`110 \pm 53`
    Purkinje cell;:math:`50 \pm 2.6`; :math:`20 \pm 1.0`
    Basket cell;:math:`25 \pm 9.7`; :math:`46 \pm 17`
    Stellate cell;:math:`28 \pm 18`; :math:`57 \pm 40`

References
^^^^^^^^^^

.. include:: ../../../cerebellum/nest_models/eglif_cond_alpha_multisyn.nestml
    :start-after: start-references
    :end-before: See also

.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035
.. [#geminiani_2024] Geminiani, A., Casellato, C., Boele, H. J., Pedrocchi, A., De Zeeuw, C. I., &
   D’Angelo, E. (2024). Mesoscale simulations predict the role of synergistic cerebellar plasticity
   during classical eyeblink conditioning. PLOS Computational Biology, 20(4), e1011277.
   https://doi.org/10.1371/journal.pcbi.1011277