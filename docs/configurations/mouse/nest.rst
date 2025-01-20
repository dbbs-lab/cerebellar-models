NEST
~~~~

Installation
^^^^^^^^^^^^
To reproduce the experiments presented below, you should install the NEST simulator (see
instructions :doc:`here </getting-started/installation>`).

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

.. warning::
   The model described here is not matching the other LIF based models because of the
   sign in the membrane potential equations: the leak current should drive the membrane potential
   towards the resting state and not the opposite.

   Additionally, the model' s proof written in the paper seems to indicate that the model may
   produce unstable oscillations because of the freedom left to the adaptation parameters.
   For instance, the granule cell model described below is spiking without any input after several
   seconds.

   We have not being able to reproduce all the results of the Geminiani et al. (2018 and 2019)
   [#geminiani_2018]_ [#geminiani_2019]_, so we advise you to be careful when using it.

Neuron parameters
+++++++++++++++++
We consider here two types of mouse experiments that the `canonical circuit` model is reproducing
in simulation:

- The ``in-vitro`` experiments of a slice of isolated mouse cerebellar tissue.
- The ``in-vivo`` experiments of awake mouse.

Additionally, the ``in-vivo`` state is derived from the ``in-vitro`` state, as only some cells and
connections receive modifications.

`In-vitro` state
----------------

The parameters for the eglif models are extracted from Table 2 and Table 3 in Geminiani et al.
(2019) [#geminiani_2019]_. First, the default LIF parameters are as follows:

.. csv-table:: LIF neuron parameters
   :header-rows: 1
   :delim: ;

   Cell name;:math:`C_m\ (pF)`;:math:`\tau_m\ (ms)`;:math:`E_L\ (mV)`;:math:`t_{ref}\ (ms)`;:math:`V_{reset}\ (mV)`;:math:`V_{th}\ (mV)`
   Granule Cell; 7 (5.5 :math:`\pm` 0.5); 24.15 (24.15 :math:`\pm` 2); -62 (-62 :math:`\pm` 11); 1.5 (1.5 :math:`\pm` 0.4); -70 (-70); -41 (-41 :math:`\pm` 3)
   Golgi Cell; 145 (145 :math:`\pm` 73); 44 (44 :math:`\pm` 22); -62 (-62); 2 (2 :math:`\pm` 0.4); -75 (-75); -55 (-55 :math:`\pm` 1)
   Purkinje Cell; 334 (334 :math:`\pm` 106); 47 (47 :math:`\pm` 32); -59 (-59 :math:`\pm` 6); 0.5 (0.5 :math:`\pm` 0.1); -69 (-69); -43 (-43 :math:`\pm` 2)
   Basket and Stellate Cells; 14.6 (14.6); 9.125 (9.125); -68 (-68); 1.59 (1.59); -78 (-78); -53 (-53)

Then, the following parameters are optimized according to the method described in Geminiani et al.
(2018) [#geminiani_2018]_ :

.. csv-table:: EGLIF neuron parameters
   :header-rows: 1
   :delim: ;

    Cell name;:math:`k_{adap}\ (nS \cdot ms^{-1})`;:math:`k_1\ (ms^{-1})`;:math:`k_2\ (ms^{-1})`;:math:`A_1\ (pA)`;:math:`A_2\ (pA)`;:math:`I_e\ (pA)`
    Granule Cell; 0.022; 0.311; 0.041; 0.01; -0.94; -0.888
    Golgi Cell; 0.217; 0.031; 0.023; 259.988; 178.01; 16.214
    Purkinje Cell; 1.491; 0.195; 0.041; 157.622; 172.622; 590.0
    Basket and Stellate Cells; 2.025; 1.887; 1.096; 5.953; 5.863; 3.711

.. warning::
   It is not clear how the spiking parameters (i.e :math:`\lambda_0` and :math:`\tau_V` and initial :math:`V_m`)
   are obtained in the Geminiani et al. (2019) paper [#geminiani_2019]_ .
   These parameters were manually set to reproduce the F/I curves from the Figure 4 and Figure 3 from
   respectively Geminiani et al. (2018 and 2019) papers [#geminiani_2018]_ [#geminiani_2019]_.

.. warning::
   For the PC, we modified also the :math:`I_e` value so that the tonic firing rate of PC is ~45 Hz but
   maintained the F/I curve slope from the paper.

The postsynaptic currents are integrated to the soma with alpha exponential functions. Each function
is defined with a reversal potential parameter :math:`E_{rev}` and a time constant :math:`\tau_{syn}`.

These parameters depend on the connection types. In Nest, they are defined in the neuron equations.

The postsynaptic receptor parameters are listed in Table 2 of Geminiani et al. (2019b)
[#geminiani_2019b]_ :

.. _table-receptor:
.. csv-table:: Neuron Postsynaptic receptor parameters
   :header-rows: 1
   :delim: ;

   Cell name; Receptor id; :math:`E_{rev,i}\ (mV)`; :math:`\tau_{syn,i}\ (ms)`; Type
   Granule Cell; 1; 0; 5.8; exc.
   Granule Cell; 2; -80; 13.6; inh.
   Golgi Cell; 1; 0; 0.23; exc.
   Golgi Cell; 2; -80; 10; inh.
   Golgi Cell; 3; 0; 0.5; exc.
   Purkinje Cell; 1; 0; 1.1; exc.
   Purkinje Cell; 2; -80; 2.8; inh.
   Basket Cell; 1; 0; 0.64; exc.
   Basket Cell; 2; -80; 2; inh.

.. warning::
   The :math:`k_2` parameter should be greater than :math:`\dfrac{1}{\tau_m}` to prevent unstable
   oscillations of the membrane potential but the authors rounded down the values for the Granule
   cells which resulted in an unstable behavior. In `basal.yaml`, we therefore rounded up this
   value.

   On a side note, in the optimization section of the Geminiani et al. (2018) paper
   [#geminiani_2018]_, the authors wrote that the :math:`k_2` parameter should not be optimized but
   set to :math:`\dfrac{1}{\tau_m}` to have stable oscillations but this is not the case for most of
   the :math:`k_2` parameters listed in  Geminiani et al. (2019) paper [#geminiani_2019]_ .

`In-vivo` state
---------------

The parameters for the `in-vivo` state are the same as the `in-vitro` state, except for the PC for which
the endogenous current :math:`I_e` is set to 700 pA and :math:`\lambda_0` :math:`\tau_V` were changed to
increase the F/I curve slope. We target here ~80 Hz of tonic firing rate to match the range of Table 1
from Geminiani et al. 2024 [#geminiani_2024]_.

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

`In-vitro` state
----------------

The synaptic parameters used for the `canonical circuit` corresponds to the one listed in Table B of
supplementary document 1 in Geminiani et al. (2024) [#geminiani_2024]_. The receptor id corresponds
to the postsynaptic receptor used for the connection (see table :ref:`table-receptor`).

.. csv-table:: Presynaptic parameters
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    Mf-glom;1;1;1
    glom-GrC;0.23;1;1
    glom-GoC;0.24;1;1
    GoC-GrC;0.24;2;2
    GoC-GoC;0.007;4;2
    GrC(aa)-GoC ;0.82;2;3
    GrC(aa)-PC;0.2;2;1
    GrC(pf)-GoC;0.05;5;3
    GrC(pf)-PC;0.05;5;1
    GrC(pf)-SC;0.05;5;1
    GrC(pf)-BC;0.04;5;1
    BC-PC;0.44;4;2
    SC-PC;0.17;5;2
    BC-BC ;0.006;4;2
    SC-SC;0.005;4;2

.. warning::
   It is currently unclear from the paper, how the synaptic parameters were optimized, or which features were targeted.

.. warning::
   In our experiments, we decreased the weights for the pf-SC, pf-BC so that the activity of MLI lies around
   ~15 Hz for both BC and SC. Then aa-PC, pf-PC were decreased to maintain the PC in a stable low activity ~50Hz.
   Finally, the SC-PC was scaled to take into account the increase of synapses from the connectivity rule.

`In-vivo` state
---------------

The parameters for the `in-vivo` state are the same as the `in-vitro` state, except for some of the connections:

.. csv-table:: In-vivo Presynaptic parameters
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    GrC(pf)-PC;0.14;5;1
    GrC(aa)-PC;0.41;2;1
    GrC(pf)-SC;0.08;5;1
    GrC(pf)-BC;0.06;5;1
    BC-PC;0.8;4;2

Our target was to match the Table 1 from Geminiani et al. 2024 [#geminiani_2024]_.

Simulation paradigms
^^^^^^^^^^^^^^^^^^^^

Different configuration files are available to reproduce experiments with the cerebellar cortex
circuit. As for the circuit reconstructions, the simulations are based on a basic paradigm (see the
following section) and can be extended with simulation extensions.

.. _basal-activity:

Basal activity
##############

The basal activity configuration file
`basal.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/nest/basal.yaml>`_
implements to the default activity of the cerebellar cortex circuit. Neurons are represented as
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

* the mean firing rate as the mean of each of its (spiking at least once) neurons' total number of
  spikes over the simulation time, expressed in Hz.
* the mean Inter-Spike intervals (ISI) as the mean of each of its (spiking at least twice) neurons
  mean duration between each of its pair of consecutive spike, expressed in ms.

For this simulation, the mean firing rates and mean ISI obtained for each neuron population are as
follows (expressed in mean :math:`\pm` standard deviation):

`In-vitro` state
++++++++++++++++

.. csv-table:: Results of the canonical circuit in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.0 \pm 1.4`; :math:`250 \pm 140`
    Granule cell; :math:`3.4 \pm 3.4`; :math:`330 \pm 270`
    Golgi cell;:math:`11 \pm 5.4`; :math:`120 \pm 94`
    Purkinje cell;:math:`47 \pm 1.1`; :math:`22 \pm 0.48`
    Basket cell;:math:`16 \pm 7.7`; :math:`77 \pm 45`
    Stellate cell;:math:`13 \pm 11`; :math:`180 \pm 240`

`In-vivo` state
+++++++++++++++

.. csv-table:: Results of the canonical circuit in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.0 \pm 1.4`; :math:`250 \pm 140`
    Granule cell; :math:`3.4 \pm 3.4`; :math:`330 \pm 270`
    Golgi cell;:math:`11 \pm 5.4`; :math:`120 \pm 94`
    Purkinje cell;:math:`91 \pm 2.2`; :math:`11 \pm 0.27`
    Basket cell;:math:`21 \pm 9.0`; :math:`56 \pm 25`
    Stellate cell;:math:`18 \pm 14`; :math:`120 \pm 140`

Mossy fiber stimulus
####################

The mossy fiber stimulus configuration file
`stimulus_mossy.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/nest/stimulus_mossy.yaml>`_
adds to the basal activity configuration file another simulation configuration ``mf_stimulus``. ``mf_stimulus`` has the
same parameters as the ``basal_activity`` but with a stimulus of the mossy fibers (see :ref:`basal-activity`).

On top of the basal paradigm, we introduce here a ``stimulus`` represented as a ``150 Hz``
``Poisson spike generator`` between ``1200`` and ``1250`` ms. This latter targets the
``mossy_fibers`` that are within a ``sphere`` of radius ``90`` :math:`\mu m` and a center at
``(150.0, 65.0, 100.0)``. This mimic the integration of an input on the mossy fibers.

For this simulation, **during the stimulus**, the mean firing rates and mean ISI obtained for each
neuron population are as follows (expressed in mean :math:`\pm` standard deviation):

`In-vitro` state
++++++++++++++++

.. csv-table:: Results of the canonical circuit during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`45 \pm 74`; :math:`6.5 \pm 3.`
    Granule cell; :math:`21 \pm 46`; :math:`9.1 \pm 7.3`
    Golgi cell;:math:`52 \pm 45`; :math:`17 \pm 12`
    Purkinje cell;:math:`59 \pm 11`; :math:`17 \pm 3.7`
    Basket cell;:math:`59 \pm 88`; :math:`12 \pm 7.3`
    Stellate cell;:math:`55 \pm 98`; :math:`14 \pm 7.5`

`In-vivo` state
+++++++++++++++

.. csv-table:: Results of the canonical circuit during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`45 \pm 74`; :math:`6.5 \pm 3.`
    Granule cell; :math:`21 \pm 46`; :math:`9.1 \pm 7.3`
    Golgi cell;:math:`52 \pm 45`; :math:`17 \pm 12`
    Purkinje cell;:math:`140 \pm 33`; :math:`7.1 \pm 1.9`
    Basket cell;:math:`77 \pm 54`; :math:`8.7 \pm 5.9`
    Stellate cell;:math:`71 \pm 55`; :math:`9.6 \pm 6.2`

References
^^^^^^^^^^

.. include:: ../../../cerebellum/nest_models/eglif_cond_alpha_multisyn.nestml
    :start-after: start-references
    :end-before: See also

.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035
.. [#geminiani_2019b] Geminiani, A., Pedrocchi, A., D’Angelo, E., & Casellato, C. (2019). Response
   dynamics in an olivocerebellar spiking neural network with non-linear neuron properties.
   Frontiers in computational neuroscience, 13, 68.
   https://doi.org/10.3389/fncom.2019.00068
.. [#geminiani_2024] Geminiani, A., Casellato, C., Boele, H. J., Pedrocchi, A., De Zeeuw, C. I., &
   D’Angelo, E. (2024). Mesoscale simulations predict the role of synergistic cerebellar plasticity
   during classical eyeblink conditioning. PLOS Computational Biology, 20(4), e1011277.
   https://doi.org/10.1371/journal.pcbi.1011277
