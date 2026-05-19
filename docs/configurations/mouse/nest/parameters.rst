NEST models and simulation paradigms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

States
^^^^^^
We consider here two types of mouse experiments that the `canonical circuit` and each derived
model is reproducing in simulation:

- The ``in-vitro`` experiments of a slice of isolated mouse cerebellar tissue.
- The ``awake`` experiments of behaving mouse.

Additionally, the ``awake`` state is derived from the ``in-vitro`` state, as only some cells and
connections receive modifications.

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

Geminiani EGLIF
###############

This model is based on the work of Geminiani et al. [#geminiani_2018]_ [#geminiani_2019]_ [#geminiani_2024]_.
Details about the model are available in :doc:`this page <eglif_cond_alpha_multisyn>`.


De Grazia EGLIF
###############

This model is based on the work of De Grazia et al. [#degrazia_2026]_.
Details about the model are available in :doc:`this page <eglif_cond_alpha_multisyn>`.

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
to the postsynaptic receptor used for the connection.
However, we made the following adjustments to obtain results close to the paper and literature:

- In our experiments, we decreased the weights for the pf-SC, pf-BC so that the activity of MLI lies around
  ~15 Hz for both BC and SC [#kim_2021]_. Then aa-PC, pf-PC were decreased to maintain the PC in a stable
  low activity ~50Hz [#telgkamp_2002]_.
- Synaptic parameter values of DCNp, DCNi and IO were manually adjusted through trial and error to ensure a
  reasonable excitation/inhibition activity ratio.
- Finally, the SC-PC was scaled to take into account the increase of synapses from the connectivity rule.

.. warning::
   It is currently unclear from the paper, how the synaptic parameters were optimized, or which features were targeted.

Awake state
-----------

The parameters for the awake state are the same as the `in-vitro` state, except for the connections
GrC(pf)-PC, GrC(aa)-PC, GrC(pf)-SC, GrC(pf)-BC, BC-PC for which we adapted the weights to match the Table 1
from Geminiani et al. 2024 [#geminiani_2024]_.

Tsodyks Markram Synapse
#######################

Description
+++++++++++

The `Tsodyks-Markram synapse <https://nest-simulator.readthedocs.io/en/latest/models/tsodyks2_synapse.html>`_ synapse
model implements synaptic short-term depression and short-term facilitation according to
Tsodyks et al. [#tsodyks_1997]_ and Fuhrman et al. [#fuhrman_2002]_.
This connection model merely scales the synaptic weight, based on the spike history parameters of the kinetic model.

Synapse parameters
++++++++++++++++++

For each synapse of the `canonical circuit`, the initial value of ``u`` was set to ``0`` and ``x`` to ``1.0``.

`In-vitro` state
----------------

The synaptic parameters used for the canonical circuit correspond to those 
listed in the table below, obtained from Masoli et al. (2022) [#masoli_2022]_ .
The receptor ID corresponds to the postsynaptic receptor used for the connection.
The weights have been rescaled under the assumption
that the first peak of the postsynaptic conductance (:math:`g_{syn_0}`) for the
Tsodyks–Markram synapse must have the same amplitude as the ones obtained with a static_synapse model.

:math:`weight_{tsodyks} = \dfrac{{weight_{static}}^2}{g_{syn_0}}`

.. Note::
   The connections mf-glom and GoC-GoC are both considered static since, for these two connections,
   we do not have Tsodyks-Markram parameters.
   Moreover, for the pf-SC connection, the weight was adjusted manually to keep the firing rate
   within the desired range.


Awake state
-----------

The parameters for the awake state are the same as the in-vitro state, except for the following connections:
GrC(pf)-PC, GrC(aa)-PC, GrC(pf)-SC, GrC(pf)-BC, BC-PC for which we adapted the weights to match the Table 1
from Geminiani et al. 2024 [#geminiani_2024]_.

.. Note:: 
   For the simulations using Tsodyks-Markram synapse, the mean firing rates and mean interspike intervals (ISI)
   obtained for each neuron population from both in-vitro and awake states are expected to be the same,
   as the ones obtained with static synapses.
   For pf-SC connection weight was adjusted manually to keep the firing rate in the desired range.

.. _nest-paradigms:

Simulation paradigms
^^^^^^^^^^^^^^^^^^^^

Different configuration files are available to reproduce experiments with the cerebellar cortex
circuit. As for the circuit reconstructions, the simulations are based on a basic paradigm (see the
following section) and can be extended with simulation extensions.

To analyze the spiking results of the following simulations, for each neuron population, we define:

* the mean firing rate as the mean of each of its (spiking at least once) neurons' total number of
  spikes over the simulation time, expressed in Hz.
* the mean Inter-Spike intervals (ISI) as the mean of each of its (spiking at least twice) neurons
  mean duration between each of its pair of consecutive spike, expressed in ms.

.. _basal-activity:

Basal activity
##############

The basal activity configuration file
`basal.yaml <https://github.com/dbbs-lab/cerebellar-models/blob/master/configurations/mouse/in-vitro/nest/basal_vitro.yaml>`_
implements to the default activity of the cerebellar cortex circuit.

This simulation is set to last ``5000 ms`` (with a ``0.1 ms`` timestep) during which the neurons
are only stimulated with ``background noise`` represented as a ``4 Hz`` ``Poisson spike generator``
on the mossy fibers population. Each neuron spiking activity is additionally recorded.

This simulation should demonstrate the activity of the network in a stable state. The results of
this simulation serves as a baseline for the following ones.

Mossy fiber stimulus
####################

The mossy fiber stimulus configuration file
`stimulus_mossy.yaml <https://github.com/dbbs-lab/cerebellar-models/blob/master/configurations/mouse/in-vitro/nest/mf_stimulus.yaml>`_
is derived from the ``basal_activity`` but adds a stimulus of the mossy fibers (see :ref:`basal-activity`).

On top of the basal paradigm, we introduce here a ``stimulus`` represented as a ``150 Hz``
``Poisson spike generator`` between ``1200`` and ``1250`` ms. This latter targets the
``mossy_fibers`` that are within a ``sphere`` of radius ``90`` :math:`\mu m` and a center at
``(150.0, 65.0, 100.0)``. This mimic the integration of an input on the mossy fibers.

Eyeblink Classical Conditioning
###############################

To test the functionality of the entire olivocerebellar network
(cerebellar cortex + deep cerebellar nuclei + inferior olive), another stimulation protocol was used. It
simulates the Eyeblink Classical Conditioning, a Pavlovian conditioning consisting in a conditioned stimulus (CS),
typically a light, paired with an unconditioned stimulus (US), usually an air puff to the eye.
According to Geminiani et al., 2024 [#geminiani_2024]_, a CS of ``40 Hz``  arrives on ``mossy_fibers`` in the
interval ``[1000, 1250] ms``, while a US of ``500 Hz`` arrives as a burst on ``io`` in the interval ``[1250, 1260] ms``.

Simulation results
^^^^^^^^^^^^^^^^^^
Simulation results depends on"which cells types were included in the `canonical circuit`,
on the state, the simulation paradigm and models used for the neuron and the synapses  parameters.
The results for the `eglif_cond_alpha_multisyn` with `static synapses` are available in :doc:`this page <simulation-results>`

References
^^^^^^^^^^

.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_cond_alpha_multisyn.nestml
    :start-after: start-references
    :end-before: See also

.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_multirec.nestml
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
.. [#telgkamp_2002] Telgkamp, P., & Raman, I. M. (2002).
   Depression of inhibitory synaptic transmission between Purkinje cells and neurons of the cerebellar nuclei.
   Journal of Neuroscience, 22(19), 8447-8457.
   https://doi.org/10.1523/JNEUROSCI.22-19-08447.2002.
.. [#kim_2021] Kim, J., & Augustine, G. J. (2021).
   Molecular layer interneurons: key elements of cerebellar network computation and behavior.
   Neuroscience, 462, 22-35.
   https://doi.org/10.1016/j.neuroscience.2020.10.008
.. [#tsodyks_1997] Tsodyks MV,  Markram H (1997). The neural code between neocortical
   pyramidal neurons depends on neurotransmitter release probability.
   PNAS, 94(2):719-23.
   https://doi.org/10.1073/pnas.94.2.719
.. [#fuhrman_2002] Fuhrman, G, Segev I, Markram H, Tsodyks MV (2002). Coding of
   temporal information by activity-dependent synapses. Journal of
   Neurophysiology, 87(1):140-8.
   https://doi.org/10.1152/jn.00258.2001
.. [#masoli_2022] Masoli, S., Rizza, M. F., Tognolina, M., Prestori, F., & D’Angelo, E. (2022).
   Computational models of neurotransmission at cerebellar synapses unveil the impact on network computation. 
   Frontiers in Computational Neuroscience, 16, 1006989.
   https://doi.org/10.3389/fncom.2022.1006989
