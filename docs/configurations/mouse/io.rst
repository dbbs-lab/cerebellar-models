Inferior olivary
~~~~~~~~~~~~~~~~
The inferior olive (IO) is a nucleus in the brainstem with neurons that exhibit continuous sub-threshold activity.
It provides one of the two main inputs to the cerebellum: the so-called climbing fibers (cf). Activation of the cf
is generally thought to be involved in the timing of motor commands and/or motor learning. cf activation triggers
large all-or-none action potentials in cerebellar Purkinje cells (PC), which override any ongoing activity and
temporarily silence the cells. Empirical evidence indicates that cf can transmit a short burst of spikes following
an IO cell somatic spike, potentially increasing the amount of information transferred to the cerebellum with each
activation [#de_gruijl_2012]_.
The default configuration with IO is implemented in
`dcn_io.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/dcn-io/dcn_io.yaml>`_.

Configuration
^^^^^^^^^^^^^
In `dcn_io.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/dcn-io/dcn_io.yaml>`_ ,
a new region called ``inferior_olivary`` was added to the ``canonical circuit + DCN`` model
(see `dcn.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/dcn-io/dcn.yaml>`_
and :doc:`DCN section <dcn>` for more microcircuital info).
This region contains only one ``Layer`` Partition: ``io layer``. ``io layer`` has a thickness of ``100 µm`` .
Additionally, to ensure that ``inferior_olivary`` are placed under the ``cerebellar_nuclei``, the ``origin``
of the ``dcn_layer`` was set to ``[0, 0, 100]`` and the ``origin`` of the ``granular_layer`` was updated to ``[0, 0, 300]``.

Cell types
++++++++++
No morphologies are currently available for IO neurons, so they are modelled as point neurons.
We considered a single population of IO neurons.
The number of IO neurons to be placed is estimated based on the ratio between the number of PC and
IO itself, which is reported to be ``5:1`` in Blatt and Eisenman [#blatt_1985]_.

Placement
+++++++++
IO neurons are assumed to be uniformly distributed in their own layer, hence the bsb ``RandomPlacement``
strategy is chosen to place them.

Connectivity
++++++++++++

.. csv-table::
   :header-rows: 1
   :delim: ;

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; Specifics; References
   22; IO; /; PC; / ; :ref:`fix_in`;``indegree`` = 1; Geminiani et al. (2024) [#geminiani_2024]_
   23; IO; /; SC; / ; :ref:`io_mli`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   23; IO; /; BC; / ; :ref:`io_mli`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   24; IO; / ; DCNp ; / ; :ref:`fix_out`; ``outdegree`` = 60 ; Geminiani et al. (2024) [#geminiani_2024]_
   25; IO; /; DCNi; / ; :ref:`fix_out`; ``outdegree`` = 33 ; Geminiani et al. (2024) [#geminiani_2024]_
   26; DCNi; / ; IO ; / ; :ref:`fix_in`; ``indegree`` = 33 ; Geminiani et al. (2024) [#geminiani_2024]_

NEST simulation
^^^^^^^^^^^^^^^

As for the cerebellar cortex, we differentiate parameters for the ``in-vitro`` and ``awake`` states.

Neuron parameters
+++++++++++++++++
IO population was represented as an EGLIF point neuron model (see :doc:`NEST section <nest>`).
Parameters sets for IO neurons are taken from Geminiani et al (2019) [#geminiani_2019]_.
The IO neuron parameters are the same for the `in-vitro` and awake state because we do not have a reference
parameter set for this cell.
The default LIF parameters are reported below:

.. csv-table:: LIF neuron parameters for IO
   :header-rows: 1
   :delim: ;

   Cell name;:math:`C_m\ (pF)`;:math:`\tau_m\ (ms)`;:math:`E_L\ (mV)`;:math:`t_{ref}\ (ms)`;:math:`V_{reset}\ (mV)`;:math:`V_{th}\ (mV)`
   IO; 189 (189 :math:`\pm` 12); 11 (11 :math:`\pm` 4); -45 (-45); 1 (1); -45 (-45); -35 (-35)

Then, the following parameters are optimized according to the method described in Geminiani et al. (2019) [#geminiani_2019]_ :

.. csv-table:: EGLIF neuron parameters for IO
   :header-rows: 1
   :delim: ;

    Cell name;:math:`k_{adap}\ (nS \cdot ms^{-1})`;:math:`k_1\ (ms^{-1})`;:math:`k_2\ (ms^{-1})`;:math:`A_1\ (pA)`;:math:`A_2\ (pA)`;:math:`I_e\ (pA)`
    IO; 1.928; 0.191; 0.091; 1810.923; 1358.197; -18.101

.. warning::
    It is not clear how the spiking parameters are obtained in the Geminiani et al. (2019) paper [#geminiani_2019]_.
    The values were extracted from a BSB configuration provided by the authors.

The postsynaptic receptors are defined as listed in Table 2 of Geminiani et al. (2019b) [#geminiani_2019b]_:

.. _io-table-receptor:
.. csv-table:: IO Postsynaptic receptor parameters
   :header-rows: 1
   :delim: ;

   Cell name; Receptor id; :math:`E_{rev,i}\ (mV)`; :math:`\tau_{syn,i}\ (ms)`; Type
   IO; 1; 0; 1; exc.
   IO; 2; -80; 60; inh.

Synapse parameters
++++++++++++++++++
IO connections are represented as ``static synapses`` (see :doc:`NEST section <nest>`). The receptor ids correspond to
the postsynaptic receptors used for the connections.
It is still unclear from the references how these parameters were optimized.

.. warning::
   The reported values were manually adjusted through trial and error to ensure a reasonable excitation/inhibition ratio
   on IO target populations.

`In-vitro` state
----------------

.. csv-table:: Presynaptic parameters for IO connections
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    IO-PC; 0.6; 4;3
    IO-BC; 5.0; 40 ; 3
    IO-SC; 6.5; 40 ; 3
    IO-DCNp; 0.5; 4; 1
    IO-DCNi; 0.25; 5; 1
    DCNi-IO; 0.45; 25; 2

Awake state
-----------

.. csv-table:: Presynaptic parameters for IO connections
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    IO-PC; 1.6; 4;3
    IO-BC; 5.0; 40 ; 3
    IO-SC; 5.0; 40 ; 3
    IO-DCNp; 0.4; 4; 1
    IO-DCNi; 0.25; 5; 1
    DCNi-IO; 0.45; 25; 2

Simulation paradigms
++++++++++++++++++++

The `dcn_io_nest.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn_io_nest.yaml>`_ are
including all the simulation paradigms described in the :doc:`NEST section <nest>`) but include the DCN and IO cells in the
circuit.

Basal activity
--------------

No basal activity changes are observed in the cerebellar network beacause IO presents no autorhythm [#de_gruijl_2012]_
[#lefler_2013]_.

Stimulation protocol
--------------------

To test the functionality of the entire olivocerebellar network, another stimulation protocol was used. It
simulates the Eyeblink Classical Conditioning, a Pavlovian conditioning consisting in a conditioned stimulus (CS),
typically a light, paired with an unconditioned stimulus (US), usually an air puff to the eye.
According to Geminiani et al., 2024 [#geminiani_2024]_, a CS of ``40 Hz``  arrives on ``mossy_fibers`` in the
interval ``[1000, 1250] ms``, while a US of ``500 Hz`` arrives as a burst on ``io`` in the interval ``[1250, 1260] ms``.

`In-vitro` state
################

.. csv-table:: Results of the canonical circuit in `in-vitro` state with DCN and IO during stimulus of the mf and the IO
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz) [CS]; Mean ISI (ms) [CS]; Mean Firing rate (Hz) [CS+US]; Mean ISI (ms) [CS+US]
    Mossy cell; :math:`44 \pm 13`; :math:`22 \pm 8.0`; :math:`38 \pm 60`; :math:`3.1 \pm 2.3`
    Granule cell; :math:`21 \pm 30`; :math:`36 \pm 36`; :math:`20 \pm 48`; :math:`6.1 \pm 1.1`
    Golgi cell;:math:`42 \pm 17`; :math:`28 \pm 13`; :math:`26 \pm 50`; :math:`8.7 \pm 0.2`
    Purkinje cell;:math:`55 \pm 5.9`; :math:`18.0 \pm 2.0`; :math:`140 \pm 100`; :math:`3.1 \pm 2.8`
    Basket cell;:math:`66 \pm 24`; :math:`19 \pm 9.5`; :math:`61 \pm 49`; not enough spikes per neuron
    Stellate cell;:math:`46 \pm 43`; :math:`31 \pm 31`; :math:`48 \pm 57`; :math:`7.2 \pm 1.1`
    DCNp; :math:`30 \pm 4.7`; :math:`33 \pm 5.2`; :math:`33 \pm 47`; not enough spikes per neuron
    DCNi; :math:`11 \pm 2.3`; :math:`93 \pm 11`;  :math:`12 \pm 32`; not enough spikes per neuron
    IO; 0; no spikes; :math:`270 \pm 120`; :math:`1.9 \pm 0.5`

Awake state
###########

.. csv-table:: Results of the canonical circuit in awake state with DCN and IO during stimulus of the mf and the IO
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz) [CS]; Mean ISI (ms) [CS]; Mean Firing rate (Hz) [CS+US]; Mean ISI (ms) [CS+US]
    Mossy cell; :math:`43 \pm 13`; :math:`23 \pm 8.1`; :math:`46 \pm 68`; :math:`2.8 \pm 1.8`
    Granule cell; :math:`22 \pm 35`; :math:`34 \pm 37`; :math:`22 \pm 50`; :math:`6.4 \pm 1.2`
    Golgi cell;:math:`42 \pm 21`; :math:`29 \pm 15`; :math:`22 \pm 45`; :math:`7.0 \pm 0.0`
    Purkinje cell;:math:`140 \pm 14`; :math:`7.1 \pm 0.71`; :math:`300 \pm 150`; :math:`2.7 \pm 1.9`
    Basket cell;:math:`83 \pm 49`; :math:`19 \pm 21`; :math:`55 \pm 54`; :math:`6.3 \pm 1.1`
    Stellate cell;:math:`74 \pm 68`; :math:`22 \pm 24`; :math:`57 \pm 81`; :math:`5.6 \pm 1.3`
    DCNp; :math:`59 \pm 7.6`; :math:`17 \pm 2.1`; :math:`65 \pm 56`; :math:`6.0 \pm 2.7`
    DCNi; :math:`11 \pm 2.0`; :math:`90 \pm 13`; :math:`9.1 \pm 29`; not enough spikes per neuron
    IO; 0; no spikes; :math:`350 \pm 150`; :math:`1.9 \pm 0.78`

References
^^^^^^^^^^

.. [#de_gruijl_2012] De Gruijl, J. R., Bazzigaluppi, P., de Jeu, M. T., & De Zeeuw, C. I. (2012).
   "Climbing fiber burst size and olivary sub-threshold oscillations in a network setting."
   PLoS computational biology, 8(12), e1002814.
   https://doi.org/10.1371/journal.pcbi.1002814
.. [#blatt_1985] Blatt, G. J., & Eisenman, L. M. (1985).
   "A qualitative and quantitative light microscopic study of the inferior olivary complex of normal, reeler,
   and weaver mutant mice." Journal of Comparative Neurology, 232(1), 117-128.
   https://doi.org/10.1002/cne.902320110
.. [#geminiani_2024] Geminiani, Alice, et al.
   "Mesoscale simulations predict the role of synergistic cerebellar plasticity during classical eyeblink conditioning."
   PLOS Computational Biology 20.4 (2024): e1011277.
   https://doi.org/10.1371/journal.pcbi.1011277.
.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   "Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models."
   Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035
.. [#geminiani_2019b] Geminiani, A., Pedrocchi, A., D’Angelo, E., & Casellato, C. (2019).
   "Response dynamics in an olivocerebellar spiking neural network with non-linear neuron properties."
   Frontiers in computational neuroscience, 13, 68.
   https://doi.org/10.3389/fncom.2019.00068
.. [#lefler_2013] Lefler, Y., Torben-Nielsen, B., & Yarom, Y. (2013).
   "Oscillatory activity, phase differences, and phase resetting in the inferior olivary nucleus."
   Frontiers in systems neuroscience, 7, 22.
   https://doi.org/10.3389/fnsys.2013.00022
