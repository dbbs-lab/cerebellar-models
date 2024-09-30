Inferior olivary
~~~~~~~~~~~~~~~~
The inferior olive is a nucleus in the brainstem with neurons that exhibit continuous sub-threshold activity.
It provides one of the two main inputs to the cerebellum: the so-called climbing fibers. Activation of these climbing
fibers is generally thought to be involved in the timing of motor commands and/or motor learning. Climbing fiber
activation triggers large all-or-none action potentials in cerebellar Purkinje cells, which override any ongoing
activity and temporarily silence the cells. Empirical evidence indicates that climbing fibers can transmit a short burst
of spikes following an olivary cell somatic spike, potentially increasing the amount of information transferred to the
cerebellum with each activation [#de_gruijl_2012]_.
The default configuration with IO is implemented in `dcn_io.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn_io.yaml>`_.

Configuration
^^^^^^^^^^^^^
In `dcn_io.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn_io.yaml>`_ , a new region called ``inferior_olivary`` was added to the ``canonical circuit + DCN``
(see `dcn.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn.yaml>`_ for more microcircuital info).
This region was defined of type ``group`` because in it there is only the ``io layer``. ``io layer`` has a thickness of ``100 µm`` .
Additionally, in order to ensure that ``inferior_olivary`` are placed under the ``cerebellar_nuclei``, the ``origin``
of the ``dcn_layer`` was set to ``[0,0,100]`` and the ``origin`` of the ``granular_layer`` was updated to ``[0,0,300]``.

Cell types
++++++++++
No morphologies are currently available for inferior olivary (IO) neurons, so they are represented only by their soma.
We considered a single population of IO neurons.
The number of IO neurons to be placed is estimated based on the ratio between the number of Purkinje cells and
IO itself, which is reported to be ``5:1`` in Blatt and Eisenman [#blatt_1985]_.

Placement
+++++++++
IO neurons are assumed to be uniformly distributed in their own layer, hence the bsb ``RandomPlacement`` strategy is
chosen to place them.

Connectivity
++++++++++++

.. csv-table::
   :header-rows: 1
   :delim: ;

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; Specifics; References
   22; IO; /; PC; / ; :ref:`fix_in`;``indegree`` =1; Geminiani et al. (2024) [#geminiani_2024]_
   23; IO; /; SC; / ; :ref:`io_mli`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   23; IO; /; BC; / ; :ref:`io_mli`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   24; IO; / ; DCNp ; / ; :ref:`all_to_all`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   25; IO; /; DCNi; / ; :ref:`all_to_all`; / ; Geminiani et al. (2024) [#geminiani_2024]_
   26; DCNi; / ; IO ; / ; :ref:`all_to_all`; / ; Geminiani et al. (2024) [#geminiani_2024]_

NEST simulation
^^^^^^^^^^^^^^^

Neuron parameters
+++++++++++++++++
IO population was represented as an EGLIF point neuron model (see :doc:`NEST section <nest>`). Parameters sets for IO neurons
are taken from Geminiani et al (2019) [#geminiani_2019]_. The default LIF parameters are the following:

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

It is not clear how the spiking parameters are obtained in the Geminiani et al. (2019) paper [#geminiani_2019]_. The values were extracted from a BSB configuration provided by the authors.

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

.. csv-table:: Presynaptic parameters for IO connections
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    IO-PC; 300; 4;3
    IO-MLI; 2.5; 40 ; 3
    IO-DCNp; 2.5; 4; 1
    IO-DCNi; 0.1; 5; 1
    DCNi-IO; 0.75; 25; 2

Simulation paradigms
++++++++++++++++++++

The `dcn_io_nest.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn_io_nest.yaml>`_ are
including all the simulation paradigms described in the :doc:`NEST section <nest>`) but include the DCN and IO cells in the
circuit.

Basal activity
##############

No basal activity changes are observed in the cerebellar network beacause IO presents no autorhythm [#de_gruijl_2012]_
[#lefler_2013]_.

Stimulation protocol
####################



References
^^^^^^^^^^

.. [#de_gruijl_2012] De Gruijl, J. R., Bazzigaluppi, P., de Jeu, M. T., & De Zeeuw, C. I. (2012). Climbing fiber burst size and olivary sub-threshold oscillations in a network setting. PLoS computational biology, 8(12), e1002814.
.. [#blatt_1985] Blatt, G. J., & Eisenman, L. M. (1985). A qualitative and quantitative light microscopic study of the inferior olivary complex of normal, reeler, and weaver mutant mice. Journal of Comparative Neurology, 232(1), 117-128.
.. [#geminiani_2024] Geminiani, Alice, et al. "Mesoscale simulations predict the role of synergistic cerebellar plasticity during classical eyeblink conditioning." PLOS Computational Biology 20.4 (2024): e1011277. https://doi.org/10.1371/journal.pcbi.1011277.
.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035
.. [#geminiani_2019b] Geminiani, A., Pedrocchi, A., D’Angelo, E., & Casellato, C. (2019). Response
   dynamics in an olivocerebellar spiking neural network with non-linear neuron properties.
   Frontiers in computational neuroscience, 13, 68.
.. [#lefler_2013] Lefler, Y., Torben-Nielsen, B., & Yarom, Y. (2013). Oscillatory activity, phase differences, and phase resetting in the
   inferior olivary nucleus. Frontiers in systems neuroscience, 7, 22.