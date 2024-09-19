Deep Cerebellar Nuclei
~~~~~~~~~~~~~~~~~~~~~~

Deep cerebellar nuclei (DCN) are the primary output structures of the cerebellar cortex [#dangelo_2018]_.
They receive inhibitory input from Purkinje cells and excitatory input from mossy fibers and climbing fibers (from IO).
DCN are composed by three distinct nuclei: dentate nucleus, fastigial nucleus and interposed nucleus.
The default configuration with DCN is implemented in `dcn.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn.yaml>`_.


Configuration
^^^^^^^^^^^^^
In `dcn.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn.yaml>`_ , a new region called ``cerebellar_nuclei`` was added to the ``canonical circuit``. This region was defined of type ``group`` because in it there is only the ``dcn layer``. ``dcn layer`` has a thickness of ``200 µm`` . Additionally, in order to ensure that ``cerebellar_nuclei`` are placed under the ``cerebellar_cortex``, the ``origin`` of the ``granular_layer`` was set to ``[0,0,200]``.

Cell types
++++++++++
For DCN, two types of neurons are considered [#uusisaari_2008]_ [#geminiani_2019b]_:

* **DCNp**: they are the excitatory neurons projecting outside the cerebellum to various brain regions, including the thalamus, the red nucleus, the vestibular nuclei, and the reticular formation;
* **DCNi**: they are GABAergic interneurons which send inhibitory feedback to IO.

.. warning::
   In Geminiani et al (2019) [#geminiani_2019]_, DCN populations are defined in a different way:

   * **DCNnL** for the excitatory population (DCNp)
   * **DCNp** for the inhibitory one (DCNi)

   Then, starting from Geminiani et al (2019b) [#geminiani_2019b]_ , names for DCN populations were redefined as reported here.

No morphologies are currently available for DCN neurons, so they are represented only by their soma.
Densities were estimated from `Blue Brain Cell Atlas <https://portal.bluebrain.epfl.ch/resources/models/cell-atlas/>`_ (version 2018 [#ero_2018]_), considering the ratio :math:`\frac{n_{GrC}}{n_{DCN}} = \frac{33 \times 10^6}{230 \times 10^3} ≈ 143` between the total number of granule cells and the total number of neurons in the cerebellar nuclei. From this ratio, considering the total amout of GrC placed in the ``canonical circuit``, it is possible to estimate the number of DCN to be placed. Literature data reported that DCNp are around the 57% of the total number of neurons in the cerebellar nuclei, while DCNi around the 32% [#baumel_2009]_ [#batini_1992]_. Taking into account these percentages and dividing by the volume of the DCN layer (set to :math:`200 \times 200 \times 300` µm), the values reported in the following table were obtained.

.. csv-table::
   :header-rows: 1
   :delim: ;

   Cell name;Type;Radius (:math:`µm`);Density (:math:`µm^{-3}`)
   DCNp ; Exc.; 9.5 [#baumel_2009]_; 9.92
   DCNi ; Inh.; 7.0 [#baumel_2009]_; 5.58

Placement
+++++++++
DCN are assumed to be uniformly distributed in their own layer, hence the bsb ``RandomPlacement`` strategy is chosen to place them.

Connectivity
++++++++++++

.. csv-table::
   :header-rows: 1
   :delim: ;

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; Specifics; References
   19; PC; axon; DCNp; / ; :ref:`fix_out`;``outdegree`` =45; Geminiani et al. (2024) [#geminiani_2024]_
   20; PC; axon; DCNi; / ; :ref:`fix_out`;``outdegree`` =12; Geminiani et al. (2024) [#geminiani_2024]_
   21; mf; / ; DCNp ; / ; :ref:`fix_in`; ``indegree`` =48; Geminiani et al. (2024) [#geminiani_2024]_


NEST simulation
^^^^^^^^^^^^^^^

Neuron parameters
+++++++++++++++++
DCN populations were represented as a EGLIF point neuron models (see :doc:`NEST section <nest>`). Parameters sets for both DCNp and DCNi are taken from Geminiani et al (2019) [#geminiani_2019]_.
The default LIF parameters are the following:

.. csv-table:: LIF neuron parameters for DCN
   :header-rows: 1
   :delim: ;

   Cell name;:math:`C_m\ (pF)`;:math:`\tau_m\ (ms)`;:math:`E_L\ (mV)`;:math:`t_{ref}\ (ms)`;:math:`V_{reset}\ (mV)`;:math:`V_{th}\ (mV)`
   DCNp; 142 (142 :math:`\pm` 31); 33 (33 :math:`\pm` 18); -45 (-45 :math:`\pm` 13); 1.5 (1.5 :math:`\pm` 0.2); -55 (-55); -36 (-36 :math:`\pm` 7)
   DCNi; 56 (56 :math:`\pm` 26); 56 (56 :math:`\pm` 30); -40 (-40 :math:`\pm` 13); 3.02 (3.02 :math:`\pm` 0.3); -55 (-55); -39 (-39 :math:`\pm` 8)

Then, the following parameters are optimized according to the method described in Geminiani et al. (2019) [#geminiani_2019]_ :

.. csv-table:: EGLIF neuron parameters for DCN
   :header-rows: 1
   :delim: ;

    Cell name;:math:`k_{adap}\ (nS \cdot ms^{-1})`;:math:`k_1\ (ms^{-1})`;:math:`k_2\ (ms^{-1})`;:math:`A_1\ (pA)`;:math:`A_2\ (pA)`;:math:`I_e\ (pA)`
    DCNp; 0.408; 0.697; 0.047; 13.857; 3.477; 75.385
    DCNi; 0.079; 0.041; 0.044; 176.358; 176.358; 2.384

It is not clear how the spiking parameters are obtained in the Geminiani et al. (2019) paper [#geminiani_2019]_. The values were extracted from a BSB configuration provided by the authors.

The postsynaptic receptors are defined as listed in Table 2 of Geminiani et al. (2019b) [#geminiani_2019b]_:

.. _dcn-table-receptor:
.. csv-table:: DCN Postsynaptic receptor parameters
   :header-rows: 1
   :delim: ;

   Cell name; Receptor id; :math:`E_{rev,i}\ (mV)`; :math:`\tau_{syn,i}\ (ms)`; Type
   DCNp; 1; 0; 1.0; exc.
   DCNp; 2; -80; 0.7; inh.
   DCNi; 1; 0; 3.64; exc.
   DCNi; 2; -80; 1.14; inh.

Synapse parameters
++++++++++++++++++
DCN connections are represented as ``static synapses`` (see :doc:`NEST section <nest>`). The receptor ids correspond to the postsynaptic receptors used for the connections (see table :ref:`dcn-table-receptor`).
It is still unclear from the references how these parameters were optimized.

.. csv-table:: Presynaptic parameters for DCN connections
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    mf-DCNp; 0.05; 4.0; 1
    PC-DCNp; 0.4; 4.0; 2
    PC-DCNi; 0.12; 4.0; 2

Simulation paradigms
++++++++++++++++++++

The `dcn_nest.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn_nest.yaml>`_ are
including all the simulation paradigms described in the :doc:`NEST section <nest>`) but include the DCN cells in the
circuit.

Basal activity
##############
For this simulation paradigm, the mean firing rates and mean ISI obtained for each neuron population are as
follows (expressed in mean :math:`\pm` standard deviation):

.. csv-table:: Results of the canonical circuit with DCN in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.1 \pm 0.93`; :math:`240 \pm 61`
    Granule cell; :math:`3.7 \pm 3.2`; :math:`430 \pm 410`
    Golgi cell;:math:`13 \pm 4.8`; :math:`92 \pm 46`
    Purkinje cell;:math:`49 \pm 3.7`; :math:`21 \pm 1.7`
    Basket cell;:math:`33 \pm 11`; :math:`34 \pm 13`
    Stellate cell;:math:`38 \pm 23`; :math:`42 \pm 32`
    DCNp; :math:`21 \pm 0.9`; :math:`47 \pm 2`
    DCNi; :math:`4.4 \pm 3.2`; :math:`130 \pm 15`

Mossy fiber stimulus
####################

For this simulation paradigm, **during the stimulus**, the mean firing rates and mean ISI obtained for each
neuron population are as follows (expressed in mean :math:`\pm` standard deviation):

.. csv-table:: Results of the canonical circuit with DCN during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`44 \pm 7`; :math:`7.0 \pm 3.7`
    Granule cell; :math:`24 \pm 48`; :math:`7.9 \pm 5.6`
    Golgi cell;:math:`48 \pm 38`; :math:`10.0 \pm 7.8`
    Purkinje cell;:math:`79 \pm 24`; :math:`12.0 \pm 4.4`
    Basket cell;:math:`110 \pm 83`; :math:`7.7 \pm 5.4`
    Stellate cell;:math:`150 \pm 100`; :math:`6.0 \pm 5.4`
    DCNp; :math:`23 \pm 7`; :math:`44.0 \pm 2.5`
    DCNi; :math:`1.2 \pm 4.8`; not enough spikes per neuron



References
^^^^^^^^^^

.. [#dangelo_2018] D'Angelo, Egidio. "Physiology of the cerebellum." Handbook of clinical neurology 154 (2018): 85-108. https://doi.org/10.1016/B978-0-444-63956-1.00006-0
.. [#uusisaari_2008] Uusisaari, M., and T. Knöpfel. "GABAergic synaptic communication in the GABAergic and non-GABAergic cells in the deep cerebellar nuclei." Neuroscience 156.3 (2008): 537-549. https://doi.org/10.1016/j.neuroscience.2008.07.060
.. [#geminiani_2019b] Geminiani, A., Pedrocchi, A., D’Angelo, E., & Casellato, C. (2019). Response
   dynamics in an olivocerebellar spiking neural network with non-linear neuron properties.
   Frontiers in computational neuroscience, 13, 68.
   https://doi.org/10.3389/fncom.2019.00068
.. [#baumel_2009] Baumel, Yuval, Gilad A. Jacobson, and Dana Cohen. "Implications of functional anatomy on information processing in the deep cerebellar nuclei." Frontiers in cellular neuroscience 3 (2009): 795. https://doi.org/10.3389/neuro.03.014.2009
.. [#batini_1992] Batini, Cesira, et al. "Cerebellar nuclei and the nucleocortical projections in the rat: retrograde tracing coupled to GABA and glutamate immunohistochemistry." Journal of Comparative Neurology 315.1 (1992): 74-84.  https://doi.org/10.1002/cne.903150106
.. [#ero_2018] Erö, Csaba, et al. "A cell atlas for the mouse brain." Frontiers in neuroinformatics 12 (2018): 84. https://doi.org/10.3389/fninf.2018.00084
.. [#geminiani_2024] Geminiani, Alice, et al. "Mesoscale simulations predict the role of synergistic cerebellar plasticity during classical eyeblink conditioning." PLOS Computational Biology 20.4 (2024): e1011277. https://doi.org/10.1371/journal.pcbi.1011277
.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035