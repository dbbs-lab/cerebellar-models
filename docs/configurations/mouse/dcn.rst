Deep Cerebellar Nuclei
~~~~~~~~~~~~~~~~~~~~~~

Deep cerebellar nuclei (DCN) are the primary output structures of the cerebellar cortex [#dangelo_2018]_.
They receive inhibitory input from Purkinje cells and excitatory input from mossy fibers and climbing fibers (from IO).
DCN are composed by three distinct nuclei: dentate nucleus, fastigial nucleus and interposed nucleus.
The default configuration with DCN is implemented in `dcn.yaml <https://github.com/dbbs-lab/cerebellum/blob/feature/dcn-io/configurations/mouse/dcn-io/dcn.yaml>`_.


Configuration
^^^^^^^^^^^^^

Cell types
++++++++++
For DCN, two types of neurons are considered [#uusisaari_2008]_ [#geminiani_2019b]_:

* **DCNp**: they are the excitatory neurons projecting outside the cerebellum to various brain regions, including the thalamus, the red nucleus, the vestibular nuclei, and the reticular formation;
* **DCNi**: they are GABAergic interneurons which send inhibitory feedback to IO.

No morphologies are currently available for DCN neurons, so they are represented only by their soma.
Densities were estimated from `Blue Brain Cell Atlas <https://portal.bluebrain.epfl.ch/resources/models/cell-atlas/>`_ (version 2018 [#ero_2018]_), considering the ratio :math:`\frac{n_{GrC}}{n_{DCN}} = \frac{33 \times 10^6}{230 \times 10^3} ≈ 143` between the total number of granule cells and the total number of neurons in the cerebellar nuclei. From this ratio, considering the total amout of GrC placed in the canonical reconstruction, it is possible to estimate the number of DCN to be placed. Literature data reported that DCNp are around the 57% of the total number of neurons in the cerebellar nuclei, while DCNi around the 32% [#baumel_2009]_ [#batini_1992]_. Taking into account these percentages and dividing by the volume of the DCN layer (set to :math:`200 \times 200 \times 300` µm), the values reported in the following table were obtained.

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

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; References
   19; PC; axon; DCNp; \ ; cerebellum.connectome.purkinje_dcn.FixedOutdegree; \ ;




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