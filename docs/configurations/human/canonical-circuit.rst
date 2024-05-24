Human cerebellar cortex configuration
=====================================

Introduction
------------

The ``human`` folder contains the configurations for the reconstruction and simulation of the mouse
cerebellar cortex with BSB v4.

These reconstructions are based on the iterative work of many researchers distributed in many
papers. The role of this file is to make explicit the origin of each value and strategy extracted
from the literature and integrated into these configurations. All the configurations present in
this folder are based on the
`human_cerebellar_cortex.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/human/human_cerebellar_cortex.yaml>`_
file. This circuit configuration will be later referred to as the ``canonical circuit``.

We will follow the structure of the BSB configuration files to present each of their sections and
the data they leverage.

.. include:: ../../getting-started/biological-context.rst
   :end-before: figure

.. Do not include the Figure to prevent double referencing

.. include:: ../../getting-started/biological-context.rst
   :start-after: end-figure
   :end-before: end-bio-context

In the human ``canonical circuit`` case, the UBC, DCN, IO and related connections are not included.

Circuit configuration
---------------------

.. include:: ../mouse/canonical-circuit.rst
   :start-after: begin_coord_framework
   :end-before: end_coord_framework

Network dimensions
~~~~~~~~~~~~~~~~~~

The ``canonical circuit`` is built in a cubic volume of :math:`1000 \times 1000 \times 900` µm in the
``(x,y,z)`` convention (see ``network``, ``regions`` and ``partitions`` in the configuration). The
thickness of each of its layer has been determined according to literature findings and to match the
size and shape of the available morphologies:

- The Purkinje layer corresponds to a one cell thick sheet of Purkinje cells. The Purkinje cell soma
  diameters determine therefore the thickness of this layer. According to
  Buchholz et al. (2020) [#Buchholz_2020]_, Purkinje cell’s soma diameters have been
  estimated between ``15`` and ``30`` in humans. We chose here ``30µm``.

- The molecular layer total thickness in the ``canonical circuit`` is ``330µm`` (Lobule VI molecular layer 
  thickness from Zheng et al. (2023) [#Zheng_2023]_). The molecular layer is itself divided
  in the ``canonical circuit`` into two sub-layers based on their neuronal composition. Here, the
  bottom part of the molecular layer (hence the part stacked right on top of the Purkinje layer)
  contains the Basket cells (``b_molecular_layer`` in the configuration); while the top part hold
  the Stellate cells (``s_molecular_layer``). In fact, according to literature data
  (Whitney et al. (2009) [#Whitney_2009]_; J. Kim & Augustine, 2021 [#kim_2021]_; Sultan & Bower, 1998 [#sultan_1998]_), 
  SCs are more likely located in the outer two-third of the Molecular layer. While this distribution of cells is closer
  to gradient in real humans, we assumed a clear separation between the populations. The
  ``basket layer`` is therefore ``110µm`` thick while the ``stellate layer`` is ``220µm`` thick (Whitney et al. (2009) [#Whitney_2009]_).

- The granular layer’s thickness has in in the ``canonical circuit`` is ``610µm`` (Lobule VI, Zheng et al. (2023) [#Zheng_2023]_).

Cellular composition
~~~~~~~~~~~~~~~~~~~~

.. include:: ../mouse/canonical-circuit.rst
   :start-after: begin_cellular_composition_intro
   :end-before: end_cellular_composition_intro

.. csv-table::
   :header-rows: 1
   :delim: ;

   Layer;Cell name;Type;Radius (:math:`µm`);Density (:math:`µm^{-3}`);References
   Granular layer;Glomerulus (glom);Exc.;3.0;0.00015;Solinas et al. (2010) [#solinas_2010]_
   Granular layer;Mossy fibers (mf);Exc.;/;count relative to glom. ratio=0.05;Billings et al. (2014) [#billings_2014]_
   Granular layer;Granule Cell (GrC);Exc.;2.9;0.001704;Andersen et al. (2012) [#Andersen_2012]_
   Granular layer;Golgi Cell (GoC);Inh.;8.0;0.00000142;Solinas et al. (2010) [#solinas_2010]_
   Purkinje layer;Purkinje cell (PC);Inh.;15;planar density: 0.00028;Andersen et al. (2012) [#Andersen_2012]_
   Molecular layer;Basket cell (BC);Inh.;12.;0.000023751;Whitney et al. (2009) [#Whitney_2009]_
   Molecular layer;Stellate cell (SC);Inh.;8.;0.000014794;Whitney et al. (2009) [#Whitney_2009]_

The density of glom have been obtained scaling down the corresponding mouse value by a factor 1/2 according to the 
correlation between the number of cells per :math:`mm^3` and the weight of the cerebellum in Lange 1975 [#Lange_1975]_. 
The gloms density for mouse is taken from Solinas et al. (2010) [#solinas_2010]_ based on the
glomerulus to granule convergence and divergence ratios (derived from values in Korbo et al., 1993
[#korbo_1993]_ and Jakab and Hámari, 1988 [#jakab_1988]_).

.. include:: ../morphologies.rst

Placement
~~~~~~~~~

.. include:: ../mouse/canonical-circuit.rst
   :start-after: begin_pc_placement
   :end-before: end_pc_placement

PC are placed in arrays, ``350 μm`` apart from each other along the
para-sagittal plane ``(xz)`` to guarantee that their dendritic
arborizations do not overlap. Furthermore, each row of PC somas is
shifted with respect to its predecessor to form a ``70`` degree angle on
the ``(xy)`` plane.

Connectivity
~~~~~~~~~~~~

.. include:: ../mouse/canonical-circuit.rst
   :start-after: begin_connectivity_common
   :end-before: end_connectivity_common

GoC basolateral arborizations spread across ``100μm`` in P25 rat according to Kanichay and Silver (2008)
[#kanichay_2008]_. This has been scaled up to ``200μm`` for the human case. Similarly to the mouse canonical circuit, 
for the GoC to glom connectivity, the GoC basolateral arborizations has been simplified to a sphere centered at the GoC soma,
in this case of radius ``100μm``.

Barmack and Yakhnitsa (2008) [#barmack_2008]_ reported that the mean mediolateral extent of the GoC
axon is :math:`180 \pm 40` μm, and that it spreads along the parasagittal plane. In our connection
from GoC to GrC (through glom), we used a ``200μm`` sphere surrounding the GoC soma to find
potential glom targets. The maximum number of glom target (divergence) for each GoC was set to ``40``
in Solinas et al. (2010) [#solinas_2010]_. However, the rationale behind this particular value is
unclear but probably to balance the granular layer excitation and inhibition.

For the rest of the connection rules, we leveraged each neuron morphologies to detect appositions of
their neurites. Fiber intersection methods require a lot of computational power. For this reason, we used
BSB :ref:`voxel_int` strategy, as it simplifies this detection representing morphologies using
voxels.

The ``affinity`` and ``distributions`` of ``contact`` points parameters of these connections were
tuned to obtain the same surface synaptic densities as in De Schepper et al. (2022) [#de_schepper_2022]_.


References
----------

.. [#Zheng_2023] Zheng, J. et al. (2023). 
   Three-Dimensional Digital Reconstruction of the Cerebellar Cortex: 
   Lobule Thickness, Surface Area Measurements, and Layer Architecture. 
   The Cerebellum, 22(2), 249-260.

.. [#Buchholz_2020] Buchholz D.E., Carroll T.S., Kocabas A., Zhu X., Behesti H., Faust P.L., Stalbow L., Fang Y. & Hatten M.E. (2020)
   Novel genetic features of human and mouse Purkinje cell differentiation defined by comparative transcriptomics. 
   Proc Natl Acad Sci USA. 117(26):15085-15095.

.. [#Whitney_2009] Whitney, E.R., Kemper, T.L., Rosene, D.L., Bauman, M.L. and Blatt, G.J. (2009).
   Density of cerebellar basket and stellate cells in autism: Evidence for a late developmental loss of Purkinje cells.
   J. Neurosci. Res., 87: 2245-2254.

.. [#Lange_1975] Lange W. (1975). Cell number and cell density in the cerebellar cortex of man and some other mammals. 
   Cell Tissue Res. 1975;157(1):115-24.

.. [#Andersen_2012] Andersen K., Andersen B.B. & Pakkenberg B. (2012).
   Stereological quantification of the cerebellum in patients with Alzheimer's disease. 
   Neurobiol Aging. 33(1):197.e11-20.

.. [#de_schepper_2022] De Schepper, R., Geminiani, A., Masoli, S., Rizza, M. F., Antonietti, A., &
   Casellato, C. (2022). Model simulations unveil the structure-function-dynamics relationship of
   the cerebellar cortical microcircuit. Communications Biology, 5(1), 1-19.
   https://doi.org/10.1038/s42003-022-04213-y

.. [#hendelman_1980] Hendelman, W. J., & Aggerwal, A. S. (1980). The Purkinje neuron: I. A Golgi
   study of its development in the mouse and in culture. Journal of Comparative Neurology, 193(4),
   1063–1079. https://doi.org/10.1002/cne.901930417

.. [#kim_2021] Kim, J., & Augustine, G. J. (2021). Molecular Layer Interneurons: Key Elements of
   Cerebellar Network Computation and Behavior. Neuroscience, 462, 22-35.
   https://doi.org/10.1016/j.neuroscience.2020.10.008

.. [#sultan_1998] Sultan, F., & Bower, J. M. (1998). Quantitative Golgi study of the rat cerebellar
   molecular layer interneurons using principal component analysis. Journal of Comparative
   Neurology, 393(3), 353-373. PMID: 9548555.
   https://doi.org/10.1002/(SICI)1096-9861(19980413)393:3<353::AID-CNE7>3.0.CO;2-0

.. [#solinas_2010] Solinas, S., Nieus, T., & D‘Angelo, E. (2010). A realistic large-scale model of
   the cerebellum granular layer predicts circuit spatio-temporal filtering properties. Frontiers in
   cellular neuroscience, 4, 903.  doi: 10.3389/fncel.2010.00012. PMID: 20508743; PMCID: PMC2876868.

.. [#billings_2014] Billings, G., Piasini, E., Lőrincz, A., Nusser, Z., & Silver, R. A. (2014).
   Network structure within the cerebellar input layer enables lossless sparse encoding. Neuron,
   83(4), 960-974. https://doi.org/10.1016/j.neuron.2014.07.020

.. [#casali_2019] Casali, S., Marenzi, E., Medini, C., Casellato, C., & D'Angelo, E. (2019).
   Reconstruction and simulation of a scaffold model of the cerebellar network. Frontiers in
   neuroinformatics, 13, 444802.https://doi.org/10.3389/fninf.2019.00037

.. [#korbo_1993] Korbo, L., Andersen, B. B., Ladefoged, O., & Møller, A. (1993). Total numbers of
   various cell types in rat cerebellar cortex estimated using an unbiased stereological method.
   Brain research, 609(1-2), 262-268. https://doi.org/10.1016/0006-8993(93)90881-M

.. [#jakab_1988] Jakab, R. L., & Hamori, J. (1988). Quantitative morphology and synaptology of
   cerebellar glomeruli in the rat. Anatomy and embryology, 179, 81-88.
   https://doi.org/10.1007/BF00305102

.. [#sultan_2001] Sultan, F. (2001). Distribution of mossy fibre rosettes in the cerebellum of cat
   and mice: Evidence for a parasagittal organization at the single fibre level. European Journal of
   Neuroscience, 13(11), 2123-2130. https://doi.org/10.1046/j.0953-816x.2001.01593.x

.. [#kanichay_2008] Kanichay, R. T., & Silver, R. A. (2008). Synaptic and cellular properties of the
   feedforward inhibitory circuit within the input layer of the cerebellar cortex. Journal of
   Neuroscience, 28(36), 8955-8967. https://doi.org/10.1523/JNEUROSCI.5469-07.2008

.. [#houston_2017] Houston, C. M., Diamanti, E., Diamantaki, M., Kutsarova, E., Cook, A., Sultan, F.,
   & Brickley, S. G. (2017). Exploring the significance of morphological diversity for cerebellar
   granule cell excitability. Scientific Reports, 7(1), 1-16. https://doi.org/10.1038/srep46147

.. [#barmack_2008] Barmack, N. H., & Yakhnitsa, V. (2008). Functions of interneurons in mouse
   cerebellum. Journal of Neuroscience, 28(5), 1140-1152.
   https://doi.org/10.1523/JNEUROSCI.3942-07.2008

.. [#hull_2012] Hull, C., Regehr, W. G. (2012). Identification of an inhibitory circuit that
   regulates cerebellar Golgi cell activity. Neuron, 73(1), 149-158.
   https://doi.org/10.1016/j.neuron.2011.10.030

.. [#cesana_2013] Cesana, E., Pietrajtis, K., Bidoret, C., Isope, P., D'Angelo, E., Dieudonné, S.,
   & Forti, L. (2013). Granule cell ascending axon excitatory synapses onto Golgi cells implement a
   potent feedback circuit in the cerebellar granular layer. Journal of Neuroscience, 33(30),
   12430-12446. https://doi.org/10.1523/JNEUROSCI.4897-11.2013

.. [#wang_2006] Wang, L., & Huang, R. H. (2006). Cerebellar granule cell: Ascending axon and
   parallel fiber. European Journal of Neuroscience, 23(7), 1731-1737.
   https://doi.org/10.1111/j.1460-9568.2006.04690.x

.. [#jorntell_2010] Jörntell, H., Bengtsson, F., Schonewille, M., & De Zeeuw, C. I. (2010).
   Cerebellar molecular layer interneurons–computational properties and roles in learning. Trends in
   neurosciences, 33(11), 524-532.https://doi.org/10.1016/j.tins.2010.08.004

.. [#ito_2013] Ito, M. (2013). Cerebellar Microcircuitry. Reference Module in Biomedical Sciences.
   https://doi.org/10.1016/B978-0-12-801238-3.04544-X
