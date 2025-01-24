Mouse cerebellar cortex configuration
=====================================

Introduction
------------

The ``mouse`` folder contains the configurations for the reconstruction and simulation of the mouse
cerebellar cortex with BSB v4.

These reconstructions are based on the iterative work of many researchers distributed in many
papers. The role of this file is to make explicit the origin of each value and strategy extracted
from the literature and integrated into these configurations. All the configurations present in
this folder are based on the
`mouse_cerebellar_cortex.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/mouse_cerebellar_cortex.yaml>`_
file. It corresponds to the configuration file written for the reconstruction of the cerebellar
circuit presented in the De Schepper et al. (2022) [#de_schepper_2022]_ paper. This circuit
configuration will be later referred to as the ``canonical circuit``.

We will follow the structure of the BSB configuration files to present each of their sections and
the data they leverage.

.. include:: ../../getting-started/biological-context.rst
   :end-before: figure

.. Do not include the Figure to prevent double referencing

.. include:: ../../getting-started/biological-context.rst
   :start-after: end-figure
   :end-before: end-bio-context

In the mouse ``canonical circuit`` case, the UBC, DCN, IO and related connections are not included.

Circuit configuration
---------------------

Coordinate framework
~~~~~~~~~~~~~~~~~~~~

By convention, the circuit is oriented so that its layers are stacked vertically with the granular
layer at the bottom and the molecular layer at the top. We derived a coordinate framework
``(x,y,z)`` from this layout, based on the right-hand orientation convention. Its origin is set at
the bottom of the circuit, the z axis pointing to the top of the molecular layer. The ``(y-z)``
plane corresponds to the para-sagittal sections that are co-planar with the Purkinje dendritic trees
and normal to the granule cells parallel fibers. Finally, unit of distance in the configurations are
expressed in micrometers ``µm``. Note that the morphologies provided are oriented by default to
match this convention.

Network dimensions
~~~~~~~~~~~~~~~~~~

The ``canonical circuit`` is built in a cubic volume of :math:`300 \times 200 \times 295` µm in the
``(x,y,z)`` convention (see ``network``, ``regions`` and ``partitions`` in the configuration). The
thickness of each of its layer has been determined according to literature findings and to match the
size and shape of the available morphologies:

- The Purkinje layer corresponds to a one cell thick sheet of Purkinje cells. The Purkinje cell soma
  diameters determine therefore the thickness of this layer. According to
  Hendelman & Aggerwal (1980) [#hendelman_1980]_, Purkinje cell’s soma diameters have been
  estimated to less than 20µm in mice. We chose here ``15µm``.

- The molecular layer total thickness has been calculated to fit the size of the dendritic
  arborization of the Purkinje cell’s morphology as ``150µm``. The molecular layer is itself divided
  in the ``canonical circuit`` into two sub-layers based on their neuronal composition. Here, the
  bottom part of the molecular layer (hence the part stacked right on top of the Purkinje layer)
  contains the Basket cells (``b_molecular_layer`` in the configuration); while the top part hold
  the Stellate cells (``s_molecular_layer``). In fact, according to literature data
  (J. Kim & Augustine, 2021 [#kim_2021]_; Sultan & Bower, 1998 [#sultan_1998]_), SCs are more likely
  located in the outer two-third of the Molecular layer. While this distribution of cells is closer
  to gradient in real mice, we assumed a clear separation between the populations. The
  ``basket layer`` is therefore ``50µm`` thick while the ``stellate layer`` is ``100µm`` thick.

- The granular layer’s thickness has been similarly fitted to match the size of the Golgi cell basal
  dendritic tree, here ``130µm``. Note that the size of the granule cell ascending axons have been
  set to this constraint.

Cellular composition
~~~~~~~~~~~~~~~~~~~~

The cellular composition of the circuit is determined in the ``cell_types`` section of the
configuration file. Each cell type is linked here to a partition in the circuit, and a morphology is
assigned. Additionally, we introduce here two components used to innervate the circuit: the mossy
fibers originating from various other brain regions, and the glomeruli that form at their
terminals. These components are used as building entities to relay stimuli from other regions into
the circuit and have therefore no morphology attached.

We will describe here the spatial parameters used in ``canonical circuit``:

.. csv-table::
   :header-rows: 1
   :delim: ;

   Layer;Cell name;Type;Radius (:math:`µm`);Density (:math:`µm^{-3}`);References
   Granular layer;Glomerulus (glom);Exc.;1.5;0.0003;Solinas et al. (2010) [#solinas_2010]_
   Granular layer;Mossy fibers (mf);Exc.;/;count relative to glom. ratio=0.05;Billings et al. (2014) [#billings_2014]_
   Granular layer;Granule Cell (GrC);Exc.;2.5;0.0039;Casali et al. (2019) [#casali_2019]_
   Granular layer;Golgi Cell (GoC);Inh.;4.0;0.000009;Casali et al. (2019) [#casali_2019]_
   Purkinje layer;Purkinje cell (PC);Inh.;7.5;planar density: 0.001166;Keller et al. (2018) [#keller_2018]_
   Molecular layer;Basket cell (BC);Inh.;6.;0.00005;Casali et al. (2019) [#casali_2019]_
   Molecular layer;Stellate cell (SC);Inh.;4.;0.00005;Casali et al. (2019) [#casali_2019]_

.. warning::
   Note that most literature data in this table comes from rat data.

The density of glom have been calculated in Solinas et al. (2010) [#solinas_2010]_ based on the
glomerulus to granule convergence and divergence ratios (derived from values in Korbo et al., 1993
[#korbo_1993]_ and Jakab and Hámari, 1988 [#jakab_1988]_).

For PC cells, the planar density was calculated to obtain approximately 70 PC in our Purkinje Layer,
which results in a planar density of 1166 :math:`PCs/mm^{2}`.
This density is consistent with the data from Keller et al. [#keller_2018]_

The densities of GrC, GoC, BC and SC are reported in Table 1 of Casali et al. (2019) [#casali_2019]_.
The authors cite Korbo et al. (1993) [#korbo_1993]_ for the values in this table, however, no
equivalent was found in the cited paper. These values might have been optimized to improve
simulation results. They could also have been obtained through geometric constrained placement
to minimize the overlap of somas.

.. include:: ../morphologies.rst

Placement
~~~~~~~~~

Except for Purkinje cells (PC), every entity is supposed to be uniformly distributed in their own
layer.The bsb ``RandomPlacement`` strategy is chosen here to place them. In short, this strategy
chose a random position for each entity within their sub-partition. Note that this does not take
into account any potential overlapping of cells’ soma unlike the ``ParticlePlacement``.
However, comparative analysis conducted in ``[CITATION]`` have shown that the latter strategy have a
limited impact on connectivity and simulation results, while the computational cost of checking soma
overlapping is not negligible.

PC are placed in arrays, ``130 μm`` apart from each other along the
para-sagittal plane ``(xz)`` to guarantee that their dendritic
arborizations do not overlap. Furthermore, each row of PC somas is
shifted with respect to its predecessor to form a ``80`` degree angle on
the ``(xy)`` plane.

Connectivity
~~~~~~~~~~~~
The following table list all the connections present in the model. The connection id of the first column corresponds to
the numbers reported in :numref:`fig-network`.

.. _table-connectivity:
.. csv-table:: Connectivity rules of the cerebellar cortex circuit
   :header-rows: 1
   :delim: ;

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; References
   1; mf ; /; glom ; /; :ref:`mossy_glom`; Sultan (2001) [#sultan_2001]_
   2; glom ; /; GrC; dendrites; :ref:`glom_grc`; Houston et al. (2017) [#houston_2017]_
   3; glom ; /; GoC; basal dendrites; :ref:`glom_goc`; Kanichay and Silver (2008) [#kanichay_2008]_
   4; GoC; axon ; GrC; same as through glom; :ref:`goc_glom`; Barmack and Yakhnitsa (2008) [#barmack_2008]_
   5; GoC; axon ; GoC; basal dendrites; :ref:`voxel_int` ; Hull and Regehr (2012) [#hull_2012]_
   6; GrC; ascending axon ; GoC; basal dendrites; :ref:`voxel_int` ; Cesana et al. (2013) [#cesana_2013]_
   7; GrC; parallel fiber ; GoC; apical dendrites ; :ref:`voxel_int` ; Kanichay and Silver (2008) [#kanichay_2008]_
   8; GrC; ascending axon ; PC ; ascending axon targets ; :ref:`voxel_int` ; Wang and Huang (2006) [#wang_2006]_
   9; GrC; parallel fiber ; PC ; parallel fiber targets ; :ref:`voxel_int` ; Wang and Huang (2006) [#wang_2006]_
   10; GrC; parallel fiber ; BC ; dendrites; :ref:`voxel_int` ; Jörntell et al. (2010) [#jorntell_2010]_
   11; GrC; parallel fiber ; SC ; dendrites; :ref:`voxel_int` ; Jörntell et al. (2010) [#jorntell_2010]_
   12; BC ; axon ; PC ; soma ; :ref:`voxel_int` ; Jörntell et al. (2010) [#jorntell_2010]_
   13; SC ; axon ; PC ; stellate cell targets; :ref:`voxel_int` ; Jörntell et al. (2010) [#jorntell_2010]_
   14; BC ; axon ; BC ; dendrites; :ref:`voxel_int` ; Ito (2013) [#ito_2013]_
   15; SC ; axon ; SC ; dendrites; :ref:`voxel_int` ; Ito (2013) [#ito_2013]_

Parameters explanation:
^^^^^^^^^^^^^^^^^^^^^^^

We currently have no morphologies for the mf and glom, which makes it impossible to use fiber or
voxel intersection techniques to implement their related connection rules. We therefore simplified
the geometry of the neurites involved in these connections.

Sultan provides ranges of distance between gloms and their respective mf in their paper [#sultan_2001]_:
:math:`57.6 \pm 60 \times 19.6 \pm 18.8` µm along respectively the x and y axes in our coordinate
system. In our model, we used the rounded mean values as the maximum distance between mf and their
gloms.

GrC of the adult mouse cerebellar cortex has :math:`3.9 \pm 0.1` dendrites that spreads for ~40μm in
each direction, as reported in Houston et al. (2017) [#houston_2017]_ (see Figure 2G).
In our model, we therefore assumed that each GrC has ``4`` dendrites of ``40μm``, to
match also the number of branches of the respective morphology. These values are used in the glom to
GrC connectivity rule. The convergence value of this connection pair is set here to the number of
dendrites.

GoC basolateral arborizations spread across 100μm in P25 rat according to Kanichay and Silver (2008)
[#kanichay_2008]_. This has been simplified to a sphere of ``50μm`` radius surrounding their soma
for the GoC to glom connectivity.

Barmack and Yakhnitsa (2008) [#barmack_2008]_ reported that the mean mediolateral extent of the GoC
axon is :math:`180 \pm 40` μm, and that it spreads along the parasagittal plane. In our connection
from GoC to GrC (through glom), we used a ``150μm`` sphere surrounding the GoC soma to find
potential glom targets. The maximum number of glom target (divergence) for each GoC was set to ``40``
in Solinas et al. (2010) [#solinas_2010]_. However, the rationale behind this particular value is
unclear but probably to balance the granular layer excitation and inhibition.

For the rest of the connection rules, we leveraged each neuron morphologies to detect appositions of
their neurites. Fiber intersection methods require a lot of computational power. For this reason, we used
BSB :ref:`voxel_int` strategy, as it simplifies this detection representing morphologies using
voxels.

The ``affinity`` and ``distributions`` of ``contact`` points parameters of these connections were
tuned to match connectivity divergence and convergence values from De Schepper et al. (2022)
[#de_schepper_2022]_.

We introduce a correction for the SC-PC connectivity which consists in a normal distribution of synapse points
per SC-PC connections (parameters ``loc`` = 10, ``scale`` = 0.3).
The idea was to have ~50 afferent synapses per PC to match the range of inhibition efficiency described in
Rizza et al. 2021 [#rizza_2021]_

Extensions to the canonical model
---------------------------------

See the :doc:`Mouse extensions section <extensions>`.

References
----------

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

.. [#keller_2018] Keller, D. (2018). Cell Densities in the Mouse Brain: A Systematic Review.
   Frontiers in Neuroanatomy, 12. https://doi.org/10.3389/fnana.2018.00083

.. [#rizza_2021] Rizza, M. F., Locatelli, F., Masoli, S., Sánchez-Ponce, D., Muñoz, A., Prestori, F.,
   & D’Angelo, E. (2021). Stellate cell computational modeling predicts signal filtering in the molecular
   layer circuit of cerebellum. Scientific Reports, 11(1), 3873.
   https://doi.org/10.1038/s41598-021-83209-w
