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
`io.yaml <https://github.com/dbbs-lab/cerebellar-models/blob/master/configurations/mouse/cell_types/io.yaml>`_.

Configuration
^^^^^^^^^^^^^
In `io.yaml <https://github.com/dbbs-lab/cerebellar-models/blob/master/configurations/mouse/cell_types/io.yaml>`_ ,
a new region called ``inferior_olivary`` was added to the ``canonical circuit + DCN`` model
(see :doc:`DCN section <dcn>` for more microcircuit info).
This region contains only one ``Layer`` Partition: ``io layer``. ``io layer`` has a thickness of :math:`100  \mu m` .
Additionally, to ensure that ``inferior_olivary`` are placed under the ``cerebellar_nuclei``, its ``origin``
was set to ``[0, 0, -300]``.

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
   24; IO; / ; DCNp ; / ; :ref:`all_to_all`; ``affinity`` = 0.5 ; Geminiani et al. (2024) [#geminiani_2024]_
   25; IO; /; DCNi; / ; :ref:`all_to_all`; ``affinity`` = 0.5 ; Geminiani et al. (2024) [#geminiani_2024]_
   26; DCNi; / ; IO ; / ; :ref:`all_to_all`; ``affinity`` = 0.5 ; Geminiani et al. (2024) [#geminiani_2024]_

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
