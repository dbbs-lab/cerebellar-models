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
This region was defined of type ``group`` because in it there is only the ``io layer``. ``io layer`` has a thickness of ``200 µm`` .
Additionally, in order to ensure that ``inferior_olivary`` are placed under the ``cerebellar_nuclei``, the ``origin``
of the ``dcn_layer`` was set to ``[0,0,200]`` and the ``origin`` of the ``granular_layer`` was updated to ``[0,0,400]``.

Cell types
++++++++++
No morphologies are currently available for inferior olivary neurons, so they are represented only by their soma.





References
^^^^^^^^^^

.. [#de_gruijl_2012] De Gruijl, J. R., Bazzigaluppi, P., de Jeu, M. T., & De Zeeuw, C. I. (2012). Climbing fiber burst size and olivary sub-threshold oscillations in a network setting. PLoS computational biology, 8(12), e1002814.