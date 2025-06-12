Placement strategies
--------------------

:class:`LabelMicrozones<.placement.microzones.LabelMicrozones>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BSB :class:`AfterPlacementHook<bsb.postprocessing.AfterPlacementHook>` to subdivide
cell populations into labelled subpopulations of same cell counts based on their
position along a provided axis.
The number of labels defines the number of subpopulations.

* ``cell_types``: List of CellType of the population to subdivide, used as references;

* ``axis``: Index of the axis used to split the population, by default 0 (x axis);

* ``labels``: List of labels to assign to each population respectively, defines also the number of subpopulations;
