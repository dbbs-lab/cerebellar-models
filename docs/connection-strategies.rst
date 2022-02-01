
:class:`ConnectomeGlomerulusGranule <.connectivity.ConnectomeGlomerulusGranule>`
================================================================================

Inherits from TouchingConvergenceDivergence. No additional configuration.
Uses the dendrite length configured in the granule cell morphology.

:class:`ConnectomeGlomerulusGolgi <.connectivity.ConnectomeGlomerulusGolgi>`
============================================================================

Inherits from TouchingConvergenceDivergence. No additional configuration.
Uses the dendrite radius configured in the Golgi cell morphology.

:class:`ConnectomeGolgiGlomerulus <.connectivity.ConnectomeGolgiGlomerulus>`
============================================================================

Inherits from TouchingConvergenceDivergence. No additional configuration.
Uses the ``axon_x``, ``axon_y``, ``axon_z`` from the Golgi cell morphology
to intersect a parallelopipid Golgi axonal region with the glomeruli.

:class:`ConnectomeGranuleGolgi <.connectivity.ConnectomeGranuleGolgi>`
======================================================================

Creates 2 connectivity sets by default *ascending_axon_to_golgi* and
*parallel_fiber_to_golgi* but these can be overwritten by providing ``tag_aa``
and/or ``tag_pf`` respectively.

Calculates the distance in the XZ plane between granule cells and Golgi cells and
uses the Golgi cell morphology's dendrite radius to decide on the intersection.

Also creates an ascending axon height for each granule cell.

* ``aa_convergence``: Preferred amount of ascending axon synapses on 1 Golgi cell.
* ``pf_convergence``: Preferred amount of parallel fiber synapses on 1 Golgi cell.

:class:`ConnectomeGolgiGranule <.connectivity.ConnectomeGolgiGranule>`
======================================================================

No configuration, it connects each Golgi to each granule cell that it shares a
connected glomerules with.

:class:`ConnectomeAscAxonPurkinje <.connectivity.ConnectomeAscAxonPurkinje>`
============================================================================

Intersects the rectangular extension of the Purkinje dendritic tree with the granule
cells in the XZ plane, uses the Purkinje cell's placement attributes ``extension_x``
and ``extension_z``.

* ``extension_x``: Extension of the dendritic tree in the X plane
* ``extension_z``: Extension of the dendritic tree in the Z plane

:class:`ConnectomePFPurkinje <.connectivity.ConnectomePFPurkinje>`
==================================================================

No configuration. Uses the Purkinje cell's placement attribute ``extension_x``.
Intersects Purkinje cell dendritic tree extension along the x axis with the x position
of the granule cells, as the length of a parallel fiber far exceeds the simulation
volume.
