Connection strategies
---------------------

.. _mossy_glom:

:class:`ConnectomeMossyGlomerulus <.connectome.to_glomerulus.ConnectomeMossyGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects one presynaptic cell (usually mossy fiber) within the
:math:`x\_length \times z\_length \times y\_size` μm box surrounding each postsynaptic cell
(usually glomerulus). Here ``y_size`` is the size of the partition where the postsynaptic cell are
(i.e. no limit). This selection is random and performed with a truncated exponential distribution.

* ``x_length``: Length of the box along the x axis surrounding the postsynaptic cell soma in which
  the presynaptic cell can be connected.

* ``y_length``: Length of the box along the y axis surrounding the postsynaptic cell soma in which
  the presynaptic cell can be connected.

.. note::
    Since the placement of mossy fibers and glomerulus is uniformed within their partition,
    the convergence and divergence ratios corresponds to the density ratios.

.. _glom_goc:

:class:`ConnectomeGlomerulusGolgi <.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects all presynaptic cells (usually glomeruli) within the sphere surrounding each
(postsynaptic cell soma (usually golgi). For each unique postsynaptic cell to connect, the tip of a
branch from the postsynaptic cell morphology is selected. This selection is random and performed
with a truncated exponential distribution based on the distance between the tip of each branch and
the glomerulus to connect.

* ``radius``: Radius of the sphere to filter the presynaptic cells within it.

For the glomerulus to golgi connectivity, the ``basal dendrites`` of the golgi cells are selected.

.. code-block:: yaml

      glomerulus_to_golgi:
        strategy: cerebellum.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi
        presynaptic:
          cell_types:
            - glomerulus
        postsynaptic:
          cell_types:
            - golgi_cell
          morphology_labels:
            - basal_dendrites
        radius: 50

.. note::
    As for the glomeruli and golgi cells, the presynaptic cell of this connection is supposed to
    have no morphology attached while the postsynaptic cell must have one.

.. _glom_grc:

:class:`ConnectomeGlomerulusGranule <.connectome.glomerulus_granule.ConnectomeGlomerulusGranule>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects ``convergence`` unique presynaptic cells (usually glomeruli) within the sphere
surrounding each postsynaptic cell soma (usually granule cell). Moreover, each presynaptic cell
should belong to a unique pre-presynaptic cluster (usually refer to unique mossy fiber cluster).
In other words, for a postsynaptic cell, each selected presynaptic cell should be connected
through a dependency strategy (usually :ref:`mossy_glom`) to a different pre-presynaptic cell.
The pre-presynaptic cluster, the presynaptic cell and the postsynaptic cell branch are all randomly
chosen. If not enough presynaptic cells could be found in the sphere surrounding the postsynaptic
cell soma, the closest presynaptic connected to the remaining cluster are selected to connect.

* ``radius``: Radius of the sphere to filter the presynaptic cells within it.

* ``convergence``: Number of presynaptic cell per postsynaptic cell.

* ``mf_glom_strat``: ConnectionStrategy that links the pre-presynaptic cell to the presynaptic cell.

* ``mf_cell_type``: CellType of the pre-presynaptic cell, used as a reference

For the glomerulus to granule connectivity, the ``dendrites`` of the granule cells are selected.

.. code-block:: yaml

      glomerulus_to_granule:
        strategy: cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule
        presynaptic:
          cell_types:
            - glomerulus
        postsynaptic:
          cell_types:
            - granule_cell
          morphology_labels:
            - dendrites
        mf_glom_strat: mossy_fibers_to_glomerulus
        mf_cell_type: mossy_fibers
        radius: 40
        convergence: 4

.. note::
    As for the glomeruli and granule cells, the presynaptic cell of this connection is supposed to
    have no morphology attached while the postsynaptic cell must have one.

.. _goc_glom:

:class:`ConnectomeGolgiGlomerulus <.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This connection strategy links a presynaptic cell (usually golgi) to all postsynaptic cells
(usually granule or ubc) connected to an intermediate cell (usually a glomerulus). The connections
between the intermediate cell and the postsynaptic cell need are defined in a separate connection
strategy (usually :ref:`glom_grc`).

The algorithm selects here the closest intermediate cell (maximum ``divergence``) that are within
a sphere surrounding each presynaptic soma. For each unique intermediate cell selected,
the tip of a branch from the presynaptic cell is randomly selected. All postsynaptic cells connected
to the selected intermediate cell through the dependent strategy are also connected to the selected
presynaptic cell. The target points of the postsynaptic cell (branch point) are copied from the
dependent strategy.

* ``divergence``: Number of postsynaptic cell per presynaptic cell.

* ``radius``: Radius of the sphere to filter the intermediate cells within it.

* ``glom_post_strat``: ConnectionStrategy that links the intermediate cell to the postsynaptic cell.

* ``glom_cell_type``: CellType of the intermediate cell, used as a reference

For the golgi to granule connectivity, the ``axon`` of the golgi cells are selected.

.. code-block:: yaml

      golgi_to_granule:
        strategy: cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus
        presynaptic:
          cell_types:
            - golgi_cell
          morphology_labels:
            - axon
        postsynaptic:
          cell_types:
            - granule_cell
        glom_post_strat: glomerulus_to_granule
        glom_cell_type: glomerulus
        radius: 150
        divergence: 40

.. note::
    Note that here the presynaptic cell is directly connected to the postsynaptic cell and not to
    the intermediate cell. The latter serves only as a reference in the dependent connection
    strategy.

.. _voxel_int:

:doc:`VoxelIntersection <bsb:bsb/bsb.connectivity.detailed>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See bsb :doc:`documentation <bsb:connectivity/connection-strategies>`.
