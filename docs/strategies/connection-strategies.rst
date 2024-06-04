Connection strategies
---------------------

.. _mossy_glom:

:class:`ConnectomeMossyGlomerulus <.connectome.to_glomerulus.ConnectomeMossyGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects one mossy fiber within the :math:`x\_length \times y\_length \times z\_size` μm box surrounding
each glomerulus. Here ``z_size`` is the size of the partition where the glomeruli are (i.e. no limit).
This selection is random and performed with a truncated exponential distribution.

* ``x_length``: Length of the box along the x axis surrounding the glomerulus in which the mossy fiber can be connected.

* ``y_length``: Length of the box along the y axis surrounding the glomerulus in which the mossy fiber can be connected.

.. note::
    Since the placement of mossy fibers and glomerulus is uniform within their partition,
    the convergence and divergence ratios correspond to the density ratios.

.. _glom_goc:

:class:`ConnectomeGlomerulusGolgi <.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects all glomeruli within the sphere surrounding each Golgi cell soma. For each unique Golgi cell to
connect, the tip of one of its ``basal dendrites`` is selected. This selection is random and performed with a truncated
exponential distribution based on the distance between the tip of each dendrite and the glomerulus to connect.

* ``radius``: Radius of the sphere to filter the glomeruli within it.

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
    As for the glomeruli and Golgi cells, the presynaptic cell of this connection is supposed to
    have no morphology attached, while the postsynaptic cell must have one.

.. _glom_grc:

:class:`ConnectomeGlomerulusGranule <.connectome.glomerulus_granule.ConnectomeGlomerulusGranule>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects n unique glomeruli within the sphere surrounding each granule cell soma, where n is equal to the
``convergence``. Moreover, for each granule cell, each selected glomerulus should be connected through :ref:`mossy_glom`
to a different mossy fiber.
The n unique mossy fibers, their connected glomerulus, and the granule cell ``dendrite`` they will target are all
randomly chosen. If not enough glomeruli could be found in the sphere surrounding the granule cell soma, the closest
glomerulus connected to one of the not yet selected mossy fibers is chosen.

* ``radius``: Radius of the sphere to filter the glomeruli within it.

* ``convergence``: Number of glomeruli per granule cell.

* ``mf_glom_strat``: ConnectionStrategy that links the mossy fibers to the glomeruli.

* ``mf_cell_type``: CellType of the mossy fiber, used as a reference.

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
    have no morphology attached, while the postsynaptic cell must have one.

.. _goc_glom:

:class:`ConnectomeGolgiGlomerulus <.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This connection strategy links a Golgi cell to all granule cells connected to a glomerulus.
The connections between the glomerulus and the postsynaptic cells are defined in :ref:`glom_grc` used as a
``reference strategy``.

The algorithm selects here the closest glomeruli (maximum ``divergence``) that are within
a sphere surrounding each Golgi cell. For each unique glomerulus selected, the tip of an axon's branch from the Golgi
cell is randomly selected. All granule cells that are connected to the selected glomerulus with the
``reference strategy``, are also connected to the selected Golgi cell. The target points of the granule cell (i.e.
dendrite selected) are copied from the ``reference strategy``.

* ``divergence``: Number of glomeruli per Golgi cell.

* ``radius``: Radius of the sphere to filter the glomeruli within it.

* ``glom_post_strat``: ConnectionStrategy that links the glomeruli to the granule cells.

* ``glom_cell_type``: CellType of the glomerulus, used as a reference

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
    Note that here the golgi cell is directly connected to the granule cell and not to the glomerulus.
    The latter serves only as a intermediate point to look for in the ``reference strategy``.

.. _voxel_int:

:doc:`VoxelIntersection <bsb:bsb/bsb.connectivity.detailed>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See bsb :doc:`documentation <bsb:connectivity/connection-strategies>`.

.. _ubc_glom:

:class:`ConnectomeUBCGlomerulus <.connectome.to_glomerulus.ConnectomeUBCGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The algorithm selects one UBC within the sphere surrounding each UBC glomerulus. This selection is random and performed
with a truncated exponential distribution.

* ``radius``: Radius of the sphere to filter the UBC within it.

.. _glom_ubc:

:class:`ConnectomeGlomerulusUBC <.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm connects a population of UBC to a list of glomeruli populations (for instance mf or ubc based glomeruli),
maintaining their relative indegree ratios.
Glomeruli are selected within the sphere surrounding each UBC chunk.
For each UBC, the algorithm selects the closest unconnected glomerulus (if possible) to limit the number of UBCs per
glomerulus.

* ``radius``: Radius of the sphere to filter the glomeruli within it.

* ``ratios_ubc``: Positive relative ratios of indegree for the different types of glomeruli. The ratios will be
  normalized so that their sum equals to 1.

.. code-block:: yaml

    glomerulus_to_ubc:
        strategy: cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulus_to_UBC
        presynaptic:
          cell_types:
            - ubc_glomerulus
            - glomerulus
        postsynaptic:
          cell_types:
            - unipolar_brush_cell
        ratios_ubc:
          ubc_glomerulus: 1. # will be interpreted as 1 / (1+2)
          glomerulus: 2.  # will be interpreted as 2 / (1+2)
        radius: 50
