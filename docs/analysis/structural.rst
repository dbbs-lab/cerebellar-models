Structural analysis and report
==============================

The cerebellum repository provides a list of analysis and plot that can be
performed during or after the BSB reconstruction to extract structural
properties of the produced circuit.

The following section will describe in more details all the plots and analysis
available through the cerebellum reporting.

.. _placement_table:

:class:`PlacementTable <.analysis.structure_analysis.PlacementTable>`
---------------------------------------------------------------------
This class computes the number and density of each cell type of the BSB
scaffold and plot them in a table.

Constructor parameters:

* ``scaffold``: BSB scaffold object.
* ``dict_colors``: A Dictionary linking the name of the elements (e.g.: cells,
  fibers) to plot to their RGBA color
* ``dict_abv``: Dictionary that links each cell type name to an abbreviation
  to display


.. _connectivity_table:

:class:`ConnectivityTable <.analysis.structure_analysis.ConnectivityTable>`
---------------------------------------------------------------------------
This class computes and plots a Table of the results of the connectivity for
BSB Scaffold. This includes for each pair of connected cell types:

- the number of synapses formed
- the number of synapses per unique pair of cell.
- the convergence ratio defined as the mean number of afferent
  connections created with a single postsynaptic cell
- the divergence ratio defined as the mean number of efferent
  connections created with a single presynaptic cell

Constructor parameters:

* ``scaffold``: BSB scaffold object.
* ``dict_colors``: A Dictionary linking the name of the elements (e.g.: cells,
  fibers) to plot to their RGBA color
* ``dict_abv``: Dictionary that links each cell type name to an abbreviation
  to display


.. _cell_placement:

:class:`CellPlacement3D <.analysis.structure_analysis.CellPlacement3D>`
-----------------------------------------------------------------------
This class extracts and plots the position of the cells in the scaffold in
3D space.
Each soma diameter is based on its cell type radius as defined in the BSB
configuration.

Constructor parameters:

* ``scaffold``: BSB scaffold object.
* ``dict_colors``: A Dictionary linking the name of the elements (e.g.: cells,
  fibers) to plot to their RGBA color
* ``ignored_ct``: List of cell type names to ignore from the scaffold
  (these cells will not be displayed).


.. _structure_report:

:class:`StructureReport <.analysis.structure_analysis.StructureReport>`
-----------------------------------------------------------------------
This class extends :class:`BSBReport <.analysis.report.BSBReport>` and produces
a report containing 3 plots with a legend:

- A table containing the number and densities of placed cells in the final
  circuit (see section :ref:`placement_table`)
- A table containing the number of synapses, convergence, divergence ratios for
  each of the connected cell type pairs (see section :ref:`connectivity_table`).
- A 3D plot of the circuit, showing neuron somas placed in the circuit space by
  BSB (see section :ref:`cell_placement`)

All these plots are saved in a single pdf file.

You can also create this report automatically through the BSB configuration like
so (see :ref:`run_structure_report`).


.. _run_structure_report:

:class:`RunStructureReport <.analysis.structure_analysis.RunStructureReport>`
-----------------------------------------------------------------------------
BSB postprocessing node to generate a scaffold :ref:`structure_report` after
running the connectivity jobs. To automatically produce this during BSB scaffold
compilation add the following section to your configuration:

.. code-block::yaml

    after_connectivity:
        print_structure_report:
            strategy: cerebellum.analysis.structure_analysis.RunStructureReport
            output_filename: bsb_report_structure.pdf
