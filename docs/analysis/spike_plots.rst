Plots and reports of spiking results
====================================

This section describes in details all the spiking plots
available through the cerebellar-models reporting, based on the results of
BSB simulations.

.. _sim_plot:

:class:`SpikePlot <.analysis.spike_plots.SpikePlot>`
----------------------------------------------------
This abstract class provides an interface to plot spike-event based simulation
results for the duration of the simulation. It extends from ``ScaffoldPlot``.

Constructor parameters:

* ``fig_size``: Tuple giving the size of the figure in inches.
* ``spiking_results``: Spiking activity loaded from BSB-NEST simulation (see :ref:`spiking_results`).
* ``dict_colors``: A Dictionary linking the name of the elements (e.g.: cells, fibers) to plot to
  their RGBA color

``spiking_results`` can be obtained from the ``SpikeSimulationReport`` (see :ref:`simulation_report`).

.. _simulation_report:

:class:`SpikeSimulationReport <.analysis.spike_plots.SpikeSimulationReport>`
----------------------------------------------------------------------------
This abstract class provides an interface to create reports for
spike-event based simulation results. It extends from ``BSBReport``.

Constructor parameters:

* ``scaffold``: Scaffold instance or path to the BSB Scaffold file to load
* ``simulation_name``: Name of the simulation as defined in the scaffold
  configuration.
* ``folder_nio``: Folder containing the simulation results stored as nio files.
* ``time_from``: The starting time from which the analysis will be performed
* ``time_to``: The end time at which the analysis will end.
* ``ignored_ct``: List of cell type names to ignore from the nio files
  (results from these cells will not be displayed).
* ``cell_type_info``: List of :class:`PlotTypeInfo <.analysis.report.PlotTypeInfo>`.
  This gives for each element to plot, its name, abbreviation and color.

This class will load the results from nio files produced by the BSB simulation
and store them in the ``spiking_results`` attribute (see :ref:`spiking_results`)

.. note::
   Any modification to the ``time_from`` and ``time_to`` values will be automatically
   forwarded to the report's ``SpikePlots``

.. _raster_psth:

:class:`RasterPSTHPlot <.analysis.spike_plots.RasterPSTHPlot>`
--------------------------------------------------------------
This class extracts and plots a simulation spike times as a raster plot
for the duration of the simulation, as well as the equivalent PSTH for each of
the cell types of the simulation.

On top of the constructor parameters of ``SpikePlot``:

* ``nb_bins``: Numbers of bins for the PSTH plot part for each cell type.


.. _simulation_table:

:class:`SimResultsTable <.analysis.spike_plots.SimResultsTable>`
----------------------------------------------------------------
This class computes the mean firing rate and mean inter-spike interval (ISI)
of each cell type during the simulation time and plot it in a table.
The firing rate value of a cell type corresponds to the mean number of
spike over the simulation time interval, while its ISI
corresponds to the mean of all mean ISI values computed
for each of its neuron.

On top of the constructor parameters of ``SpikePlot``:

* ``dict_abv``: Dictionary that links each cell type name to an abbreviation
  to display

.. _firing_rates_plot:

:class:`FiringRatesPlot <.analysis.spike_plots.FiringRatesPlot>`
----------------------------------------------------------------
This class plots the mean instantaneous firing rate :math:`\lambda (t)`
of each population, expressed according to time, for the duration of
the simulation (see :ref:`firing_rates`).

The final displayed signal :math:`\lambda (t)` corresponds to the mean of
the neurons' :math:`\lambda _m (t)`.

On top of the constructor parameters of ``SpikePlot``:

* ``kernel``: Elephant kernel to use to smoothen the firing rate signals.

.. _isi_plot:

:class:`ISIPlot <.analysis.spike_plots.ISIPlot>`
------------------------------------------------
This class generates the Inter-spike interval (ISI) histogram plot for
each cell type (see also :ref:`isis_distrib`).

On top of the constructor parameters of ``SpikePlot``:

* ``nb_bins``: Numbers of bins for the ISI histogram for each cell type.


.. _frequency_plot:

:class:`FrequencyPlot <.analysis.spike_plots.FrequencyPlot>`
------------------------------------------------------------
This class plots the frequency distribution analysis of the instantaneous
firing rate signal for each cell type (see :ref:`fft_analysis`).

Separators for the major bands of frequencies for neural activity can also
be plotted on top of each panel:

- Delta band: :math:`[0.5; 4]` Hz
- Theta band: :math:`[4; 8]` Hz
- Alpha band: :math:`[8; 12]` Hz
- Beta band: :math:`[12; 30]` Hz
- Gamma band: :math:`[30; 100]` Hz


.. _corr_coef_plot:

:class:`SpikeCorrelationPlot <.analysis.spike_plots.SpikeCorrelationPlot>`
--------------------------------------------------------------------------

Spike cross-correlation matrix plot for each cell type.
Spike trains will be time binned before computing the pairwise coefficient
(see :ref:`pearson_coef`).
On top of the constructor parameters of ``SpikePlot``:

* ``bin_size``: Size of the time bins in ms.
* ``dict_abv``: Dictionary that links each cell type name to an abbreviation
  to display

:class:`SortedPSTH <.analysis.spike_plots.SortedPSTH>`
------------------------------------------------------

This class plots the time-binned instantaneous firing rate for each cell population.
Each population firing rate is obtained from a subsample of the neurons within the population,
and are displayed in rows in the same panel.
The populations are sorted based on the time at which the largest positive change of firing rate
has been detected. This allows you to see how this event is progressing within the circuit during the
simulation.

On top of the constructor parameters of ``SpikePlot``:

* ``nb_bins``: Number of bins for the PSTH.
* ``sample_size``: Maximum number of neurons to subsample for each population


.. _basic_sim_report:

:class:`BasicSimulationReport <.analysis.spike_plots.BasicSimulationReport>`
----------------------------------------------------------------------------
This class extends
:class:`SpikeSimulationReport <.analysis.spike_plots.SpikeSimulationReport>`
and produces a report containing the following
:class:`SpikePlot <.analysis.spike_plots.SpikePlot>` (see section
:ref:`sim_plot`) with a legend:

- A plot showing both the raster plot and Peristimulus Time Histogram (PSTH) for
  the duration of the simulation (see section :ref:`raster_psth`)
- A table containing the mean firing rate and mean InterSpike Intervals (ISIs) for
  each cell type (see section :ref:`simulation_table`).
- A plot showing the mean firing rate according to time of each cell type
  (see section :ref:`firing_rates_plot`)
- A plot showing the ISIs distribution of each cell type (see section
  :ref:`isi_plot`)
- A plot showing the frequency spectrum of each cell type (see section
  :ref:`frequency_plot`)
- A plot showing the cross-correlation of the spiking activity between each cell type
  (see section :ref:`corr_coef_plot`)

All these plots are saved in a single pdf file.
