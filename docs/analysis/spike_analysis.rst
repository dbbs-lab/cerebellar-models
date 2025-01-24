Spike analysis and report of simulation results
===============================================

This section describes in details all the spiking analysis and plots
available through the cerebellum reporting, based on the results of
BSB simulations.

.. _sim_plot:

:class:`SpikePlot <.analysis.spiking_results.SpikePlot>`
--------------------------------------------------------
This abstract class provides an interface to plot spike-event based simulation
results for the duration of the simulation. It extends from ``ScaffoldPlot``.

Constructor parameters:

* ``scaffold``: BSB scaffold object.
* ``simulation_name``: Name of the simulation as defined in the scaffold
  configuration.
* ``time_from``: The starting time from which the analysis will be performed
* ``time_to``: The end time at which the analysis will end.
* ``all_spikes``: A Boolean numpy array of shape (N*M) storing spike events for
  each time step. N corresponds to the number of time steps, M to the number of
  neuron. Neurons are sorted by neuron type.
* ``nb_neurons``: A list containing the number of neuron spiking during the
  simulation for each cell type with the same order as ``all_spikes``.
* ``populations``: The list of the cell type names producing spikes during the
  simulation.

It is advised to obtain the ``all_spikes``, ``nb_neurons``, and ``populations``
parameters from a ``SimulationReport`` (see :ref:`simulation_report`).

.. note::
   The ``time_from`` and ``time_to`` values can differ from the simulation start
   and end time, but they should remain within the simulation time interval.

.. _simulation_report:

:class:`SpikeSimulationReport <.analysis.spiking_results.SpikeSimulationReport>`
--------------------------------------------------------------------------------
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
and store them  in the following attributes:

* ``all_spikes``: A Boolean numpy array of shape (N*M) storing spike events for
  each time step. N corresponds to the number of time steps, M to the number of
  neuron. Neurons are sorted by neuron type.
* ``nb_neurons``: A list containing the number of neuron spiking during the
  simulation for each cell type with the same order as ``all_spikes``.
* ``populations``: The list of the cell type names producing spikes during the
  simulation.

The latter parameters will be automatically forwarded to the report's ``SpikePlots``

.. note::
   The ``time_from`` and ``time_to`` values can differ from the simulation start
   and end time, but they should remain within the simulation time interval.
   Any modification to the ``time_from`` and ``time_to`` values will be automatically
   forwarded to the report's ``SpikePlots``

.. _raster_psth:

:class:`RasterPSTHPlot <.analysis.spiking_results.RasterPSTHPlot>`
------------------------------------------------------------------
This class extracts and plots a simulation spike times as a raster plot
for the duration of the simulation, as well as the equivalent PSTH for each of
the cell types of the simulation.

On top of the constructor parameters of ``SpikePlot``:

* ``nb_bins``: Numbers of bins for the PSTH plot part for each cell type.


.. _simulation_table:

:class:`SimResultsTable <.analysis.spiking_results.SimResultsTable>`
--------------------------------------------------------------------
This class computes the mean firing rate and mean ISI of each cell type
during the simulation time and plot it in a table.
The firing rate value of a cell type corresponds to the mean number of
spike over the simulation time interval, while its inter-spike interval
corresponds to the mean of all mean inter-spike interval values computed
for each of its neuron.

On top of the constructor parameters of ``SpikePlot``:

* ``dict_abv``: Dictionary that links each cell type name to an abbreviation
  to display


.. _firing_rates:

:class:`FiringRatesPlot <.analysis.spiking_results.FiringRatesPlot>`
--------------------------------------------------------------------
This class plots the mean instantaneous firing rate :math:`\lambda (t)`
of each population, expressed according to time, for the duration of
the simulation.

The instantaneous firing rate function of a single neuron :math:`m`,
:math:`\lambda _m (t)` corresponds to the one described in Nawrot et al.
(1999) [#nawrot_1999]_ :
:math:`\lambda _m (t) = \displaystyle\sum_{i=1} ^{n} K(t-t_i)`

where :math:`\{t_0,t_1, ..., t_n\}` are the time of the spike events of
the neuron :math:`m` and :math:`K(t)` is a kernel function with the
following properties:

.. math::
   \begin{cases}
      K(t) \ge 0 \\
      \displaystyle\int_{-\infty}^{+\infty} K(t) \,dt = 1 \\
      \displaystyle\int_{-\infty}^{+\infty} t \cdot K(t) \,dt = 0 \\
   \end{cases}

Additionally, we define :math:`\sigma` the width of the kernel (in ms)
:math:`K` as:
:math:`\sigma = \sqrt{ \displaystyle\int_{-\infty}^{+\infty} t^2 \cdot K(t) \,dt }`

Different kernel functions would have different smoothening properties.
Here we are using a normalized version of the
:doc:`triangle <scipy:reference/generated/scipy.signal.windows.triang>`
function from scipy.

To avoid the edge effects of the kernel convolution with the spike train
(i.e. the time where the kernel can not fully overlap the spike train
because of its width), we extract the computed :math:`\lambda _m (t)`
values on the interval :math:`[time\_from + \sigma; time\_to - \sigma]`.

The final displayed signal :math:`\lambda (t)` corresponds to the mean of
the neurons' :math:`\lambda _m (t)` surrounded by its standard deviation
(clamped at 0). We also display the mean and std of :math:`\lambda (t)`
as a text.

On top of the constructor parameters of ``SpikePlot``:

* ``w_single``: Width of the kernel :math:`\sigma` expressed as number
  of time steps.
* ``max_neuron_sampled``: Maximum number of neurons used to compute
  the firing rate signal. It is used to limit the time it takes to
  complete the kernel convolution operation.


.. _isis_distrib:

:class:`ISIPlot <.analysis.spiking_results.ISIPlot>`
----------------------------------------------------
This class generates the Inter-spike interval (ISI) histogram plot for
each cell type.

An ISI corresponds to the time (in ms) between two consecutive spikes.
For each neuron type, the values extracted for the histogram corresponds
to the mean ISI value of each of its neuron. Only the neurons spiking
at least two times during the simulation interval will be used.

On top of the constructor parameters of ``SpikePlot``:

* ``nb_bins``: Numbers of bins for the ISI histogram for each cell type.


.. _frequency_plot:

:class:`FrequencyPlot <.analysis.spiking_results.FrequencyPlot>`
----------------------------------------------------------------
This class plots the frequency distribution analysis of the instantaneous
firing rate signal for each cell type.
The analysis performs a
:doc:`Fast Fourier Transform <scipy:reference/generated/scipy.fftpack.fftfreq>`
on the instantaneous firing rate calculated as defined in :ref:`firing_rates`.

This class uses the same constructor parameters as in :ref:`firing_rates`.

Separators for the major bands of frequencies for neural activity can also
be plotted on top of each panel:

- Delta band: :math:`[0.5; 4]` Hz
- Theta band: :math:`[4; 8]` Hz
- Alpha band: :math:`[8; 12]` Hz
- Beta band: :math:`[12; 30]` Hz


.. _basic_sim_report:

:class:`BasicSimulationReport <.analysis.spiking_results.BasicSimulationReport>`
--------------------------------------------------------------------------------
This class extends
:class:`spike simulation report <.analysis.spiking_results.SpikeSimulationReport>`
and produces a report containing 5
:class:`SpikePlot <.analysis.spiking_results.SpikePlot>` (see section
:ref:`sim_plot`) with a legend:

- A plot showing both the raster plot and Peristimulus Time Histogram (PSTH) for
  the duration of the simulation (see section :ref:`raster_psth`)
- A table containing the mean firing rate and mean InterSpike Intervals (ISIs) for
  each cell type (see section :ref:`simulation_table`).
- A plot showing the mean firing rate according to time of each cell type
  (see section :ref:`firing_rates`)
- A plot showing the ISIs distribution of each cell type (see section
  :ref:`isis_distrib`)
- A plot showing the frequency spectrum of each cell type (see section
  :ref:`frequency_plot`)

All these plots are saved in a single pdf file.


References
----------

.. [#nawrot_1999] Nawrot, M., Aertsen, A., & Rotter, S. (1999). Single-trial estimation of neuronal firing rates:
   from single-neuron spike trains to population activity. Journal of neuroscience methods, 94(1), 81-92.
   https://doi.org/10.1016/S0165-0270(99)00127-2