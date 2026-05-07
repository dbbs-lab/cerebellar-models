Spike results extraction and analysis
=====================================

This section describes the methods to extract spiking results from
BSB NEST simulations and to analyze these results using the
`Elephant <https://elephant.readthedocs.io/en/latest/>`_ library.

.. _spiking_results:

:class:`SpikingResults <.analysis.spiking_results.SpikingResults>`
------------------------------------------------------------------

This class allows you to extract and manipulate spiking results coming from BSB-NEST simulations.

Class parameters
^^^^^^^^^^^^^^^^

* :attr:`scaffold<.analysis.spiking_results.SpikingResults.scaffold>`: BSB scaffold object.
* :attr:`simulation_name<.analysis.spiking_results.SpikingResults.simulation_name>`: Name of the
  simulation as defined in the scaffold configuration.
* :attr:`time_from<.analysis.spiking_results.SpikingResults.time_from>`: The starting time from which
  the analysis will be performed
* :attr:`time_to<.analysis.spiking_results.SpikingResults.time_to>`: The end time at which the analysis will end.
* :attr:`folder_nio<.analysis.spiking_results.SpikingResults.folder_nio>`: Path to folder containing
  the `.nio` files produced by the BSB-NEST simulation
* :attr:`ignored_ct<.analysis.spiking_results.SpikingResults.ignored_ct>`: List of cell type names to
  ignore from the nio files (results from these cells will not be returned).

.. note::
   The ``time_from`` and ``time_to`` values can differ from the simulation start
   and end time, but they should remain within the simulation time interval.

The class will automatically load the spikes from the files for you. You can then use it to analyze your
results.

.. note::

    Note that you can change **at any time** all the parameters used to filter the results.

Class properties
^^^^^^^^^^^^^^^^
The following properties are all lists providing different information for each neuron population
(lists have the same size). They will be updated if the class parameters are modified.

* :attr:`filt_spikes<.analysis.spiking_results.SpikingResults.filt_spikes>`: List of
  :class:`SpikeTrain <neo.core.SpikeTrain>`, one for each population, filtered according to the class parameters.
* :attr:`nb_neurons<.analysis.spiking_results.SpikingResults.nb_neurons>`: Number of neurons for each population.
* :attr:`populations<.analysis.spiking_results.SpikingResults.populations>`: The name of each neuron population.


Methods to analyze the spiking results
--------------------------------------

.. _firing_rates:

:func:`get_firing_rates<.analysis.spiking_results.get_firing_rates>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Computes the instantaneous firing rate for each cell population.
Parameters:

* ``spiking_results``: SpikingResults instance.
* ``kernel``: Size of time bin


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

To avoid the edge effects of the kernel convolution with the spike train
(i.e. the time where the kernel can not fully overlap the spike train
because of its width), a compensation effect is calculated.

Different kernel functions would have different smoothening properties.
Here you can use any of the
`Kernels <https://elephant.readthedocs.io/en/latest/reference/kernels.html>`_
provided by Elephant. By default, if no kernel is provided, we use the
optimized gaussian kernel of Elephant.
See also the function
:doc:`instantaneous_rate <elephant:reference/_toctree/statistics/elephant.statistics.instantaneous_rate>`
for more information.

:func:`get_spike_matrix<.analysis.spiking_results.get_spike_matrix>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract the 2D boolean matrix of the spiking activity for each neuron of a SpikeTrain object.
Each line corresponds to a neuron, while the columns correspond to a time step in the simulation.
Matrix box is True if the neuron of this line spiked at the corresponding time step.
Neurons are sorted according to their NEST id.
Parameters:

* ``spikes``: :class:`SpikeTrain <neo.core.SpikeTrain>` instance of a cell population.
* ``dt``: Simulation time step.

.. _isis_distrib:

:func:`extract_isis<.analysis.spiking_results.extract_isis>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract the Inter Spike Intervals (ISI) from a SpikeTrain instance.
Parameters:

* ``spikes``: :class:`SpikeTrain <neo.core.SpikeTrain>` instance of a cell population.
* ``dt``: Simulation time step.

An ISI corresponds to the time (in ms) between two consecutive spikes.
The values extracted corresponds to the mean ISI value of each of its neuron.
Only the neurons spiking at least two times during the SpikeTrain interval are used.

.. _fft_analysis:

:func:`get_frequencies<.analysis.spiking_results.get_frequencies>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This function performs a
:doc:`Fast Fourier Transform <scipy:reference/generated/scipy.fftpack.fftfreq>`
on a instantaneous firing rate signal.
Parameters:

* ``spiking_results``: SpikingResults instance.
* ``firing_rates``: instantaneous firing rates for each neuron population, for each time step,
  corresponds to the output of the :ref:`firing_rates` function.

.. _pearson_coef:

:func:`get_correlation_coefficients<.analysis.spiking_results.get_correlation_coefficients>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This method returns the 2D matrix of the pairwise
:doc:`Pearson’s correlation coefficients <elephant:reference/_toctree/spike_train_correlation/elephant.spike_train_correlation.correlation_coefficient>`
between each cell type.
Parameters:

* ``spiking_results``: SpikingResults instance.
* ``bin_size``: Size of time bin

References
----------

.. [#nawrot_1999] Nawrot, M., Aertsen, A., & Rotter, S. (1999). Single-trial estimation of neuronal firing rates:
   from single-neuron spike trains to population activity. Journal of neuroscience methods, 94(1), 81-92.
   https://doi.org/10.1016/S0165-0270(99)00127-2