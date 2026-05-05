"""
Module for the plots and reports related to the simulation analysis of BSB scaffold.
"""

from typing import List, Tuple, Union

import numpy as np
from bsb import Scaffold
from elephant.kernels import GaussianKernel, Kernel
from matplotlib import gridspec as gs
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from quantities import ms

from .plots import (
    Legend,
    Plot,
    ScaffoldPlot,
)
from .report import (
    BSBReport,
    PlotTypeInfo,
)
from .spiking_results import (
    SpikingResults,
    extract_isis,
    get_correlation_coefficients,
    get_firing_rates,
    get_frequencies,
)
from .structure_analysis import TablePlot


class SpikePlot(ScaffoldPlot):
    """
    Abstract class for plotting the spiking simulation results of a BSB scaffold.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(fig_size, scaffold, dict_colors=dict_colors, **kwargs)
        self.spiking_results = spiking_results
        """Loaded spiking results."""

    @property
    def time_to(self):
        return self.spiking_results.time_to

    @time_to.setter
    def time_to(self, value):
        self.spiking_results.time_to = value

    @property
    def time_from(self):
        return self.spiking_results.time_from

    @time_from.setter
    def time_from(self, value):
        self.spiking_results.time_from = value

    @property
    def populations(self):
        return self.spiking_results.populations

    @property
    def nb_neurons(self):
        return self.spiking_results.nb_neurons

    @property
    def filt_spikes(self):
        return self.spiking_results.filt_spikes

    @property
    def dt(self):
        return self.spiking_results.dt


class SpikeSimulationReport(BSBReport):
    """
    Abstract class for reports of simulation results of BSB scaffold.
    """

    def __init__(
        self,
        scaffold: Union[str, Scaffold],
        simulation_name: str,
        folder_nio: str,
        time_from: float = 0,
        time_to: float = None,
        ignored_ct=None,
        cell_types_info: List[PlotTypeInfo] = None,
    ):
        super().__init__(scaffold, cell_types_info)
        self.spiking_results = SpikingResults(
            scaffold=scaffold,
            simulation_name=simulation_name,
            time_from=time_from,
            time_to=time_to,
            folder_nio=folder_nio,
            ignored_ct=ignored_ct,
        )

    @property
    def time_to(self):
        return self.spiking_results.time_to

    @time_to.setter
    def time_to(self, value):
        self.spiking_results.time_to = value
        for plot in self.plots.values():
            if isinstance(plot, SpikePlot):
                plot.time_to = value

    @property
    def time_from(self):
        return self.spiking_results.time_from

    @time_from.setter
    def time_from(self, value):
        self.spiking_results.time_from = value
        for plot in self.plots.values():
            if isinstance(plot, SpikePlot):
                plot.time_from = value

    @property
    def nb_neurons(self):
        return self.spiking_results.nb_neurons

    @property
    def populations(self):
        return self.spiking_results.populations

    @property
    def dt(self):
        return self.spiking_results.dt

    def add_plot(self, name: str, plot: Plot):
        super().add_plot(name, plot)
        if isinstance(plot, SpikePlot):
            plot.spiking_results = self.spiking_results
            plot.is_updated = False
            if plot.is_plotted:
                plot.clear()

    @property
    def filt_spikes(self):
        """
        Filter the spike events for the time of the analysis.

        :return: Sliced List of SpikeTrain.
        :rtype: List[neo.core.SpikeTrain]
        """
        return self.spiking_results.filt_spikes


class RasterPSTHPlot(SpikePlot):
    """
    Combined raster plot and PSTH plot of the spiking activity results for each neuron type.
    The subplots are split in two columns.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        nb_bins: int = 30,
        dict_colors: dict = None,
        **kwargs,
    ):
        self.spiking_results = spiking_results
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        if nb_bins <= 0:
            raise ValueError("nb_bins must be greater than 0.")
        self.nb_bins = nb_bins
        """Number of bins for the PSTH subplot."""

    def init_plot(self, **kwargs):
        if self.is_initialized:
            plt.close(self.figure)
        self.is_initialized = True
        self.is_plotted = False
        self.nb_cols = 2
        num_filter = len(self.populations)
        self.nb_rows = int(np.ceil(num_filter / 2.0))  # nb rows
        self.figure = plt.figure(figsize=self.fig_size, **kwargs)
        if self.nb_rows > 0:
            global_gsp = gs.GridSpec(self.nb_rows, 2)
        self.axes = [[] for _ in range(self.nb_rows)]
        for i in range(num_filter):
            local_gsp = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=global_gsp[i], hspace=0)
            ax1 = plt.Subplot(self.figure, local_gsp[0])
            self.figure.add_subplot(ax1)

            ax2 = plt.Subplot(self.figure, local_gsp[1])
            self.figure.add_subplot(ax2)
            self.axes[i // 2].append([ax1, ax2])

    def clear(self):
        for ax in self.get_axes():
            ax[0].clear()
            ax[1].clear()
        self.is_plotted = False

    def plot(self, relative_time=False, params_raster: dict = None, params_psth: dict = None):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param bool relative_time: If True, the x-axis values will be relative to the time interval.
        :param params_raster: Dictionary of parameters for the raster plot (see matplotlib scatter).
        :param params_psth: Dictionary of parameters for the PSTH plot (see matplotlib hist).
        """
        super().plot()

        # extract dict params
        loc_params_raster = {"marker": "o", "alpha": 1, "rasterized": True}
        if params_raster is not None:
            loc_params_raster.update(params_raster)
        params_psth = params_psth if params_psth is not None else {}

        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)

        bin_times = np.linspace(self.time_from, self.time_to, self.nb_bins)
        loc_spikes = self.filt_spikes
        for i, ct in enumerate(self.populations):
            times = loc_spikes[i].magnitude
            _, newIds = np.unique(loc_spikes[i].array_annotations["senders"], return_inverse=True)
            cell_params = loc_params_raster.copy()
            if "s" not in cell_params and self.nb_neurons[i] > 0:
                cell_params["s"] = 50.0 / self.nb_neurons[i]
            color = (
                self.labelled_dict_colors[ct][:3]
                if ct in self.labelled_dict_colors
                else [0.6, 0.6, 0.6]
            )
            ax = self.get_ax(i)[0]
            if self.nb_neurons[i] > 0:
                ax.scatter(
                    times,
                    newIds,
                    color=color,
                    **cell_params,
                )
            ax.invert_yaxis()
            ax.set_xlim(
                [0, self.time_to - self.time_from]
                if relative_time
                else [self.time_from, self.time_to]
            )
            ax.get_xaxis().set_visible(False)
            ax.set_ylabel("Neuron id")
            ax.set_title(f"{ct}")

            ax = self.get_ax(i)[1]
            if self.nb_neurons[i] > 0:
                ax.hist(times, bin_times, color=color, **params_psth)
            ax.set_xlabel("Time in ms")
            ax.set_xlim(
                [0, self.time_to - self.time_from]
                if relative_time
                else [self.time_from, self.time_to]
            )
            ax.set_ylabel("Spike counts")


class Spike2Columns(SpikePlot):
    """
    Utility class to plot simulation results for each neuron type in a 2 columns fashion.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        dict_colors: dict = None,
        **kwargs,
    ):
        # population needs to be set before the super.__init__ because it is used in init_plot
        self.spiking_results = spiking_results
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )

    def init_plot(self, **kwargs):
        if self.is_initialized:
            plt.close(self.figure)
        self.is_initialized = True
        self.is_plotted = False
        self.nb_cols = 2
        num_filter = len(self.populations)
        self.nb_rows = int(np.ceil(num_filter / 2.0))  # nb rows
        self.figure = plt.figure(figsize=self.fig_size, **kwargs)
        self.axes = [[] for _ in range(self.nb_rows)]
        for i in range(num_filter):
            self.axes[i // 2].append(
                plt.subplot2grid((self.nb_rows, 2), (i // 2, i % 2), rowspan=1, fig=self.figure)
            )


class FiringRatesPlot(Spike2Columns):
    """
    Instantaneous firing rate plot for each cell type based on a time kernel.
    A firing rate signal is computed as the mean of the convolution of spike times
    for each neuron with the time kernel.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        kernel=None,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        if kernel is not None and not isinstance(kernel, Kernel):
            raise TypeError("Kernel must be an instance of elephant Kernel or None")
        self.kernel = kernel or "auto"
        """Elephant kernel to filter the spike trains"""

    def update(self):
        super().update()
        self.firing_rates = get_firing_rates(self.spiking_results, self.kernel)

    def plot(self, relative_time=False, **kwargs):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param bool relative_time: If True, the x-axis values will be relative to the time interval.
        """
        super().plot()
        time_interval = np.arange(
            self.time_from,
            self.time_to,
            self.dt,
        )
        for i, ct in enumerate(self.populations):
            ax = self.get_ax(i)
            ax.plot(
                time_interval,
                self.firing_rates[:, i],
                color=self.labelled_dict_colors[ct][:3],
                **kwargs,
            )
            ax.set_xlabel("Time in ms")
            ax.set_ylabel("Rate in Hz")
            kernel_text = (
                f" (kernel width = {self.kernel.sigma})" if isinstance(self.kernel, Kernel) else ""
            )
            ax.set_title(f"Mean estimated firing rate for {ct}{kernel_text}")
            ax.set_xlim(
                [0, time_interval[-1] - time_interval[0]]
                if relative_time
                else [time_interval[0], time_interval[-1]]
            )
            ax.text(
                0.01,
                0.95,
                r"FR: {:.2} $\pm$ {:.2}".format(
                    np.mean(self.firing_rates[:, i]), np.std(self.firing_rates[:, i])
                ),
                ha="left",
                va="top",
                transform=ax.transAxes,
            )


class ISIPlot(Spike2Columns):
    """
    Inter-spike interval histogram plot for each cell type.
    For each neuron type, one mean inter-spike interval value is computed for each of its neuron.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        nb_bins: int = 50,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        if nb_bins <= 0:
            raise ValueError("nb_bins must be greater than 0.")
        self.nb_bins = nb_bins
        """Number of bins of the histogram."""

    def plot(self, **kwargs):
        super().plot()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)
        isis_dist = [extract_isis(self.filt_spikes[i], self.dt) for i in range(num_filter)]
        for i, ct in enumerate(self.populations):
            ax2 = self.get_ax(i)
            if len(isis_dist[i]) > 0:
                ax2.hist(
                    isis_dist[i], self.nb_bins, color=self.labelled_dict_colors[ct][:3], **kwargs
                )
            ax2.set_xlabel("ISIs bins in ms")
            ax2.set_yscale("log")
            ax2.set_title(f"Distribution of {ct} ISIs")


class FrequencyPlot(FiringRatesPlot):
    """
    Plot of the frequency distribution analysis of the instantaneous firing rate signal.
    """

    def update(self):
        super().update()
        self.frequencies, self.freq_powers = get_frequencies(
            self.spiking_results, self.firing_rates
        )

    def plot(self, max_freq=30.0, plot_bands=True, **kwargs):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param float max_freq: maximum frequency (in Hz).
        :param bool plot_bands: if True, plot the frequency bands.
        """
        super(FiringRatesPlot, self).plot()
        dict_plot = {}
        dict_plot.update(kwargs)
        for i, (fr, pw, ct) in enumerate(zip(self.frequencies, self.freq_powers, self.populations)):
            ax = self.get_ax(i)
            ax.plot(fr[1:], pw[1:], color=self.labelled_dict_colors[ct], label=ct, **dict_plot)
            ax.set_xlim([0.0, max_freq])
            ax.set_xlabel("Frequency [Hz]")
            ax.set_ylabel("Power [dB]")
            ax.set_title(f"Frequency spectrum for {ct}")
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
            if plot_bands:
                ax.axvline(4.0, ls="--", color="black")
                ax.axvline(8.0, ls="--", color="black")
                ax.axvline(12.0, ls="--", color="black")
                ax.axvline(30.0, ls="--", color="black")


class SimResultsTable(TablePlot, SpikePlot):
    """
    Table of the firing rates and inter-spike intervals for each cell type.
    The firing rate value of a cell type corresponds to the mean number of spike over the time interval,
    while its inter-spike interval corresponds to the mean of all mean inter-spike interval values
    computed for each of its neuron.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        dict_colors: dict = None,
        dict_abv=None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        self.columns = ["Firing rate [Hz]", "Inter Spike Intervals [ms]"]
        self.dict_abv = dict_abv or {}
        """Dictionary of abbreviations for cell types"""

    def plot(self, **kwargs):
        super().plot()
        self.plot_table(**kwargs)

    def update(self):
        super().update()
        self.reset_table()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)
        loc_spikes = self.filt_spikes
        for i in range(num_filter):
            spikes = loc_spikes[i]

            unique_counts = np.unique(spikes.array_annotations["senders"], return_counts=True)[1]
            unique_counts = np.concatenate(
                [unique_counts, np.zeros(self.nb_neurons[i] - len(unique_counts))]
            )
            all_fr = unique_counts / ((self.time_to - self.time_from) / 1000.0)
            isi = extract_isis(spikes, self.dt)

            self._values.append([all_fr, isi])
            self.table_values.append(
                [
                    (
                        "{:.2} pm {:.2}".format(np.mean(all_fr), np.std(all_fr))
                        if len(all_fr) > 0
                        else "/"
                    ),
                    "{:.2} pm {:.2}".format(np.mean(isi), np.std(isi)) if len(isi) > 0 else "/",
                ]
            )
        self.rows = [(self.dict_abv[ct] if ct in self.dict_abv else ct) for ct in self.populations]

    def get_firing_rates(self):
        """
        Return a dictionary which gives for each cell type the firing rate
        of each neuron spiking.
        The plot needs to be updated.

        :rtype: Dict[str, int]
        """
        return {ct: line[0] for ct, line in zip(self.rows, self._values)}

    def get_isis_values(self):
        return {ct: line[1] for ct, line in zip(self.rows, self._values)}


class SpikeCorrelation(SpikePlot):
    """
    Spike cross-correlation matrix plot for each cell type.
    Spike trains will be time binned before computing the pairwise
    Pearson’s correlation coefficients.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        bin_size: float = 5 * ms,
        dict_colors: dict = None,
        dict_abv=None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        self.bin_size = bin_size
        """Size of the time bins used to group spikes before computing correlation coefficients."""
        self.dict_abv = dict_abv or {}
        """Dictionary of abbreviations for cell types"""

    def update(self):
        super().update()
        self.corrcoef = get_correlation_coefficients(self.spiking_results, self.bin_size)

    def plot(self):
        super().plot()
        ax = self.get_ax()
        len_ = len(self.populations)
        im = np.copy(self.corrcoef)
        im[np.tri(len_) > 0] = np.nan
        im = ax.imshow(im, interpolation="nearest")
        ax.set_xticks(np.arange(len_))
        ax.set_xticklabels(
            [self.dict_abv.get(l, l) for l in self.populations],
            rotation=90,
        )
        ax.set_yticks(np.arange(len_))
        ax.set_title("Pearson correlation coef. matrix", fontsize=40)
        ax.set_yticklabels([self.dict_abv.get(l, l) for l in self.populations])
        ax.set_xlabel("Target cell type", fontsize=20)
        ax.set_ylabel("Source cell type", fontsize=20)
        ax_divider = make_axes_locatable(ax)
        cax1 = ax_divider.append_axes("right", size="5%", pad=0)
        self.figure.colorbar(im, cax=cax1)


class SortedPSTH(SpikePlot):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        spiking_results: SpikingResults,
        nb_bins: int = 50,
        sample_size: int = 100,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            spiking_results,
            dict_colors,
            **kwargs,
        )
        if nb_bins <= 0:
            raise ValueError("nb_bins must be greater than 0.")
        self.nb_bins = nb_bins
        """Number of bins for the PSTH."""
        self.sample_size = sample_size
        """Maximum number of neurons to subsample for each population."""

    def plot(self):
        super().plot()

        loc_spikes = self.filt_spikes
        ax = self.get_ax()
        bin_times = np.linspace(self.time_from, self.time_to, self.nb_bins + 1)
        img = np.zeros((self.nb_bins, 0, 4))
        max_delta_times = []
        for i, ct in enumerate(self.populations):
            u_neurons = np.unique(loc_spikes[i].array_annotations["senders"])
            if u_neurons.size > self.sample_size:
                u_neurons = u_neurons[np.random.choice(u_neurons.size, size=self.sample_size)]
            loc_img = np.zeros((self.nb_bins, u_neurons.size, 4))
            for j, neuron in enumerate(u_neurons):
                filter_ = loc_spikes[i].array_annotations["senders"] == neuron
                times = loc_spikes[i].magnitude[filter_]
                # find bin with the largest positive change in spiking activity
                counts = np.histogram(times, bins=bin_times)[0]
                loc_img[:, j, :3] = np.array(
                    self.labelled_dict_colors[ct][:3]
                    if ct in self.labelled_dict_colors
                    else [0.6, 0.6, 0.6]
                )
                loc_img[:, j, 3] += counts
                max_delta_times.append(bin_times[np.argmax(np.diff(counts))])
            img = np.concatenate([img, loc_img], axis=1)

        # sorting
        sorting = np.argsort(max_delta_times)
        img = img[:, : len(max_delta_times)][:, sorting]
        img[:, :, 3] /= np.max(img[:, :, 3], axis=0)
        if np.sum(self.nb_neurons):
            ax.imshow(np.moveaxis(img, 0, 1), interpolation="nearest", aspect="auto")
        ax.set_yticks(np.arange(len(self.populations)))
        ax.set_yticklabels(np.array(self.populations)[sorting])
        ax.set_ylabel("Sorted Neuron id")
        ax.set_xticks(np.arange(0, len(bin_times) - 1, 2))
        ax.set_xticklabels(np.int16(np.round(bin_times[:-1:2])))
        ax.set_xlabel("Time in ms")


class BasicSimulationReport(SpikeSimulationReport):
    """
    Simulation report of the spike activity containing:

    - a plot with the raster and PSTH for each cell type,
    - a table plot storing the mean firing rate and ISI value for each cell type,
    - an instantaneous firing rate plot for each cell type,
    - an inter-spike interval histogram plot for each cell type,
    - a frequency spectrum plot for each cell type,
    - a legend plot
    """

    def __init__(
        self,
        scaffold: Union[str, Scaffold],
        simulation_name: str,
        folder_nio: str,
        time_from: float = 0,
        time_to: float = None,
        ignored_ct=None,
        cell_types_info: List[PlotTypeInfo] = None,
    ):
        super().__init__(
            scaffold, simulation_name, folder_nio, time_from, time_to, ignored_ct, cell_types_info
        )
        num_labelled_ct = len(self.populations)
        raster = RasterPSTHPlot(
            (15, 3 * np.ceil(num_labelled_ct / 2)),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
        )
        table = SimResultsTable(
            (5, 0.22 * (num_labelled_ct + 1)),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
            dict_abv=self.abbreviations,
        )
        firing_rates = FiringRatesPlot(
            (15, 2 * np.ceil(num_labelled_ct / 2)),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
            kernel=GaussianKernel(sigma=20 * ms),
        )
        isis = ISIPlot(
            (15, 2 * np.ceil(num_labelled_ct / 2)),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
        )
        freq = FrequencyPlot(
            (15, 2 * np.ceil(num_labelled_ct / 2)),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
        )
        corr = SpikeCorrelation(
            (10, 10.5),
            scaffold=self.scaffold,
            spiking_results=self.spiking_results,
        )
        legend = Legend(
            (10, 0.6 * num_labelled_ct / 3.0),
            3,
            dict_legend=dict(columnspacing=2.0, handletextpad=0.1, fontsize=20, loc="lower center"),
            dict_abbreviations=self.labelled_abbreviations,
        )
        self.add_plot("raster_psth", raster)
        self.add_plot("table", table)
        self.add_plot("firing_rates", firing_rates)
        self.add_plot("isis", isis)
        self.add_plot("freq", freq)
        self.add_plot("corr", corr)
        self.add_plot("legend", legend)
        legend.dict_colors = raster.labelled_dict_colors.copy()
        legend.remove_ct(self.labelled_cell_names, self.spiking_results.ignored_ct)

    def preprocessing(self):
        self.plots["table"].set_axis_off()
        self.plots["legend"].set_axis_off()
