"""
Module for the plots and reports related to the simulation analysis of BSB scaffold.
"""

from os import listdir
from os.path import isfile, join
from typing import List, Tuple, Union

import numpy as np
from bsb import Scaffold
from elephant.kernels import GaussianKernel, Kernel
from elephant.statistics import instantaneous_rate, isi
from matplotlib import gridspec as gs
from matplotlib import pyplot as plt
from neo import SpikeTrain
from neo import io as nio
from quantities import ms

from cerebellum.analysis.plots import Legend, Plot, ScaffoldPlot
from cerebellum.analysis.report import BSBReport, PlotTypeInfo
from cerebellum.analysis.structure_analysis import TablePlot


def _check_simulation(scaffold: Scaffold, simulation_name: str):
    """
    Check if a simulation is in a Scaffold and raise an error if not.
    """
    if simulation_name not in scaffold.simulations:
        raise ValueError(f"Simulation name {simulation_name} not in the scaffold simulations")


class SpikePlot(ScaffoldPlot):
    """
    Abstract class for plotting the spiking simulation results of a BSB scaffold.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        dict_colors: dict = None,
        **kwargs,
    ):
        _check_simulation(scaffold, simulation_name)
        super().__init__(fig_size, scaffold, dict_colors=dict_colors, **kwargs)
        self.simulation_name = simulation_name
        """Name of the simulation as defined in the scaffold configuration."""
        self._time_from = time_from or 0
        """Start time of the analysis"""
        self.time_to = time_to or self.scaffold.simulations[self.simulation_name].duration
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        self.dt = self.scaffold.simulations[simulation_name].resolution
        """Time step of the simulation in ms"""
        self.all_spikes = all_spikes
        """List of SpikeTrain for each cell type"""
        if len(nb_neurons) != len(populations):
            raise ValueError("populations and nb_neurons must have the same length")
        self.nb_neurons = nb_neurons
        """Number of neuron for each neuron type"""
        self.populations = populations
        """List of neuron type names"""

    def _check_times(self, start, stop):
        if stop < 0 or start < 0:
            raise ValueError("time_from and time_to must be non-negative")
        max_time = self.scaffold.simulations[self.simulation_name].duration
        if stop > max_time:
            raise ValueError("time_to must be less than the simulation's duration")
        if start > stop:
            raise ValueError("time_from must be less than time_to")

    @property
    def time_to(self):
        return self._time_to

    @time_to.setter
    def time_to(self, value):
        self._check_times(self.time_from, value)
        self._time_to = value

    @property
    def time_from(self):
        return self._time_from

    @time_from.setter
    def time_from(self, value):
        self._check_times(value, self.time_to)
        self._time_from = value

    def _set_simulation_params(
        self,
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons,
        populations,
    ):
        is_different = (
            self.simulation_name != simulation_name
            or self.time_from != time_from
            or self.time_to != time_to
            or np.any(self.all_spikes != all_spikes)
            or np.any(self.nb_neurons != nb_neurons)
            or self.populations != populations
        )
        if is_different:
            _check_simulation(self.scaffold, simulation_name)
            self.simulation_name = simulation_name
            self.time_from = time_from
            self.time_to = time_to
            self.all_spikes = all_spikes
            self.nb_neurons = nb_neurons
            self.populations = populations
            self.is_updated = False
            if self.is_plotted:
                self.clear()
        return is_different

    def get_filt_spikes(self):
        """
        Filter the spike events for the time of the analysis.

        :return: Sliced List of SpikeTrain.
        :rtype: List[neo.core.SpikeTrain]
        """
        return [st.time_slice(self.time_from * ms, self.time_to * ms) for st in self.all_spikes]


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
        _check_simulation(self.scaffold, simulation_name)
        self.simulation_name = simulation_name
        """Name of the simulation as defined in the scaffold configuration."""
        self.folder_nio = folder_nio
        """Folder containing the simulation results stored as nio files."""
        self._time_from = time_from
        """Start time of the analysis"""
        self.time_to = time_to or self.scaffold.simulations[simulation_name].duration
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        self.dt = self.scaffold.simulations[simulation_name].resolution
        """Time step of the simulation in ms"""
        self.ignored_ct = ignored_ct or ["glomerulus", "ubc_glomerulus"]
        """List of ignored cell type names"""
        self.all_spikes: List[SpikeTrain] = None
        """List of SpikeTrain for each cell type"""
        self.nb_neurons = None
        """Number of neuron for each neuron type"""
        self.populations = None
        """List of neuron type names"""

        self.all_spikes, self.nb_neurons, self.populations = self.load_spikes()

    def _check_times(self, start, stop):
        if stop < 0 or start < 0:
            raise ValueError("time_from and time_to must be non-negative")
        max_time = self.scaffold.simulations[self.simulation_name].duration
        if stop > max_time:
            raise ValueError("time_to must be less than the simulation's duration")
        if start > stop:
            raise ValueError("time_from must be less than time_to")

    @property
    def time_to(self):
        return self._time_to

    @time_to.setter
    def time_to(self, value):
        self._check_times(self.time_from, value)
        self._time_to = value
        for plot in self.plots.values():
            if isinstance(plot, SpikePlot):
                plot.time_to = value

    @property
    def time_from(self):
        return self._time_from

    @time_from.setter
    def time_from(self, value):
        self._check_times(value, self.time_to)
        self._time_from = value
        for plot in self.plots.values():
            if isinstance(plot, SpikePlot):
                plot.time_from = value

    @staticmethod
    def extract_ct_device_name(cell_type: str):
        """Extract the cell type name from its device name."""
        return cell_type.split("_rec")[0]

    def _extract_spikes_dict(self):
        """
        Extract the spike events from nio files stored in a folder and group them by neuron type.

        :return: - List of spike events grouped by neuron type.
                 - Dictionary storing for each neuron type its index and its unique list of neuron ids.
                   The index is stored under the "id" key and the neuron ids are stored under the "senders" key.
        :rtype: Tuple[List[List[float]], Dict[str, numpy.ndarray[int]]
        """
        spikes_res = []
        cell_dict = {}
        current_id = 0

        for f in listdir(self.folder_nio):
            file_ = join(self.folder_nio, f)
            if isfile(file_) and (".nio" in file_):
                block = nio.NixIO(file_, mode="ro").read_all_blocks()[0]  # assume only one block
                spiketrains = block.segments[0].spiketrains  # assume only one segment

                for st in spiketrains:
                    st.segment = None  # remove spiketrain segment to allow merging
                    cell_type = self.extract_ct_device_name(st.annotations["device"])
                    if cell_type not in self.cell_names:
                        cell_type += "_cell"
                    if cell_type in self.cell_names and cell_type not in self.ignored_ct:
                        if cell_type not in cell_dict:
                            cell_dict[cell_type] = current_id
                            current_id += 1
                            spikes_res.append([])
                        if isinstance(st.annotations["senders"], np.int64):  # pragma: nocover
                            st.annotations["senders"] = [st.annotations["senders"]]
                        if isinstance(st.annotations["senders"], list):
                            # TODO: bsb is not storing senders data in the right dictionary
                            # See https://github.com/dbbs-lab/bsb-nest/issues/18
                            st.array_annotations["senders"] = st.annotations["senders"]
                            del st.annotations["senders"]
                            spikes_res[cell_dict[cell_type]].append(st)
        return spikes_res, cell_dict

    def load_spikes(self):
        """
        Load the spike trains from nio files.

        :return: - Boolean numpy array of shape (N*M) storing spike events for each time step.
                   N corresponds to the number of time steps, M to the number of neuron. Neurons are sorted by type.
                 - List of number of unique neuron per type.
                 - List of cell type names.
        :rtype: Tuple[List[neo.core.SpikeTrain], numpy.ndarray[int], List[str]]
        """
        spikes_res, cell_dict = self._extract_spikes_dict()
        all_spikes = []
        nb_neurons = np.zeros(len(cell_dict), dtype=int)
        for i, cell_type in enumerate(cell_dict):
            sts = spikes_res[cell_dict[cell_type]]
            merged = sts[0]
            for st in sts[1:]:
                merged = merged.merge(st)
            all_spikes.append(merged)
            nb_neurons[i] = all_spikes[i].annotations["pop_size"]
        return all_spikes, nb_neurons, list(cell_dict.keys())

    def add_plot(self, name: str, plot: Plot):
        super().add_plot(name, plot)
        if isinstance(plot, SpikePlot):
            plot._set_simulation_params(
                self.simulation_name,
                self.time_from,
                self.time_to,
                self.all_spikes,
                self.nb_neurons,
                self.populations,
            )

    def get_filt_spikes(self):
        """
        Filter the spike events for the time of the analysis.

        :return: Sliced List of SpikeTrain.
        :rtype: List[neo.core.SpikeTrain]
        """
        return [st.time_slice(self.time_from * ms, self.time_to * ms) for st in self.all_spikes]


class RasterPSTHPlot(SpikePlot):
    """
    Combined raster plot and PSTH plot of the spiking activity results for each neuron type.
    The subplots are split in two columns.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        nb_bins: int = 30,
        dict_colors: dict = None,
        **kwargs,
    ):
        # population needs to be set before the super.__init__ because it is used in init_plot
        self.populations = populations
        super().__init__(
            fig_size,
            scaffold,
            simulation_name,
            time_from,
            time_to,
            all_spikes,
            nb_neurons,
            populations,
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
        loc_spikes = self.get_filt_spikes()
        for i, ct in enumerate(self.populations):
            times = loc_spikes[i].magnitude
            _, newIds = np.unique(loc_spikes[i].array_annotations["senders"], return_inverse=True)
            cell_params = loc_params_raster.copy()
            if "s" not in cell_params and self.nb_neurons[i] > 0:
                cell_params["s"] = 50.0 / self.nb_neurons[i]
            color = self.dict_colors[ct][:3] if ct in self.dict_colors else [0.6, 0.6, 0.6]
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
            ax.set_ylabel("Neuron_Masoli id")
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
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        dict_colors: dict = None,
        **kwargs,
    ):
        # population needs to be set before the super.__init__ because it is used in init_plot
        self.populations = populations
        super().__init__(
            fig_size,
            scaffold,
            simulation_name,
            time_from,
            time_to,
            all_spikes,
            nb_neurons,
            populations,
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
    Each population firing rate signal is plotted surrounding by its standard deviation
    A firing rate signal is computed as the mean of the convolution of spike times
    for each neuron with the time kernel.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        kernel=None,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            simulation_name,
            time_from,
            time_to,
            all_spikes,
            nb_neurons,
            populations,
            dict_colors,
            **kwargs,
        )
        if kernel is not None and not isinstance(kernel, Kernel):
            raise TypeError("Kernel must be an instance of elephant Kernel or None")
        self.kernel = kernel or "auto"
        """Width of the kernel expressed as number of ms"""

    def update(self):
        super().update()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)

        loc_spikes = self.get_filt_spikes()
        duration = int((self.time_to - self.time_from) / self.dt)
        self.firing_rates = np.zeros((duration, num_filter))
        for i in range(num_filter):
            if loc_spikes[i].size <= 0:
                continue  # pragma: nocover
            self.firing_rates[:, i] = (
                instantaneous_rate(
                    loc_spikes[i],
                    sampling_period=self.dt * ms,
                    kernel=self.kernel,
                    border_correction=True,
                ).magnitude[:, 0]
                / self.nb_neurons[i]
            )

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
                color=self.dict_colors[ct][:3],
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
                "FR: {:.2} $\pm$ {:.2}".format(
                    np.mean(self.firing_rates[:, i]), np.std(self.firing_rates[:, i])
                ),
                ha="left",
                va="top",
                transform=ax.transAxes,
            )


def extract_isis(spikes, dt):
    """
    Extract inter-spike intervals from a list of spike trains.
    One mean inter-spike interval value is computed for each neuron.

    :param neo.core.SpikeTrain spikes: population SpikeTrain
    :param float dt: time step

    :return: list of inter-spike intervals
    :rtype: List[float]
    """

    isi_ = []
    senders = spikes.array_annotations["senders"]
    u_senders, inv = np.unique(senders, return_inverse=True)
    mat = np.zeros((int((spikes.t_stop - spikes.t_start) / dt), len(u_senders)), dtype=bool)
    mat[np.asarray(np.rint((spikes.times - spikes.t_start) / dt), dtype=int) - 1, inv] = True
    for sender in range(len(u_senders)):
        isis = isi(np.where(mat[:, sender])[0] * dt * ms)
        if len(isis) > 0:
            isi_.append(np.mean(isis))
    return isi_


class ISIPlot(Spike2Columns):
    """
    Inter-spike interval histogram plot for each cell type.
    For each neuron type, one mean inter-spike interval value is computed for each of its neuron.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        nb_bins: int = 50,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            simulation_name,
            time_from,
            time_to,
            all_spikes,
            nb_neurons,
            populations,
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
        isis_dist = [extract_isis(self.all_spikes[i], self.dt) for i in range(num_filter)]
        for i, ct in enumerate(self.populations):
            ax2 = self.get_ax(i)
            if len(isis_dist[i]) > 0:
                ax2.hist(isis_dist[i], self.nb_bins, color=self.dict_colors[ct][:3], **kwargs)
            ax2.set_xlabel("ISIs bins in ms")
            ax2.set_yscale("log")
            ax2.set_title(f"Distribution of {ct} ISIs")


class FrequencyPlot(FiringRatesPlot):
    """
    Plot of the frequency distribution analysis of the instantaneous firing rate signal.
    """

    def update(self):
        super().update()
        self.frequencies = np.zeros((self.firing_rates.shape[1], self.firing_rates.shape[0] // 2))
        self.freq_powers = np.zeros((self.firing_rates.shape[1], self.firing_rates.shape[0] // 2))
        for i, fr in enumerate(self.firing_rates.T):
            glob_fr = fr[:-1]
            t = np.abs(np.fft.fft(glob_fr))
            x = np.fft.fftfreq(t.shape[0], self.dt / 1e3)  # convert ms to s
            idx = np.argsort(x)
            self.freq_powers[i] = t[idx][t.shape[0] // 2 :] * 2
            self.frequencies[i] = x[idx][x.shape[0] // 2 :]

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
            ax.plot(fr[1:], pw[1:], color=self.dict_colors[ct], label=ct, **dict_plot)
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
        simulation_name: str,
        time_from: float,
        time_to: float,
        all_spikes,
        nb_neurons: List,
        populations: List,
        dict_colors: dict = None,
        dict_abv=None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold,
            simulation_name,
            time_from,
            time_to,
            all_spikes,
            nb_neurons,
            populations,
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
        loc_spikes = self.get_filt_spikes()
        for i in range(num_filter):
            spikes = loc_spikes[i]

            all_fr = np.unique(spikes.array_annotations["senders"], return_counts=True)[1] / (
                (self.time_to - self.time_from) / 1000.0
            )
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

        raster = RasterPSTHPlot(
            (15, 3 * np.ceil(len(self.nb_neurons) / 2)),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        table = SimResultsTable(
            (5, 0.25 * (len(self.nb_neurons) + 1)),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
            dict_abv=self.abbreviations,
        )
        firing_rates = FiringRatesPlot(
            (15, 2 * np.ceil(len(self.nb_neurons) / 2)),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
            kernel=GaussianKernel(sigma=20 * ms),
        )
        isis = ISIPlot(
            (15, 2 * np.ceil(len(self.nb_neurons) / 2)),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        freq = FrequencyPlot(
            (15, 2 * np.ceil(len(self.nb_neurons) / 2)),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        legend = Legend(
            (10, 2),
            3,
            dict_legend=dict(columnspacing=2.0, handletextpad=0.1, fontsize=20, loc="lower center"),
        )
        self.add_plot("raster_psth", raster)
        self.add_plot("table", table)
        self.add_plot("firing_rates", firing_rates)
        self.add_plot("isis", isis)
        self.add_plot("freq", freq)
        self.add_plot("legend", legend)
        legend.remove_ct(self.cell_names, self.ignored_ct)

    def preprocessing(self):
        self.plots["table"].set_axis_off()
        self.plots["legend"].set_axis_off()
