"""
    Module for the plots and reports related to the simulation analysis of BSB scaffold.
"""

from os import listdir
from os.path import isfile, join
from typing import List, Tuple

import numpy as np
from bsb import Scaffold
from bsb_nest import NestSimulation
from bsb_neuron import NeuronSimulation
from matplotlib import gridspec as gs
from matplotlib import pyplot as plt
from neo import io as nio
from scipy import signal

from cerebellum.analysis.plots import Legend, Plot, ScaffoldPlot
from cerebellum.analysis.report import LIST_CT_INFO, BSBReport, PlotTypeInfo
from cerebellum.analysis.structure_analysis import TablePlot


class SimulationPlot(ScaffoldPlot):
    """
    Abstract class for plotting the simulation results of a BSB scaffold.
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
        **kwargs,
    ):
        assert simulation_name in scaffold.simulations
        assert isinstance(scaffold.simulations[simulation_name], (NestSimulation, NeuronSimulation))
        super().__init__(fig_size, scaffold, **kwargs)
        self.simulation_name = simulation_name
        """Name of the simulation as defined in the scaffold configuration."""
        self.time_from = time_from or 0
        """Start time of the analysis"""
        self.time_to = time_to or self.scaffold.simulations[simulation_name].duration
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        self.dt = self.scaffold.simulations[simulation_name].resolution
        """Time step of the simulation in ms"""
        self.all_spikes = all_spikes
        """Boolean numpy array of shape (N*M) storing spike events for each time step. 
           N corresponds to the number of time steps, M to the number of neuron.
           Neurons are sorted by neuron type"""
        self.nb_neurons = nb_neurons
        """Number of neuron for each neuron type"""
        self.populations = populations
        """List of neuron type names"""

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
            assert simulation_name in self.scaffold.simulations
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


class SimulationReport(BSBReport):
    """
    Abstract class for reports of simulation results of BSB scaffold.
    """

    def __init__(
        self,
        scaffold: str | Scaffold,
        simulation_name: str,
        folder_nio: str,
        time_from: float = 0,
        time_to: float = None,
        ignored_ct=None,
        cell_type_info: List[PlotTypeInfo] = None,
    ):
        super().__init__(scaffold, cell_type_info)
        self.simulation_name = simulation_name
        """Name of the simulation as defined in the scaffold configuration."""
        self.folder_nio = folder_nio
        """Folder containing the simulation results stored as nio files."""
        self.time_from = time_from
        """Start time of the analysis"""
        self.time_to = time_to or self.scaffold.simulations[simulation_name].duration
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        self.dt = self.scaffold.simulations[simulation_name].resolution
        """Time step of the simulation in ms"""
        self.ignored_ct = ignored_ct or ["gloms", "ubc_glomerulus"]
        """List of ignored cell type names"""
        self.all_spikes = None
        """Boolean numpy array of shape (N*M) storing spike events for each time step. 
                   N corresponds to the number of time steps, M to the number of neuron.
                   Neurons are sorted by neuron type"""
        self.nb_neurons = None
        """Number of neuron for each neuron type"""
        self.populations = None
        """List of neuron type names"""

        self.all_spikes, self.nb_neurons, self.populations = self.load_spikes()

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

        cell_names = [ct.split("_cell")[0] for ct in self.cell_names]
        for f in listdir(self.folder_nio):
            file_ = join(self.folder_nio, f)
            if isfile(file_) and (".nio" in file_):
                block = nio.NixIO(file_, mode="ro").read_all_blocks()[0]
                spiketrains = block.segments[0].spiketrains

                for st in spiketrains:
                    cell_type = self.extract_ct_device_name(st.annotations["device"])
                    if cell_type in cell_names and cell_type not in self.ignored_ct:
                        if cell_type not in cell_dict:
                            cell_dict[cell_type] = {"id": current_id, "senders": []}
                            current_id += 1
                            spikes_res.append([])
                        if isinstance(st.annotations["senders"], np.int64):  # wtf ?
                            st.annotations["senders"] = [st.annotations["senders"]]
                        if len(st.annotations["senders"]) > 0:
                            spikes_res[cell_dict[cell_type]["id"]].append(st)
                            cell_dict[cell_type]["senders"].extend(st.annotations["senders"])
        for cell_type in cell_dict:
            cell_dict[cell_type]["senders"] = np.unique(cell_dict[cell_type]["senders"])
        return spikes_res, cell_dict

    def load_spikes(self):
        """
        Load the spike trains from nio files.

        :return: - Boolean numpy array of shape (N*M) storing spike events for each time step.
                   N corresponds to the number of time steps, M to the number of neuron. Neurons are sorted by type.
                 - List of number of unique neuron per type.
                 - List of cell type names.
        :rtype: Tuple[numpy.ndarray[numpy.ndarray[bool]], numpy.ndarray[int], List[str]]
        """
        spikes_res, cell_dict = self._extract_spikes_dict()
        u_gids = []
        u_cell_types = []
        for i, cell_type in enumerate(cell_dict):
            senders = cell_dict[cell_type]["senders"].tolist()
            u_gids.extend(senders)
            u_cell_types.extend([i] * len(senders))
        sorting = np.argsort(u_gids)
        u_gids = np.array(u_gids)[sorting]
        u_cell_types = np.array(u_cell_types)[sorting]

        inv_convert = np.full(np.max(u_gids) + 1, -1)
        for i, u_gid in enumerate(u_gids):
            inv_convert[u_gid] = i

        tot_num_neuron = len(u_gids)
        all_spikes = np.zeros(
            (int((self.time_to - self.time_from) / self.dt) + 1, tot_num_neuron), dtype=bool
        )
        for cell_type in cell_dict:
            for st in spikes_res[cell_dict[cell_type]["id"]]:
                spikes = st.magnitude
                senders = inv_convert[np.array(st.annotations["senders"])]
                filter_spikes = (spikes > self.time_from) * (spikes <= self.time_to)
                spikes = spikes[filter_spikes]
                spikes = np.asarray(np.floor((spikes - self.time_from) / self.dt), dtype=int)
                senders = np.array(senders)[filter_spikes]
                all_spikes[(spikes, senders)] = True
        nb_neurons = np.zeros(len(cell_dict), dtype=int)
        for i, uf in enumerate(cell_dict.keys()):
            nb_neurons[i] = np.where(u_cell_types == i)[0].size
        return all_spikes, nb_neurons, list(cell_dict.keys())

    def add_plot(self, name: str, plot: Plot):
        super().add_plot(name, plot)
        if isinstance(plot, SimulationPlot):
            plot._set_simulation_params(
                self.simulation_name,
                self.time_from,
                self.time_to,
                self.all_spikes,
                self.nb_neurons,
                self.populations,
            )


class RasterPSTHPlot(SimulationPlot):
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
            **kwargs,
        )
        self.nb_bins = nb_bins
        """Number of bins for the PSTH subplot."""

    def init_plot(self, **kwargs):
        self.is_plotted = False
        self.nb_cols = 2
        num_filter = len(self.populations)
        self.nb_rows = int(np.ceil(num_filter / 2.0))  # nb rows
        self.figure = plt.figure(figsize=self.fig_size, **kwargs)
        global_gsp = gs.GridSpec(self.nb_rows, 2)
        self.axes = [[] for _ in range(self.nb_rows)]
        for i in range(num_filter):
            local_gsp = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=global_gsp[i], hspace=0)
            ax1 = plt.Subplot(self.figure, local_gsp[0])
            self.figure.add_subplot(ax1)

            ax2 = plt.Subplot(self.figure, local_gsp[1])
            self.figure.add_subplot(ax2)
            self.axes[i // 2].append([ax1, ax2])

    def plot(self, relative_time=False, **kwargs):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param bool relative_time: If True, the x-axis values will be relative to the time interval.
        """
        super().plot()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)

        bin_times = np.linspace(0, self.time_to - self.time_from, self.nb_bins)

        for i, ct in enumerate(self.populations):
            times, newIds = np.where(self.all_spikes[:, int(counts[i]) : int(counts[i + 1])])

            ax = self.get_ax(i)[0]
            ax.scatter(
                times * self.dt,
                newIds,
                marker="o",
                c=np.repeat([self.dict_colors[ct][:3]], len(times), axis=0),
                s=50.0 / self.nb_neurons[i],
                alpha=1,
                rasterized=True,
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
            ax.hist(times * self.dt, bin_times, color=self.dict_colors[ct][:3])
            ax.set_xlabel("Time in ms")
            ax.set_xlim(
                [0, self.time_to - self.time_from]
                if relative_time
                else [self.time_from, self.time_to]
            )
            ax.set_ylabel("Spike counts")
            # ax.set_title(f'Spike PSTH plot for {order_ct[i]}')


class Simulation2Columns(SimulationPlot):
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
            **kwargs,
        )

    def init_plot(self, **kwargs):
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


class FiringRatesPlot(Simulation2Columns):
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
        w_single=1000,
        max_neuron_sampled=10000,
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
            **kwargs,
        )
        self.w_single = w_single
        self.max_neuron_sampled = max_neuron_sampled

    def update(self):
        super().update()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)

        kernel_single = (
            signal.windows.triang(self.w_single) * 2 / self.w_single
        )  # normalized boxcar kernel for single-trial firing rate

        self.time_interval = np.arange(
            self.time_from + self.w_single * self.dt,
            self.time_to + (-self.w_single + 1) * self.dt,
            self.dt,
        )
        self.firing_rates = np.zeros((self.all_spikes.shape[0] - self.w_single * 2, num_filter))
        self.std_rates = np.zeros((self.all_spikes.shape[0] - self.w_single * 2, num_filter))
        for i in range(num_filter):
            spikes = self.all_spikes[:, int(counts[i]) : int(counts[i + 1])]
            if self.nb_neurons[i] > self.max_neuron_sampled:
                spikes = spikes[
                    :,
                    np.linspace(
                        0, self.nb_neurons[i], self.max_neuron_sampled, endpoint=False, dtype=int
                    ),
                ]
            R = signal.lfilter(kernel_single, 1, spikes, axis=0) / self.dt * 1000.0
            self.firing_rates[:, i] = np.mean(R, axis=1)[self.w_single : -self.w_single]
            self.std_rates[:, i] = np.std(R, axis=1)[self.w_single : -self.w_single]

    def plot(self, relative_time=False, **kwargs):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param bool relative_time: If True, the x-axis values will be relative to the time interval.
        """
        super().plot()
        for i, ct in enumerate(self.populations):
            ax = self.get_ax(i)
            ax.fill_between(
                self.time_interval,
                (np.maximum(0, self.firing_rates[:, i] - self.std_rates[:, i])),
                (self.firing_rates[:, i] + self.std_rates[:, i]),
                alpha=0.5,
                color=self.dict_colors[ct][:3],
            )
            ax.plot(self.time_interval, self.firing_rates[:, i], color=self.dict_colors[ct][:3])
            ax.set_xlabel("Time in ms")
            ax.set_ylabel("Rate in Hz")
            ax.set_title(
                f"Mean estimated firing rate for {ct} (kernel width = {self.w_single * self.dt} ms)"
            )
            ax.set_xlim(
                [0, self.time_interval[-1] - self.time_interval[0]]
                if relative_time
                else [self.time_interval[0], self.time_interval[-1]]
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

    :param numpy.ndarray[numpy.ndarray[bool]] spikes: Boolean numpy array of shape (N*M)
                                                      storing spike events for each time step.
                                                      N corresponds to the number of time steps,
                                                      M to the number of neuron.
    :param float dt: time step

    :return: list of inter-spike intervals
    :rtype: List[float]
    """
    filter_ = np.where(spikes.T)
    u, idx = np.unique(filter_[0], return_index=True)
    times = np.split(filter_[1], idx[1:])

    isi = []
    for k in range(len(u)):
        isis = np.diff(times[k]) * dt
        if len(isis) > 0:
            isi.append(np.mean(isis))
    return isi


class IsisPlot(Simulation2Columns):
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
            **kwargs,
        )
        self.nb_bins = nb_bins
        """Number of bins of the histogram."""

    def plot(self, **kwargs):
        super().plot()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)
        isis_dist = [
            extract_isis(self.all_spikes[:, int(counts[i]) : int(counts[i + 1])], self.dt)
            for i in range(num_filter)
        ]
        for i, ct in enumerate(self.populations):
            ax2 = self.get_ax(i)
            if len(isis_dist[i]) > 0:
                ax2.hist(isis_dist[i], self.nb_bins, color=self.dict_colors[ct][:3], **kwargs)
                ax2.set_xlabel("ISIs bins in ms")
                ax2.set_yscale("log")
                ax2.set_title(f"Distribution of {ct} ISIs")
            else:
                self.set_axis_off([ax2])


class FrequencyPlot(FiringRatesPlot):
    """
    Plot of the frequency spectogram of the firing rate signal.
    """

    def plot(self, max_freq=30.0, plot_bands=True, *args, **kwargs):
        """
        Plot or replot the figure
        Calls the update function if needed.

        :param float max_freq: maximum frequency (in Hz).
        :param bool plot_bands: if True, plot the frequency bands.
        """
        super(FiringRatesPlot, self).plot()
        for i, (fr, ct) in enumerate(zip(self.firing_rates.T, self.populations)):
            glob_fr = fr[:-1]
            t = np.abs(np.fft.fft(glob_fr))
            x = np.fft.fftfreq(t.shape[0], self.dt / 1e3)
            idx = np.argsort(x)
            t = t[idx][t.shape[0] // 2 :] * 2
            x = x[idx][x.shape[0] // 2 :]
            ax = self.get_ax(i)
            ax.plot(x[1:], t[1:], color=self.dict_colors[ct], alpha=0.7, label=ct)
            ax.set_xlim([0.0, max_freq])
            ax.set_xlabel("Frequency [Hz]")
            ax.set_ylabel("Power [dB]")
            ax.set_title(f"Frequency spectrum for {ct}")
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
            if plot_bands:
                ax.axvline(4.0, ls="--", color="black")
                ax.axvline(8.0, ls="--", color="black")
                ax.axvline(12.0, ls="--", color="black")


class SimResultsTable(TablePlot, SimulationPlot):
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
        dict_abv: dict = None,
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
            **kwargs,
        )
        self.columns = ["Firing rate [Hz]", "Inter Spike Intervals [ms]"]
        self.dict_abv = dict_abv or {ct.name: ct.abbreviation for ct in LIST_CT_INFO}

    def plot(self, **kwargs):
        super().plot()
        self.plot_table(**kwargs)

    def update(self):
        super().update()
        self.update_values()
        num_filter = len(self.nb_neurons)
        counts = np.zeros(num_filter + 1)
        counts[1:] = np.cumsum(self.nb_neurons)
        for i in range(num_filter):
            spikes = self.all_spikes[:, int(counts[i]) : int(counts[i + 1])]
            all_fr = np.sum(spikes, axis=0) / ((self.time_to - self.time_from) / 1000.0)
            isi = extract_isis(spikes, self.dt)

            self.values.append(
                [
                    "{:.2} pm {:.2}".format(np.mean(all_fr), np.std(all_fr)),
                    "{:.2} pm {:.2}".format(np.mean(isi), np.std(isi)) if len(isi) > 0 else "/",
                ]
            )
        self.rows = self.populations


class BasicSimulationReport(SimulationReport):
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
        scaffold: str | Scaffold,
        simulation_name: str,
        folder_nio: str,
        time_from: float = 0,
        time_to: float = None,
        ignored_ct=None,
        cell_type_info: List[PlotTypeInfo] = None,
    ):
        super().__init__(
            scaffold, simulation_name, folder_nio, time_from, time_to, ignored_ct, cell_type_info
        )
        raster = RasterPSTHPlot(
            (15, 10),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        table = SimResultsTable(
            (5, 2.5),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        firing_rates = FiringRatesPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        isis = IsisPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        freq = FrequencyPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name=self.simulation_name,
            time_from=self.time_from,
            time_to=self.time_to,
            all_spikes=self.all_spikes,
            nb_neurons=self.nb_neurons,
            populations=self.populations,
        )
        legend = Legend(
            (10, 2.5),
            3,
            dict_legend=dict(columnspacing=2.0, handletextpad=0.1, fontsize=20, loc="lower center"),
        )
        self.add_plot("raster_psth", raster)
        self.add_plot("table", table)
        self.add_plot("firing_rates", firing_rates)
        self.add_plot("isis", isis)
        self.add_plot("freq", freq)
        self.add_plot("legend", legend)
        table.set_axis_off()
        legend.set_axis_off()
