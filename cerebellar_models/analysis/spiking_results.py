"""
Module for the plots and reports related to the simulation analysis of BSB scaffold.
"""

from os import listdir
from os.path import isfile, join
from typing import List, Tuple

import numpy as np
from bsb import Scaffold
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import correlation_coefficient
from elephant.statistics import instantaneous_rate, isi
from neo import SpikeTrain
from neo import io as nio
from quantities import ms

from cerebellar_models.analysis.plots import ScaffoldPlot


class SpikingResults:

    def __init__(
        self,
        scaffold: Scaffold,
        simulation_name,
        time_from: float,
        time_to: float,
        folder_nio,
        ignored_ct,
    ):
        self._scaffold = scaffold
        self.simulation_name = simulation_name
        """Name of the simulation as defined in the scaffold configuration."""
        self._time_from = time_from or 0
        """Start time of the analysis"""
        self.time_to = time_to or self.scaffold.simulations[self.simulation_name].duration
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        self._dt = self.scaffold.simulations[simulation_name].resolution
        """Time step of the simulation in ms"""
        self.folder_nio = folder_nio
        """Folder containing the simulation results stored as nio files."""
        self.ignored_ct = ignored_ct if ignored_ct is not None else ["glomerulus", "ubc_glomerulus"]
        """List of ignored cell type names"""
        self._all_spikes = []
        """List of SpikeTrain for each cell type"""
        self._nb_neurons = np.zeros(0, dtype=int)
        """Number of neuron for each neuron type"""
        self._populations = []
        """List of neuron type names"""
        self.load_spikes()

    @staticmethod
    def _check_simulation(scaffold: Scaffold, simulation_name: str):
        """
        Check if a simulation is in a Scaffold and raise an error if not.
        """
        if simulation_name not in scaffold.simulations:
            raise ValueError(f"Simulation name {simulation_name} not in the scaffold simulations")

    def _extract_ct_device_name(self, device_name: str):
        """Extract the cell type name from its device name."""
        if "_record" in device_name:
            targetting = (
                self.scaffold.simulations[self.simulation_name].devices[device_name].targetting
            )
            ct = targetting.cell_models[0].name
            labels = targetting["labels"] if "labels" in targetting else set()
            return ct, labels
        else:
            return device_name, set()

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
                    cell_type, labels = self._extract_ct_device_name(st.annotations["device"])
                    if cell_type in list(self.scaffold.cell_types.keys()):
                        cell_type_label = ScaffoldPlot.get_labelled_ct_name(cell_type, labels)
                        if cell_type_label not in cell_dict:
                            cell_dict[cell_type_label] = current_id
                            current_id += 1
                            spikes_res.append([])
                        if "senders" in st.array_annotations:
                            spikes_res[cell_dict[cell_type_label]].append(st)
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
        self._all_spikes = []
        self._nb_neurons = np.zeros(len(cell_dict), dtype=int)
        for i, cell_type in enumerate(cell_dict):
            sts = spikes_res[cell_dict[cell_type]]
            merged = sts[0]
            for st in sts[1:]:
                merged = merged.merge(st)
            self._all_spikes.append(merged)
            self._nb_neurons[i] = self._all_spikes[i].annotations["pop_size"]
        self._populations = list(cell_dict.keys())

    @property
    def filt_spikes(self) -> List[SpikeTrain]:
        """
        Filter the spike events for the time of the analysis.

        :return: Sliced List of SpikeTrain.
        :rtype: List[neo.core.SpikeTrain]
        """
        return [
            sp.time_slice(self.time_from * ms, self.time_to * ms)
            for sp, pop in zip(self._all_spikes, self._populations)
            if pop not in self.ignored_ct
        ]

    @property
    def nb_neurons(self) -> np.ndarray:
        return self._nb_neurons[~np.isin(self._populations, self.ignored_ct)]

    @property
    def populations(self) -> List[str]:
        return [pop for pop in self._populations if pop not in self.ignored_ct]

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

    @property
    def simulation_name(self):
        return self._simulation_name

    @simulation_name.setter
    def simulation_name(self, simulation_name: str):
        self._check_simulation(self.scaffold, simulation_name)
        self._simulation_name = simulation_name

    @property
    def scaffold(self):
        return self._scaffold

    @property
    def dt(self):
        return self._dt


def get_firing_rates(spiking_results, kernel=None):
    num_filter = len(spiking_results.nb_neurons)
    counts = np.zeros(num_filter + 1)
    counts[1:] = np.cumsum(spiking_results.nb_neurons)

    loc_spikes = spiking_results.filt_spikes
    duration = int((spiking_results.time_to - spiking_results.time_from) / spiking_results.dt)
    firing_rates = np.zeros((duration, num_filter))
    for i in range(num_filter):
        if loc_spikes[i].size <= 0:
            continue  # pragma: nocover
        firing_rates[:, i] = (
            instantaneous_rate(
                loc_spikes[i],
                sampling_period=spiking_results.dt * ms,
                kernel=kernel,
                border_correction=True,
            ).magnitude[:, 0]
            / spiking_results.nb_neurons[i]
        )
    return firing_rates


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


def get_frequencies(spiking_results, firing_rates):
    frequencies = np.zeros((firing_rates.shape[1], firing_rates.shape[0] // 2))
    freq_powers = np.zeros((firing_rates.shape[1], firing_rates.shape[0] // 2))
    for i, fr in enumerate(firing_rates.T):
        glob_fr = fr[:-1]
        t = np.abs(np.fft.fft(glob_fr))
        x = np.fft.fftfreq(t.shape[0], spiking_results.dt / 1e3)  # convert ms to s
        idx = np.argsort(x)
        freq_powers[i] = t[idx][t.shape[0] // 2 :] * 2
        frequencies[i] = x[idx][x.shape[0] // 2 :]
    return frequencies, freq_powers


def get_correlation_coefficients(spiking_results, bin_size):
    filt_spikes = spiking_results.filt_spikes
    return (
        correlation_coefficient(
            BinnedSpikeTrain(filt_spikes, bin_size=bin_size),
        )
        if len(filt_spikes) > 0
        else np.zeros((0, 0))
    )
