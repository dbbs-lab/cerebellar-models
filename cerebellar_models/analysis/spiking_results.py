"""
Module for loading spiking results from BSB-NEST simulations and analysing them with
Elephant.
"""

from os import listdir
from os.path import abspath, isdir, isfile, join
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
    """
    Class used to load the results of the NEST simulations
    produced from a BSB scaffold.
    """

    def __init__(
        self,
        scaffold: Scaffold,
        simulation_name,
        time_from: float,
        time_to: float,
        folder_nio: str,
        ignored_ct: list[str],
    ):
        """
        :param scaffold: BSB scaffold
        :param simulation_name: name of the simulation
        :param time_from: start time of the analysis
        :param time_to: end time of the analysis
        :param folder_nio: folder where the Neo results are stored
        :param ignored_ct: List of ignored cells names
        """
        self._scaffold = scaffold
        self._folder_nio = None  # will be initialized last
        self.simulation_name = simulation_name
        self._time_from = time_from or 0
        self.time_to = time_to or self.scaffold.simulations[self.simulation_name].duration
        self._dt = self.scaffold.simulations[simulation_name].resolution
        self.ignored_ct = ignored_ct if ignored_ct is not None else ["glomerulus", "ubc_glomerulus"]
        """List of ignored cell type names"""
        self._all_spikes = []
        self._nb_neurons = np.zeros(0, dtype=int)
        self._populations = []
        self.folder_nio = folder_nio

    @staticmethod
    def _check_simulation(scaffold: Scaffold, simulation_name: str):
        """
        Check if a simulation is in a Scaffold and raise an error if not.
        """
        if simulation_name not in scaffold.simulations:
            raise ValueError(f"Simulation name {simulation_name} not in the scaffold simulations")

    def _extract_ct_device_name(self, device_name: str):
        """
        Extract the cell type name from its device name.
        """
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

        :return: List of time-sliced SpikeTrain.
        :rtype: List[neo.core.SpikeTrain]
        """
        return [
            sp.time_slice(self.time_from * ms, self.time_to * ms)
            for sp, pop in zip(self._all_spikes, self._populations)
            if pop not in self.ignored_ct
        ]

    @property
    def nb_neurons(self) -> np.ndarray:
        """Number of neuron for each neuron type"""
        return self._nb_neurons[~np.isin(self._populations, self.ignored_ct)]

    @property
    def populations(self) -> List[str]:
        """List of neuron type names"""
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
        """End time of the analysis. By default, this corresponds to the simulation duration."""
        return self._time_to

    @time_to.setter
    def time_to(self, value: float):
        self._check_times(self.time_from, value)
        self._time_to = value

    @property
    def time_from(self):
        """Start time of the analysis"""
        return self._time_from

    @time_from.setter
    def time_from(self, value: float):
        self._check_times(value, self.time_to)
        self._time_from = value

    @property
    def simulation_name(self):
        """Name of the simulation as defined in the scaffold configuration."""
        return self._simulation_name

    @simulation_name.setter
    def simulation_name(self, simulation_name: str):
        self._simulation_name = simulation_name
        if self._scaffold is not None:
            self._check_simulation(self.scaffold, simulation_name)
            if self._folder_nio is not None:
                self.load_spikes()

    @property
    def scaffold(self):
        """BSB Scaffold used as reference for simulation results."""
        return self._scaffold

    @scaffold.setter
    def scaffold(self, scaffold: Scaffold):
        self._scaffold = scaffold
        if self._simulation_name is not None:
            self._check_simulation(self.scaffold, self.simulation_name)
            if self._folder_nio is not None:
                self.load_spikes()

    @property
    def dt(self):
        """Time step of the simulation in ms"""
        return self._dt

    @property
    def folder_nio(self):
        return self._folder_nio

    @folder_nio.setter
    def folder_nio(self, value):
        """Path to folder containing the simulation results stored as nio files."""
        if not isdir(abspath(value)):
            raise ValueError(f"The folder path to nio results cannot be reached: {abspath(value)}")
        self._folder_nio = value
        if self._scaffold is not None and self._simulation_name is not None:
            self.load_spikes()


def get_firing_rates(spiking_results: SpikingResults, kernel=None) -> np.ndarray:
    """
    Get the instantaneous firing rate for each cell type from SpikingResults
    based on a time kernel.

    :param SpikingResults spiking_results: simulation spike results
    :param kernel: Elephant kernel to filter the spike trains
    :return: numpy array storing instantaneous firing rates for each population, for each
        time step.
    :rtype: numpy.ndarray[float]
    """
    num_filter = len(spiking_results.nb_neurons)
    counts = np.zeros(num_filter + 1)
    counts[1:] = np.cumsum(spiking_results.nb_neurons)

    loc_spikes = spiking_results.filt_spikes
    duration = int((spiking_results.time_to - spiking_results.time_from) / spiking_results.dt)
    firing_rates = np.zeros((duration, num_filter))
    for i in range(num_filter):
        if loc_spikes[i].size <= 2:
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


def get_spike_matrix(spikes, dt):
    """
    Extract the 2D boolean matrix of the spiking activity for each neuron in the SpikeTrain object.
    Neurons are sorted according to their NEST id.

    :param neo.core.SpikeTrain spikes: population SpikeTrain object
    :param float dt: time step
    :return: numpy array 2D boolean matrix storing spike events for each neuron, for each time
        step.
    :rtype: numpy.ndarray[bool]
    """
    senders = spikes.array_annotations["senders"]
    u_senders, inv = np.unique(senders, return_inverse=True)
    mat = np.zeros((int((spikes.t_stop - spikes.t_start) / dt), len(u_senders)), dtype=bool)
    mat[np.asarray(np.rint((spikes.times - spikes.t_start) / dt), dtype=int) - 1, inv] = True
    return mat


def extract_isis(spikes, dt):
    """
    Extract inter-spike intervals from a SpikeTrain object.
    One mean inter-spike interval value is computed for each neuron.

    :param neo.core.SpikeTrain spikes: population SpikeTrain object
    :param float dt: time step
    :return: list of inter-spike intervals
    :rtype: List[float]
    """

    isi_ = []
    mat = get_spike_matrix(spikes, dt)
    for sender in range(mat.shape[1]):
        isis = isi(np.where(mat[:, sender])[0] * dt * ms)
        if len(isis) > 0:
            isi_.append(np.mean(isis))
    return isi_


def get_frequencies(spiking_results, firing_rates):
    """
    Get the Fast Fourier Transform on instantaneous firing rates signal.

    :param SpikingResults spiking_results: simulation spike results
    :param firing_rates: numpy array storing instantaneous firing rates for each
        population, for each time step.
    :return: Tuple of lists of frequencies and their corresponding FFT powers.
    :rtype: Tuple[numpy.ndarray[float], numpy.ndarray[float]]
    """
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
    """
    Get the spike cross-correlation matrix for each cell type.
    Spike trains will be time binned before computing the pairwise
    Pearson’s correlation coefficients.

    :param SpikingResults spiking_results: simulation spike results
    :param float bin_size: size of time bin
    :return: numpy array 2D matrix storing the Pearson correlation coefficients
        between each neuron population.
    :rtype: numpy.ndarray[float]
    """
    filt_spikes = spiking_results.filt_spikes
    return (
        correlation_coefficient(
            BinnedSpikeTrain(filt_spikes, bin_size=bin_size),
        )
        if len(filt_spikes) > 0
        else np.zeros((0, 0))
    )
