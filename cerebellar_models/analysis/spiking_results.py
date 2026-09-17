"""
Module for loading spiking results from BSB-NEST simulations and analysing them with
Elephant.
"""

from os import listdir
from os.path import abspath, isdir, isfile, join
from typing import List, Tuple

import numpy as np
from bsb import ResultsError, ResultsMismatchError, Scaffold
from bsb.simulation.results import (
    NetworkRecording,
    RecordedCell,
    ResultsReader,
    SimulationResult,
    SimulationRun,
    read_results,
)
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import correlation_coefficient
from elephant.statistics import instantaneous_rate, isi
from neo import SpikeTrain
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
        folder_nio: str = None,
        ignored_ct: list[str] = None,
        result: SimulationResult = None,
    ):
        """
        :param scaffold: BSB scaffold
        :param simulation_name: name of the simulation
        :param time_from: start time of the analysis
        :param time_to: end time of the analysis
        :param folder_nio: folder where the Neo results are stored. Takes priority
            over ``result`` when both are given.
        :param ignored_ct: List of ignored cells names
        :param result: BSB SimulationResult to load the spike trains from directly,
            without going through a results file. Used as a fallback when
            ``folder_nio`` is not given, e.g. no results file was written yet
            (such as in an ``after_simulation`` hook run in-memory).
        """
        if folder_nio is None and result is None:
            raise ValueError("At least one of 'folder_nio' or 'result' must be provided.")
        self._scaffold = scaffold
        self._folder_nio = None  # will be initialized last
        self._result = None  # will be initialized last
        self.simulation_name = simulation_name
        self.ignored_ct = ignored_ct if ignored_ct is not None else ["glomerulus", "ubc_glomerulus"]
        """List of ignored cell type names"""
        self._run: SimulationRun | None = None
        self._recordings: List[List[NetworkRecording]] = []
        self._all_spikes = []
        self._nb_neurons = np.zeros(0, dtype=int)
        self._populations = []
        # Label combinations actually used in a placement set, cached per cell type
        # name so a population of many cells only computes it once.
        self._label_combos: dict[str, list[tuple[set, np.ndarray]]] = {}
        # Placeholders, resolved for real below once the run (and its provenance) is
        # loaded: _check_times validates against the run's own recorded duration.
        self._time_from = 0.0
        self._time_to = 0.0
        self._dt = None
        # Store the fallback before triggering the initial load, so it's available to
        # _find_run() even when folder_nio (which takes priority) is set.
        self._result = result
        if folder_nio is not None:
            self.folder_nio = folder_nio  # setter triggers load_spikes()
        else:
            self.load_spikes()
        # The run's own provenance, not the scaffold's live configuration, says how
        # long it ran and at what resolution: the configuration may have changed
        # since, and the run is what these results actually came from.
        provenance = self._run.provenance or {}
        self._dt = provenance.get("resolution_ms")
        self._time_from = time_from or 0
        self.time_to = time_to if time_to is not None else provenance.get("duration_ms")

    @staticmethod
    def _check_simulation(scaffold: Scaffold, simulation_name: str):
        """
        Check if a simulation is in a Scaffold and raise an error if not.
        """
        if simulation_name not in scaffold.simulations:
            raise ValueError(f"Simulation name {simulation_name} not in the scaffold simulations")

    @staticmethod
    def _started_at(run: SimulationRun) -> str:
        """
        A run's start time, for sorting runs from most to least recent.

        ISO 8601 timestamps sort correctly as plain strings. A run that recorded
        none sorts before every run that did, rather than failing the comparison.
        """
        return (run.provenance or {}).get("started_at") or ""

    def _find_run(self) -> SimulationRun:
        """
        The most recent run of :attr:`simulation_name` matching :attr:`scaffold`'s
        storage, from ``folder_nio`` or from the in-memory ``result``.

        Every ``.nio`` file in ``folder_nio`` is read with
        :func:`~bsb.read_results`. One produced by a different network raises
        :class:`~bsb.exceptions.ResultsMismatchError`, which is caught here and the
        file simply excluded: a folder can hold results of more than one network,
        and that is the routine case, not an error. Anything else that goes wrong
        reading a file is not caught.

        :raises bsb.exceptions.ResultsError: No run of :attr:`simulation_name`
            matching the scaffold's storage was found.
        :rtype: bsb.simulation.results.SimulationRun
        """
        candidates: List[SimulationRun] = []
        if self.folder_nio is not None:
            for f in sorted(listdir(self.folder_nio)):
                file_ = join(self.folder_nio, f)
                if not (isfile(file_) and file_.endswith(".nio")):
                    continue
                try:
                    reader = read_results(self.scaffold, file_)
                except ResultsMismatchError:
                    continue
                candidates.extend(run for run in reader.runs if run.name == self.simulation_name)
        else:
            reader = ResultsReader(self.scaffold, [self._result.block])
            candidates.extend(run for run in reader.runs if run.name == self.simulation_name)
        if not candidates:
            where = f" in '{self.folder_nio}'" if self.folder_nio is not None else ""
            raise ResultsError(
                f"No results of simulation '{self.simulation_name}' matching the "
                f"scaffold's storage were found{where}."
            )
        return max(candidates, key=self._started_at)

    def _labels_of(self, target: RecordedCell) -> set:
        """
        The labels of a recorded cell, found among the label combinations actually
        used in its placement set.

        :param target: The recorded cell.
        :return: The labels of the cell, or an empty set if it carries none.
        """
        ps_name = target.cell_type.name
        combos = self._label_combos.get(ps_name)
        if combos is None:
            ps = target.placement_set
            combos = [
                (labels, ps.get_label_mask(list(labels)))
                for labels in ScaffoldPlot.get_unique_labels(ps)
            ]
            self._label_combos[ps_name] = combos
        for labels, mask in combos:
            if mask[target.id]:
                return labels
        return set()

    def _population_of(self, recording: NetworkRecording) -> str:
        """
        The population a recording belongs to: its cell type, split further by
        microzone label, or the recording device's own name when it names no
        particular cell (a device recording itself rather than one of its targets).
        """
        target = recording.target
        if recording.kind != "cell" or target is None:
            return recording.device
        return ScaffoldPlot.get_labelled_ct_name(target.cell_type.name, self._labels_of(target))

    def _extract_recordings(self, run: SimulationRun) -> Tuple[List[List[NetworkRecording]], dict]:
        """
        Group a run's spike-train recordings by population.

        :return: - List of Recordings grouped by population: one list of per-cell
                   ``NetworkRecording`` per population, since spikes are recorded
                   one train per cell.
                 - Dictionary storing for each population its index in that list.
        :rtype: Tuple[List[List[NetworkRecording]], Dict[str, int]]
        """
        recordings: List[List[NetworkRecording]] = []
        cell_dict = {}
        for recording in run.recordings():
            if not recording.is_spike_train:
                continue
            population = self._population_of(recording)
            if population not in cell_dict:
                cell_dict[population] = len(recordings)
                recordings.append([])
            recordings[cell_dict[population]].append(recording)

        return recordings, cell_dict

    def load_spikes(self):
        """
        Load the spike trains from the most recent matching run.
        """
        self._run = self._find_run()
        self._recordings, cell_dict = self._extract_recordings(self._run)
        self._populations = list(cell_dict.keys())

        self._fill_lists()

    @staticmethod
    def _sender_id(recording: NetworkRecording, index: int) -> int:
        """
        The id a recording's spikes are attributed to in its population's merged
        train: the cell's own id, or its position within the population for a
        recording that names no particular cell.
        """
        target = recording.target
        return target.id if recording.kind == "cell" and target is not None else index

    #: Per-cell annotations that make no sense on a population's merged train,
    #: since every cell in it disagrees on them; ``senders`` (an array annotation,
    #: one entry per event) is what replaces them.
    _PER_CELL_ANNOTATIONS = ("bsb_cell_id", "bsb_post_cell_id", "bsb_pre_cell_id")

    def _fill_lists(self):
        """
        Merge each population's per-cell recordings into one ``SpikeTrain``, and
        count how many cells were targeted per population.
        """
        duration = (self._run.provenance or {}).get("duration_ms") or 0.0
        self._all_spikes = []  # One Neo SpikeTrain per population
        self._nb_neurons = np.zeros(
            len(self._recordings), dtype=int
        )  # Nb of neurons per population
        for i, group in enumerate(self._recordings):
            self._nb_neurons[i] = len(group)
            times = np.concatenate(
                [recording.signal.times.rescale(ms).magnitude for recording in group]
            )
            senders = np.concatenate(
                [
                    np.full(len(recording.signal), self._sender_id(recording, j))
                    for j, recording in enumerate(group)
                ]
            )
            order = np.argsort(times)
            first = group[0].signal
            annotations = {
                key: value
                for key, value in first.annotations.items()
                if key not in self._PER_CELL_ANNOTATIONS
            }
            self._all_spikes.append(
                SpikeTrain(
                    times[order] * ms,
                    t_stop=duration,
                    name=first.name,
                    array_annotations={"senders": senders[order]},
                    **annotations,
                )
            )

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
        max_time = (self._run.provenance or {}).get("duration_ms") if self._run else None
        if max_time is not None and stop > max_time:
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
            if self._folder_nio is not None or self._result is not None:
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
            if self._folder_nio is not None or self._result is not None:
                self.load_spikes()

    @property
    def dt(self):
        """Time step of the simulation in ms"""
        return self._dt

    @property
    def result(self):
        """BSB SimulationResult the spike trains are loaded from, if any."""
        return self._result

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


def get_spike_matrix(spikes, dt):
    """
    Extract the 2D boolean matrix of the spiking activity for each neuron in the SpikeTrain object.
    Neurons are sorted according to their BSB placement (cell) id.

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
