import unittest

import numpy as np
from neo import SpikeTrain
from quantities import ms

from cerebellar_models.analysis.spiking_results import SpikingResults, extract_isis
from tests.test_spike_plots import ReportBasalSimCircuitFixture


class TestSpikingResults(
    ReportBasalSimCircuitFixture,
    unittest.TestCase,
    engine_name="hdf5",
    setup_cls=True,
):
    def setUp(self):
        super().setUp()
        self.spiking_results = SpikingResults(
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=0,
            time_to=None,
            folder_nio="./",
            ignored_ct=None,
        )

    def test_default_ignored_ct(self):
        self.assertIn("glomerulus", self.spiking_results.ignored_ct)
        self.assertIn("ubc_glomerulus", self.spiking_results.ignored_ct)

    def test_simulation_name_setter_reloads(self):
        # Reassigning simulation_name with both scaffold and folder_nio set
        # must re-trigger load_spikes() and not error out.
        expected_populations = list(self.spiking_results.populations)
        expected_nb_neurons = list(self.spiking_results.nb_neurons)
        self.spiking_results.simulation_name = "basal_activity"
        self.assertEqual(list(self.spiking_results.populations), expected_populations)
        self.assertEqual(list(self.spiking_results.nb_neurons), expected_nb_neurons)

    def test_simulation_name_setter_invalid(self):
        with self.assertRaises(ValueError):
            self.spiking_results.simulation_name = "does_not_exist"

    def test_scaffold_setter_reloads(self):
        # Reassigning scaffold (with simulation_name and folder_nio already set)
        # must re-trigger _check_simulation and load_spikes.
        expected_populations = list(self.spiking_results.populations)
        self.spiking_results.scaffold = self.scaffold
        self.assertEqual(list(self.spiking_results.populations), expected_populations)
        self.assertIs(self.spiking_results.scaffold, self.scaffold)

    def test_scaffold_setter_invalid_simulation(self):
        # If the new scaffold doesn't expose the currently-set simulation_name,
        # the setter must raise. We force the inconsistent state by mutating the
        # private attribute, which simulates a scaffold that was swapped out for
        # one missing the simulation.
        self.spiking_results._simulation_name = "does_not_exist"
        with self.assertRaises(ValueError):
            self.spiking_results.scaffold = self.scaffold

    def test_folder_nio_setter_invalid(self):
        with self.assertRaises(ValueError):
            self.spiking_results.folder_nio = "/this/path/does/not/exist"

    def test_folder_nio_setter_reloads(self):
        # Reassigning folder_nio with both scaffold and simulation_name set
        # must re-trigger load_spikes().
        expected_populations = list(self.spiking_results.populations)
        self.spiking_results.folder_nio = "./"
        self.assertEqual(list(self.spiking_results.populations), expected_populations)

    def test_dt_property(self):
        self.assertEqual(
            self.spiking_results.dt,
            self.scaffold.simulations["basal_activity"].resolution,
        )


class TestExtractISIs(unittest.TestCase):
    def test_extract_isis(self):
        spikes = np.random.random((20, 10)) >= 0.85
        spike_times = np.where(spikes)[0]
        senders = np.where(spikes)[1]
        st = SpikeTrain(
            (spike_times + 1) * 0.1,
            units="ms",
            array_annotations={"senders": senders},
            t_stop=2,
        )
        isis = extract_isis(st, 0.1)
        enough_spikes = np.zeros(10, dtype=bool)
        u, c = np.unique(senders, return_counts=True)
        enough_spikes[u] = c >= 2
        self.assertEqual(len(isis), np.count_nonzero(enough_spikes))
        loc_spikes = spikes[:, enough_spikes]
        for i in range(len(isis)):
            self.assertTrue(
                np.absolute(isis[i] - np.mean(np.diff(np.where(loc_spikes[:, i])[0] * 0.1)) * ms)
                <= 1e-7
            )

    def test_extract_isis_shift(self):
        time_shift = 0.5
        spikes = np.random.random((20, 10)) >= 0.85
        spike_times = np.where(spikes)[0]
        senders = np.where(spikes)[1]
        st = SpikeTrain(
            (spike_times + 1) * 0.1 + time_shift,
            units="ms",
            array_annotations={"senders": senders},
            t_start=time_shift,
            t_stop=2 + time_shift,
        )
        enough_spikes = np.zeros(10, dtype=bool)
        u, c = np.unique(senders, return_counts=True)
        enough_spikes[u] = c >= 2
        loc_spikes = spikes[:, enough_spikes]
        isis = extract_isis(st, 0.1)
        self.assertEqual(len(isis), np.count_nonzero(enough_spikes))
        for i in range(len(isis)):
            self.assertTrue(
                np.absolute(isis[i] - np.mean(np.diff(np.where(loc_spikes[:, i])[0] * 0.1)) * ms)
                <= 1e-7
            )
