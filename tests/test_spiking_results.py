import copy
import shutil
import tempfile
import unittest
from os.path import join
from types import SimpleNamespace

import numpy as np
from bsb import ResultsError
from bsb.simulation.results import cell_annotations, device_annotations
from bsb.storage.provenance import decode_annotation, encode_annotation
from neo import AnalogSignal, SpikeTrain
from neo import io as nio
from quantities import ms, mV

from cerebellar_models.analysis.plots import ScaffoldPlot
from cerebellar_models.analysis.report import LIST_CT_INFO
from cerebellar_models.analysis.spike_plots import SpikePlot
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

    def test_result_property(self):
        self.assertEqual(
            self.spiking_results.result,
            None,
        )

    def test_raise_without_nio_nor_results(self):
        with self.assertRaises(ValueError):
            self.spiking_results = SpikingResults(
                scaffold=self.scaffold,
                simulation_name="basal_activity",
                time_from=0,
                time_to=None,
            )

    def _write_foreign_results(self, filename):
        """
        Copy the streamed results under `filename`, as if another network produced them.
        """
        with nio.NixIO(self.streamed_nio_file, mode="ro") as reader:
            blocks = reader.read_all_blocks()
        for block in blocks:
            provenance = decode_annotation(block.annotations["bsb_provenance"])
            provenance["scaffold"]["storage_id"] = "another-network"
            block.annotations["bsb_provenance"] = encode_annotation(provenance)
        with nio.NixIO(filename, mode="ow") as writer:
            writer.write_all_blocks(blocks)

    def test_folder_skips_results_of_other_networks(self):
        folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, folder)
        self._write_foreign_results(join(folder, "other_network.nio"))
        # A folder holding only another network's results has nothing to load.
        with self.assertRaises(ResultsError):
            SpikingResults(self.scaffold, "basal_activity", 0, None, folder_nio=folder)
        # Next to this network's results, the foreign file is simply ignored.
        shutil.copy(self.streamed_nio_file, join(folder, "this_network.nio"))
        results = SpikingResults(self.scaffold, "basal_activity", 0, None, folder_nio=folder)
        self.assertEqual(results.populations, self.spiking_results.populations)
        self.assertEqual(list(results.nb_neurons), list(self.spiking_results.nb_neurons))

    def test_device_and_analog_recordings(self):
        block = copy.deepcopy(self.simulation_results[0].block)
        segment = block.segments[0]
        # A spike train of the device itself, which names no cell of the network.
        segment.spiketrains.append(
            SpikeTrain(
                [10.0, 20.0, 30.0] * ms,
                t_stop=self.simulation_duration * ms,
                bsb_device_name="device_spikes",
                **device_annotations("record"),
            )
        )
        # A membrane potential recording, which is no spike train to analyze.
        cell_model = self.basal_simulation.cell_models["purkinje_cell"]
        segment.analogsignals.append(
            AnalogSignal(
                np.zeros((10, 1)) * mV,
                sampling_period=0.1 * ms,
                bsb_device_name="purkinje_vm",
                **cell_annotations(cell_model, 0, "record"),
            )
        )
        results = SpikingResults(
            self.scaffold,
            "basal_activity",
            0,
            None,
            ignored_ct=[],
            result=SimpleNamespace(block=block),
        )
        self.assertIn("device_spikes", results.populations)
        self.assertNotIn("purkinje_vm", results.populations)
        idx = results.populations.index("device_spikes")
        self.assertEqual(results.nb_neurons[idx], 1)
        spikes = results.filt_spikes[idx]
        self.assertEqual(list(spikes.magnitude), [10.0, 20.0, 30.0])
        # Its spikes are attributed to its position in the population, as it has no cell id.
        self.assertEqual(list(spikes.array_annotations["senders"]), [0, 0, 0])

    def test_labelled_populations(self):
        ps = self.scaffold.get_placement_set("purkinje_cell")
        ids = np.arange(len(ps))
        plus, minus = ids[: len(ids) // 2], ids[len(ids) // 2 :]
        ps.label(["plus"], plus)
        ps.label(["minus"], minus)
        self.addCleanup(ps.remove_labels, ["minus"], minus)
        self.addCleanup(ps.remove_labels, ["plus"], plus)
        # Read the in-memory run: the streamed file would warn that the network was
        # written to since the run, which is the very point of this test.
        results = SpikingResults(
            self.scaffold, "basal_activity", 0, None, result=self.simulation_results[0]
        )
        self.assertNotIn("purkinje_cell", results.populations)
        nb_neurons = dict(zip(results.populations, results.nb_neurons))
        self.assertEqual(nb_neurons["purkinje_cell_plus"], len(plus))
        self.assertEqual(nb_neurons["purkinje_cell_minus"], len(minus))
        for population, cells in (("purkinje_cell_plus", plus), ("purkinje_cell_minus", minus)):
            senders = results.filt_spikes[results.populations.index(population)].array_annotations[
                "senders"
            ]
            self.assertTrue(np.isin(senders, cells).all(), f"{population} holds foreign cells")

        # The live scaffold and the recorded populations must agree on each shade.
        colors = {ct.name: ct.color for ct in LIST_CT_INFO}
        scaffold_colors = ScaffoldPlot(
            (10, 10), self.scaffold, dict_colors=colors
        ).labelled_dict_colors
        spike_colors = SpikePlot((10, 10), results, dict_colors=colors).labelled_dict_colors
        self.assertNotIn("purkinje_cell", spike_colors)
        for population in ("purkinje_cell_minus", "purkinje_cell_plus"):
            np.testing.assert_allclose(spike_colors[population], scaffold_colors[population])
        np.testing.assert_allclose(
            spike_colors["purkinje_cell_plus"],
            spike_colors["purkinje_cell_minus"] * np.array([2 / 3] * 3 + [1.0]),
        )

    def test_overlapping_labelled_populations(self):
        # The cells labeled {plus, edge} also carry `plus`, so they are found by the
        # label mask of {plus}: they must still be told apart from the {plus} cells.
        ps = self.scaffold.get_placement_set("purkinje_cell")
        ids = np.arange(len(ps))
        edge, plus, minus = np.array_split(ids, 3)
        ps.label(["plus"], np.concatenate([edge, plus]))
        ps.label(["edge"], edge)
        ps.label(["minus"], minus)
        self.addCleanup(ps.remove_labels, ["minus"], minus)
        self.addCleanup(ps.remove_labels, ["edge"], edge)
        self.addCleanup(ps.remove_labels, ["plus"], np.concatenate([edge, plus]))
        results = SpikingResults(
            self.scaffold, "basal_activity", 0, None, result=self.simulation_results[0]
        )
        expected = {}
        for combo in ScaffoldPlot.get_unique_labels(ps):
            cells = {
                frozenset({"plus", "edge"}): edge,
                frozenset({"plus"}): plus,
                frozenset({"minus"}): minus,
            }[frozenset(combo)]
            expected[ScaffoldPlot.get_labelled_ct_name("purkinje_cell", combo)] = cells
        self.assertEqual(len(expected), 3)
        nb_neurons = dict(zip(results.populations, results.nb_neurons))
        for population, cells in expected.items():
            self.assertEqual(nb_neurons.get(population), len(cells), population)
            senders = results.filt_spikes[results.populations.index(population)].array_annotations[
                "senders"
            ]
            self.assertTrue(np.isin(senders, cells).all(), f"{population} holds foreign cells")


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
