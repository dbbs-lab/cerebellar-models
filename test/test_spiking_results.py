import os
import unittest

import numpy as np
from bsb import Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture
from matplotlib import pyplot as plt

from cerebellum.analysis.plots import ScaffoldPlot
from cerebellum.analysis.spiking_results import (
    BasicSimulationReport,
    FiringRatesPlot,
    FrequencyPlot,
    ISIPlot,
    RasterPSTHPlot,
    SimResultsTable,
    SpikePlot,
    SpikeSimulationReport,
    extract_isis,
)


class MiniCerebCircuitTest(RandomStorageFixture, unittest.TestCase, engine_name="hdf5"):
    def setUp(self):
        super().setUp()
        # one third of the canonical circuit
        self.cfg = parse_configuration_file("configurations/mouse/nest/stimulus_mossy_vitro.yaml")
        self.cfg.network.x = 100
        self.cfg.network.y = 66
        self.cfg.network.z = 100
        self.cfg.partitions.granular_layer.thickness = 40
        self.cfg.partitions.purkinje_layer.thickness = 10
        self.cfg.partitions.b_molecular_layer.thickness = 17
        self.cfg.partitions.t_molecular_layer.thickness = 33
        # make sure there are enough mfs.
        self.cfg.cell_types.glomerulus.spatial.density = 0.00034

        self.simulation_duration = 1000.0
        self.cfg.simulations.basal_activity.duration = self.simulation_duration
        self.cfg.simulations.mf_stimulus.duration = self.simulation_duration

        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold.compile(skip_after_connectivity=True, clear=True)


class ReportBasalSimCircuitTest(MiniCerebCircuitTest, engine_name="hdf5"):
    def setUp(self):
        super().setUp()
        self.simulation_results = self.scaffold.run_simulation("basal_activity")
        self.simulation_results.write("test_sim_results.nio", "ow")
        self.simulationReport = SpikeSimulationReport(
            self.scaffold,
            "basal_activity",
            "./",
        )

    def tearDown(self):
        super().tearDown()
        plt.close("all")
        os.remove("test_sim_results.nio")


class TestSpikePlots(ReportBasalSimCircuitTest, NumpyTestCase, engine_name="hdf5"):
    def test_spike_reports(self):
        self.assertEqual(self.simulationReport.time_to, self.simulation_duration)
        self.assertEqual(self.simulationReport.dt, 0.1)
        self.assertEqual(
            len(self.simulationReport.nb_neurons), len(self.simulationReport.populations)
        )
        self.assertEqual(
            len(self.simulationReport.all_spikes), int(self.simulation_duration / 0.1) + 1
        )
        self.assertEqual(
            len(self.simulationReport.all_spikes[0]), sum(self.simulationReport.nb_neurons)
        )
        self.assertAll(np.sum(self.simulationReport.all_spikes, axis=0) > 0)
        self.assertTrue("mossy_fibers" in self.simulationReport.populations)
        self.assertTrue("glomerulus" not in self.simulationReport.populations)
        with self.assertRaises(ValueError):
            SpikeSimulationReport(self.scaffold, "blabla", "./")

        empty_report = SpikeSimulationReport(self.scaffold, "basal_activity", "./cerebellum")
        self.assertEqual(empty_report.all_spikes.shape[0], int(self.simulation_duration / 0.1))
        self.assertEqual(empty_report.all_spikes.shape[1], 0)
        self.assertEqual(empty_report.nb_neurons.size, 0)
        self.assertEqual(empty_report.populations, [])

        with self.assertRaises(ValueError):
            SpikePlot((10, 10), self.scaffold, "bla", None, None, None, None, None)
        with self.assertRaises(ValueError):
            SpikePlot(
                (10, 10),
                self.scaffold,
                "basal_activity",
                None,
                None,
                None,
                nb_neurons=[1],
                populations=[],
            )
        with self.assertRaises(ValueError):
            SpikePlot(
                (10, 10),
                self.scaffold,
                "basal_activity",
                -1,
                None,
                None,
                nb_neurons=[1],
                populations=["cell"],
            )
        with self.assertRaises(ValueError):
            SpikePlot(
                (10, 10),
                self.scaffold,
                "basal_activity",
                10,
                9,
                None,
                nb_neurons=[1],
                populations=["cell"],
            )
        with self.assertRaises(ValueError):
            SpikePlot(
                (10, 10),
                self.scaffold,
                "basal_activity",
                10,
                51000,
                None,
                nb_neurons=[1],
                populations=["cell"],
            )
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = -1
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = 100000
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = -1
        with self.assertRaises(ValueError):
            self.simulationReport.time_to = 100000
        plot = SpikePlot(
            (10, 10),
            self.scaffold,
            "mf_stimulus",
            None,
            None,
            np.copy(self.simulationReport.all_spikes),
            np.copy(self.simulationReport.nb_neurons),
            self.simulationReport.populations.copy(),
        )
        plot.is_updated = True

        plot2 = ScaffoldPlot((10, 10), None)
        self.simulationReport.add_plot("simulation", plot)
        self.simulationReport.add_plot("scaffold", plot2)
        self.assertAll(self.simulationReport.all_spikes == plot.all_spikes)
        self.assertAll(self.simulationReport.nb_neurons == plot.nb_neurons)
        self.assertAll(np.array(self.simulationReport.populations) == np.array(plot.populations))
        self.assertFalse(plot.is_updated)
        self.assertEqual(self.simulationReport.simulation_name, plot.simulation_name)
        self.assertEqual(self.scaffold, plot2.scaffold)
        self.simulationReport.time_to = 500.0
        self.simulationReport.time_from = 500.0
        self.assertEqual(plot.time_to, 500.0)
        self.assertEqual(plot.time_from, 500.0)

    def test_raster_psth(self):
        plot = RasterPSTHPlot(
            (15, 10),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=self.simulationReport.all_spikes,
            nb_neurons=self.simulationReport.nb_neurons,
            populations=self.simulationReport.populations,
            dict_colors=self.simulationReport.colors,
            nb_bins=31,
        )

        plot.plot()
        self.assertEqual(np.array(plot.axes).size, len(self.simulationReport.populations) * 2)
        mf_axes = plot.get_ax()
        self.assertEqual(len(mf_axes), 2)
        xlims = np.array([self.simulationReport.time_from, self.simulationReport.time_to])
        self.assertAll(np.array(plot.get_ax()[0].get_xlim()) == xlims)
        self.assertAll(np.array(plot.get_ax()[1].get_xlim()) == xlims)
        self.assertEqual(len(plot.get_ax()[0].collections), 1)
        scatter = plot.get_ax()[0].collections[0]
        self.assertEqual(len(plot.get_ax()[1].containers), 1)
        hist = plot.get_ax()[1].containers[0]
        self.assertEqual(scatter.get_sizes()[0], 50 / self.simulationReport.nb_neurons[0])
        self.assertAll(
            np.array(scatter.get_facecolor()[0][:3])
            == self.simulationReport.colors["mossy_fibers"][:3]
        )
        self.assertEqual(scatter.get_alpha(), 1)
        self.assertTrue(scatter.get_rasterized())
        mf_spike_times = (
            np.array(
                np.where(self.simulationReport.all_spikes[:, : self.simulationReport.nb_neurons[0]])
            )
            * np.array([[self.simulationReport.dt, 1.0]]).T
        )
        self.assertAll(mf_spike_times == np.array(scatter.get_offsets()).T)
        self.assertEqual(len(hist), 30)
        self.assertEqual(hist.orientation, "vertical")

        del plot.dict_colors["mossy_fibers"]
        plot.plot(
            relative_time=True,
            params_raster={"alpha": 0.8, "edgecolors": "black", "s": 5.0},
            params_psth={"orientation": "horizontal"},
        )
        self.assertAll(np.array(plot.get_ax()[0].get_xlim()) == xlims - xlims[0])
        self.assertAll(np.array(plot.get_ax()[1].get_xlim()) == xlims - xlims[0])
        self.assertEqual(len(plot.get_ax()[0].collections), 1)
        scatter = plot.get_ax()[0].collections[0]
        self.assertEqual(len(plot.get_ax()[1].containers), 1)
        hist = plot.get_ax()[1].containers[0]
        self.assertEqual(scatter.get_sizes()[0], 5.0)
        self.assertAll(np.array(scatter.get_facecolor()[0]) == np.array([0.6, 0.6, 0.6, 0.8]))
        self.assertEqual(scatter.get_alpha(), 0.8)
        self.assertTrue(scatter.get_rasterized())
        self.assertAll(scatter.get_edgecolor()[0] == np.array([0, 0, 0, 0.8]))
        self.assertAll(mf_spike_times == np.array(scatter.get_offsets()).T)
        self.assertEqual(len(hist), 30)
        self.assertEqual(hist.orientation, "horizontal")

        # Test that an empty plot does not throw error.
        plot = RasterPSTHPlot(
            (15, 10),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=np.zeros((10001, 0), dtype=bool),
            nb_neurons=np.zeros(0, dtype=int),
            populations=[],
        )
        plot.plot()
        self.assertEqual(len(plot.get_axes()), 0)
        with self.assertRaises(ValueError):
            RasterPSTHPlot(
                (15, 10),
                scaffold=self.scaffold,
                simulation_name="basal_activity",
                time_from=None,
                time_to=None,
                all_spikes=np.zeros((10001, 0), dtype=bool),
                nb_neurons=np.zeros(0, dtype=int),
                populations=[],
                nb_bins=0,
            )

    def test_firing_rates(self):
        plot = FiringRatesPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=self.simulationReport.all_spikes,
            nb_neurons=self.simulationReport.nb_neurons,
            populations=self.simulationReport.populations,
            dict_colors=self.simulationReport.colors,
            w_single=500,
            max_neuron_sampled=100,
        )
        plot.plot()
        self.assertEqual(plot.nb_cols, 2)
        self.assertEqual(plot.nb_rows, 3)
        xlims = np.array(
            [
                self.simulationReport.time_from,
                self.simulationReport.time_to - self.simulationReport.dt,
            ]
        )
        self.assertAll(np.array(plot.firing_rates.shape) == np.array([10000, 6]))
        self.assertAll(np.absolute(np.array(plot.get_ax().get_xlim()) - xlims) <= 1e-7)
        self.assertEqual(len(plot.get_ax().lines), 1)
        self.assertAll(plot.get_ax().lines[0].get_path().vertices[:, 1] == plot.firing_rates[:, 0])
        plot.plot(relative_time=True)
        self.assertAll(np.absolute(np.array(plot.get_ax().get_xlim()) - xlims + xlims[0]) <= 1e-7)
        # Test that an empty plot does not throw error.
        plot = FiringRatesPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=np.zeros((10001, 0), dtype=bool),
            nb_neurons=np.zeros(0, dtype=int),
            populations=[],
        )
        plot.plot()
        self.assertEqual(len(plot.get_axes()), 0)

        with self.assertRaises(ValueError):
            FiringRatesPlot(
                (15, 10),
                scaffold=self.scaffold,
                simulation_name="basal_activity",
                time_from=None,
                time_to=None,
                all_spikes=np.zeros((10001, 0), dtype=bool),
                nb_neurons=np.zeros(0, dtype=int),
                populations=[],
                w_single=0,
            )
        with self.assertRaises(ValueError):
            FiringRatesPlot(
                (15, 10),
                scaffold=self.scaffold,
                simulation_name="basal_activity",
                time_from=None,
                time_to=None,
                all_spikes=np.zeros((10001, 0), dtype=bool),
                nb_neurons=np.zeros(0, dtype=int),
                populations=[],
                max_neuron_sampled=0,
            )

    def test_plot_isis(self):
        plot = ISIPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=self.simulationReport.all_spikes,
            nb_neurons=self.simulationReport.nb_neurons,
            populations=self.simulationReport.populations,
            dict_colors=self.simulationReport.colors,
            nb_bins=50,
        )
        plot.plot()
        self.assertEqual(len(plot.get_ax().containers), 1)
        hist = plot.get_ax().containers[0]
        self.assertEqual(len(hist), 50)
        self.assertEqual(hist.orientation, "vertical")

        plot.plot(orientation="horizontal")
        self.assertEqual(len(plot.get_ax().containers), 1)
        hist = plot.get_ax().containers[0]
        self.assertEqual(len(hist), 50)
        self.assertEqual(hist.orientation, "horizontal")

        # Test that an empty plot does not throw error.
        plot = ISIPlot(
            (15, 10),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=np.zeros((10001, 0), dtype=bool),
            nb_neurons=np.zeros(0, dtype=int),
            populations=[],
            dict_colors=self.simulationReport.colors,
        )
        plot.plot()
        self.assertEqual(len(plot.get_axes()), 0)
        with self.assertRaises(ValueError):
            ISIPlot(
                (15, 10),
                scaffold=self.scaffold,
                simulation_name="basal_activity",
                time_from=None,
                time_to=None,
                all_spikes=np.zeros((10001, 0), dtype=bool),
                nb_neurons=np.zeros(0, dtype=int),
                populations=[],
                dict_colors=self.simulationReport.colors,
                nb_bins=0,
            )

    def test_freq_plot(self):
        plot = FrequencyPlot(
            (15, 6),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=self.simulationReport.all_spikes,
            nb_neurons=self.simulationReport.nb_neurons,
            populations=self.simulationReport.populations,
            dict_colors=self.simulationReport.colors,
            w_single=500,
            max_neuron_sampled=100,
        )
        plot.plot()
        self.assertAll(np.array(plot.firing_rates.shape) == np.array([10000, 6]))
        self.assertAll(np.array(plot.frequencies.shape) == np.array((6, 5000)))
        self.assertAll(np.array(plot.freq_powers.shape) == np.array((6, 5000)))
        self.assertEqual(
            len(plot.get_ax().lines),
            4,
            "There should be 1 line for the freq + 3 vertical lines for bands",
        )
        self.assertEqual(plot.get_ax().lines[0].get_alpha(), None)
        self.assertAll(
            np.absolute(np.array(plot.get_ax().get_xlim()) - np.array([0, 30.0])) <= 1e-7
        )

        plot.plot(max_freq=40.0, plot_bands=False, alpha=0.7)
        self.assertEqual(len(plot.get_ax().lines), 1, "There should be 1 line for the freq")
        self.assertEqual(plot.get_ax().lines[0].get_alpha(), 0.7)
        self.assertAll(
            np.absolute(np.array(plot.get_ax().get_xlim()) - np.array([0, 40.0])) <= 1e-7
        )

        # Test that an empty plot does not throw error.
        plot = FrequencyPlot(
            (15, 10),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=np.zeros((10001, 0), dtype=bool),
            nb_neurons=np.zeros(0, dtype=int),
            populations=[],
            dict_colors=self.simulationReport.colors,
        )
        plot.plot()
        self.assertEqual(len(plot.get_axes()), 0)

    def test_sim_table(self):
        plot = SimResultsTable(
            (5, 2.5),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=self.simulationReport.all_spikes,
            nb_neurons=self.simulationReport.nb_neurons,
            populations=self.simulationReport.populations,
            dict_colors=self.simulationReport.colors,
            dict_abv={"bla": "go", "granule_cell": "GrC"},
        )
        plot.plot()
        rows = self.simulationReport.populations
        rows[rows.index("granule_cell")] = "GrC"

        self.assertAll(np.array(plot.rows) == np.array(rows))
        self.assertAll(
            np.asarray(np.array(plot.table_values).shape)
            == np.array([len(self.simulationReport.populations), 2])
        )
        self.assertAll(
            np.array(list(plot.get_firing_rates().keys()))
            == np.array(self.simulationReport.populations)
        )
        for expected, tested in zip([v[0] for v in plot._values], plot.get_firing_rates().values()):
            self.assertAll(np.array(tested) == np.array(expected))
        self.assertAll(
            np.array(list(plot.get_isis_values().keys()))
            == np.array(self.simulationReport.populations)
        )
        for expected, tested in zip([v[1] for v in plot._values], plot.get_isis_values().values()):
            self.assertAll(np.array(tested) == np.array(expected))

        # Test that an empty plot does not throw error.
        plot = SimResultsTable(
            (5, 2.5),
            scaffold=self.scaffold,
            simulation_name="basal_activity",
            time_from=None,
            time_to=None,
            all_spikes=np.zeros((10001, 0), dtype=bool),
            nb_neurons=np.zeros(0, dtype=int),
            populations=[],
            dict_colors=self.simulationReport.colors,
        )
        with self.assertWarns(UserWarning):
            plot.plot()

    def test_basic_simulation_report(self):
        report = BasicSimulationReport(self.scaffold, "basal_activity", "./")
        plot_keys = np.array(["raster_psth", "table", "firing_rates", "isis", "freq", "legend"])
        self.assertAll(np.array(list(report.plots.keys())) == plot_keys)
        filename = "test_report.pdf"
        report.print_report(filename, dpi=100)
        # should be seven cell types
        self.assertEqual(len(report.plots["table"].table_values), 6)
        # Raster PSTH plot should have two sub plots for each population
        self.assertEqual(len(report.plots["raster_psth"].get_ax()[0].collections), 1)
        self.assertEqual(len(report.plots["raster_psth"].get_ax()[1].containers), 1)
        # Firing rates plot should store firing_rates
        self.assertAll(
            np.array(report.plots["firing_rates"].firing_rates.shape) == np.array([10000, 6])
        )
        # ISIs histogram should have 50 bars
        self.assertEqual(len(report.plots["isis"].get_ax().containers[0]), 50)
        # Frequency analysis plot should store the frequencies distrib.
        self.assertAll(np.array(report.plots["freq"].frequencies.shape) == np.array((6, 5000)))
        self.assertAll(np.array(report.plots["freq"].freq_powers.shape) == np.array((6, 5000)))
        # only 6 cell types in the legend
        self.assertEqual(len(report.plots["legend"].get_ax().legend_.legend_handles), 6)
        self.assertTrue(filename in os.listdir())
        os.remove(filename)


class TestExtractISIs(unittest.TestCase):
    def test_extract_isis(self):
        spikes = np.random.random((20, 10)) >= 0.85
        isis = extract_isis(spikes, 0.1)
        enough_spikes = np.zeros(10, dtype=bool)
        u, c = np.unique(np.where(spikes)[1], return_counts=True)
        enough_spikes[u] = c >= 2
        self.assertEqual(len(isis), np.count_nonzero(enough_spikes))
        loc_spikes = spikes[:, enough_spikes]
        for i in range(len(isis)):
            self.assertTrue(
                np.absolute(isis[i] - np.mean(np.diff(np.where(loc_spikes[:, i])[0] * 0.1))) <= 1e-7
            )
