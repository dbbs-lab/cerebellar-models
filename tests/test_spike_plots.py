import os
import unittest
from copy import deepcopy
from os.path import abspath, dirname, join

import numpy as np
from bsb import Scaffold, parse_configuration_content
from bsb_test import NumpyTestCase, RandomStorageFixture
from matplotlib import pyplot as plt
from neo import SpikeTrain
from quantities import ms

from cerebellar_models.analysis.plots import ScaffoldPlot
from cerebellar_models.analysis.spike_plots import (
    BasicSimulationReport,
    FiringRatesPlot,
    FrequencyPlot,
    ISIPlot,
    RasterPSTHPlot,
    SimResultsTable,
    SortedPSTH,
    SpikeCorrelationPlot,
    SpikePlot,
    SpikeSimulationReport,
)
from cerebellar_models.analysis.spiking_results import extract_isis


class MiniCerebCircuitFixture(RandomStorageFixture, engine_name="hdf5", setup_cls=True):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        ROOT_FOLDER = abspath(dirname(dirname(__file__)))
        os.chdir(ROOT_FOLDER)
        # one third of the canonical circuit
        nest_folder = "configurations/mouse/in-vitro/nest/"
        dict_cfg = {
            "components": ["cerebellar_models/nest_models/build_models.py"],
            "$import": {
                "ref": "configurations/mouse/mouse_cerebellar_cortex.yaml#/",
                "values": [
                    "storage",
                    "network",
                    "regions",
                    "partitions",
                    "morphologies",
                    "cell_types",
                    "placement",
                    "connectivity",
                ],
            },
            "simulations": {
                "basal_activity": {
                    "modules": ["cerebmodule"],
                    "$import": {
                        "ref": join(nest_folder, "basal_vitro.yaml")
                        + "#/simulations/basal_activity",
                        "values": [
                            "simulator",
                            "resolution",
                            "duration",
                            "seed",
                            "cell_models",
                            "devices",
                        ],
                    },
                    "cell_models": {
                        "$import": {
                            "ref": join(nest_folder, "cell_models/eglif_cond_alpha_multisyn.yaml")
                            + "#/simulations/basal_activity/cell_models",
                            "values": [
                                "granule_cell",
                                "golgi_cell",
                                "purkinje_cell",
                                "basket_cell",
                                "stellate_cell",
                            ],
                        },
                    },
                    "connection_models": {
                        "$import": {
                            "ref": join(nest_folder, "connection_models/static_synapse.yaml")
                            + "#/simulations/basal_activity/connection_models",
                            "values": [
                                "mossy_fibers_to_glomerulus",
                                "glomerulus_to_granule",
                                "glomerulus_to_golgi",
                                "golgi_to_glomerulus",
                                "golgi_to_golgi",
                                "ascending_axon_to_golgi",
                                "parallel_fiber_to_golgi",
                                "parallel_fiber_to_purkinje",
                                "ascending_axon_to_purkinje",
                                "parallel_fiber_to_stellate",
                                "stellate_to_stellate",
                                "stellate_to_purkinje",
                                "parallel_fiber_to_basket",
                                "basket_to_basket",
                                "basket_to_purkinje",
                            ],
                        },
                    },
                    "devices": {
                        "$import": {
                            "ref": join(nest_folder, "cell_models/eglif_cond_alpha_multisyn.yaml")
                            + "#/simulations/basal_activity/devices",
                            "values": [
                                "granule_record",
                                "golgi_record",
                                "purkinje_record",
                                "basket_record",
                                "stellate_record",
                            ],
                        }
                    },
                },
                "mf_stimulus": {
                    "$import": {
                        "ref": "#/simulations/basal_activity",
                        "values": [
                            "simulator",
                            "resolution",
                            "duration",
                            "modules",
                            "seed",
                            "cell_models",
                            "connection_models",
                            "devices",
                        ],
                    },
                    "devices": {
                        "stimulus": {
                            "device": "poisson_generator",
                            "rate": 150,
                            "start": 1200,
                            "stop": 1250,
                            "targetting": {
                                "strategy": "sphere",
                                "radius": 90,
                                "origin": [150.0, 65.0, 100.0],
                                "cell_models": ["mossy_fibers"],
                            },
                            "weight": 1.0,
                            "delay": 0.1,
                        }
                    },
                },
            },
        }
        cls.cfg = parse_configuration_content(dict_cfg, path=os.path.abspath("./config.json"))
        cls.cfg.network.x = 100
        cls.cfg.network.y = 66
        cls.cfg.network.z = 100
        cls.cfg.partitions.granular_layer.thickness = 40
        cls.cfg.partitions.purkinje_layer.thickness = 10
        cls.cfg.partitions.b_molecular_layer.thickness = 17
        cls.cfg.partitions.t_molecular_layer.thickness = 33
        # make sure there are enough mfs.
        cls.cfg.cell_types.glomerulus.spatial.density = 0.00034

        cls.scaffold = Scaffold(cls.cfg, cls.storage)
        cls.scaffold.compile(skip_after_connectivity=True, clear=True)


class ReportBasalSimCircuitFixture(MiniCerebCircuitFixture, engine_name="hdf5", setup_cls=True):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.simulation_duration = 1000.0
        cls.scaffold.simulations["basal_activity"].duration = cls.simulation_duration
        cls.scaffold.simulations["mf_stimulus"].duration = cls.simulation_duration
        cls.simulation_results = cls.scaffold.run_simulation("basal_activity")
        cls.simulation_results.write("test_sim_results.nio", "ow")
        cls.simulationReport = SpikeSimulationReport(
            cls.scaffold,
            "basal_activity",
            "./",
        )

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        plt.close("all")
        os.remove("test_sim_results.nio")


class TestSpikePlots(
    ReportBasalSimCircuitFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
    setup_cls=True,
):
    def test_spike_reports(self):
        self.assertEqual(self.simulationReport.time_to, self.simulation_duration)
        self.assertEqual(self.simulationReport.dt, 0.1)
        self.assertEqual(
            len(self.simulationReport.nb_neurons), len(self.simulationReport.populations)
        )
        self.assertEqual(
            len(self.simulationReport.filt_spikes), len(self.simulationReport.nb_neurons)
        )
        self.assertAll(
            np.array([len(st.magnitude) for st in self.simulationReport.filt_spikes]) > 0
        )
        self.assertTrue("mossy_fibers" in self.simulationReport.populations)
        self.assertTrue("glomerulus" not in self.simulationReport.populations)
        with self.assertRaises(ValueError):
            SpikeSimulationReport(self.scaffold, "blabla", "./")

        empty_report = SpikeSimulationReport(self.scaffold, "basal_activity", "./cerebellar_models")
        self.assertEqual(len(empty_report.filt_spikes), 0)
        self.assertEqual(empty_report.nb_neurons.size, 0)
        self.assertEqual(empty_report.populations, [])
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = -1
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = 100000
        with self.assertRaises(ValueError):
            self.simulationReport.time_from = -1
        with self.assertRaises(ValueError):
            self.simulationReport.time_to = 100000

    def test_update_spike_report(self):
        plot = SpikePlot(
            (10, 10),
            self.simulationReport.spiking_results,
            None,
        )
        plot.is_updated = True

        plot2 = ScaffoldPlot((10, 10), None)
        self.simulationReport.add_plot("simulation", plot)
        self.simulationReport.add_plot("scaffold", plot2)
        self.assertAll(
            np.array(
                [
                    np.all(s1 == s2)
                    for s1, s2 in zip(self.simulationReport.filt_spikes, plot.filt_spikes)
                ]
            )
        )
        self.assertAll(self.simulationReport.nb_neurons == plot.nb_neurons)
        self.assertAll(np.array(self.simulationReport.populations) == np.array(plot.populations))
        self.assertFalse(plot.is_updated)
        self.assertEqual(
            self.simulationReport.spiking_results.simulation_name,
            plot.spiking_results.simulation_name,
        )
        self.assertEqual(self.scaffold, plot2.scaffold)
        self.simulationReport.time_to = 500.0
        self.simulationReport.time_from = 500.0
        self.assertEqual(plot.time_to, 500.0)
        self.assertEqual(plot.time_from, 500.0)

    def test_raster_psth(self):
        plot = RasterPSTHPlot(
            (15, 10),
            spiking_results=self.simulationReport.spiking_results,
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
        mf_spikes = self.simulationReport.filt_spikes[0]
        mf_spike_times = (
            np.array(
                (
                    mf_spikes.magnitude / self.simulationReport.dt,
                    np.unique(mf_spikes.array_annotations["senders"], return_inverse=True)[1],
                )
            )
            * np.array([[self.simulationReport.dt, 1.0]]).T
        )
        self.assertAll(np.absolute(mf_spike_times - np.array(scatter.get_offsets()).T) <= 1e-7)
        self.assertEqual(len(hist), 30)
        self.assertEqual(hist.orientation, "vertical")

    def test_relative_time(self):
        plot = RasterPSTHPlot(
            (15, 10),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
            nb_bins=31,
        )
        xlims = np.array([self.simulationReport.time_from, self.simulationReport.time_to])
        mf_spikes = self.simulationReport.filt_spikes[0]
        mf_spike_times = (
            np.array(
                (
                    mf_spikes.magnitude / self.simulationReport.dt,
                    np.unique(mf_spikes.array_annotations["senders"], return_inverse=True)[1],
                )
            )
            * np.array([[self.simulationReport.dt, 1.0]]).T
        )
        # test also if absence of color
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
        self.assertAll(np.absolute(mf_spike_times - np.array(scatter.get_offsets()).T) <= 1e-7)
        self.assertEqual(len(hist), 30)
        self.assertEqual(hist.orientation, "horizontal")

    def test_raster_psth_subinterval(self):
        # check for sub interval
        xlims = np.array([200.0, 800.0])
        plot = RasterPSTHPlot(
            (15, 10),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
            nb_bins=31,
        )
        old_time_from = plot.time_from
        old_time_to = plot.time_to
        plot.time_from = xlims[0]
        plot.time_to = xlims[1]
        plot.plot()
        loc_mf_spikes = self.simulationReport.filt_spikes[0].time_slice(xlims[0], xlims[1])
        mf_spike_times = (
            np.array(
                (
                    loc_mf_spikes.magnitude / self.simulationReport.dt,
                    np.unique(loc_mf_spikes.array_annotations["senders"], return_inverse=True)[1],
                )
            )
            * np.array([[self.simulationReport.dt, 1.0]]).T
        )
        self.assertAll(np.array(plot.get_ax()[0].get_xlim()) == xlims)
        self.assertAll(np.array(plot.get_ax()[1].get_xlim()) == xlims)
        self.assertEqual(len(plot.get_ax()[0].collections), 1)
        scatter = plot.get_ax()[0].collections[0]
        self.assertAll(np.absolute(mf_spike_times - np.array(scatter.get_offsets()).T) <= 1e-7)
        self.assertEqual(len(plot.get_ax()[1].containers), 1)
        hist = plot.get_ax()[1].containers[0]
        self.assertEqual(len(hist), 30)
        plot.clear()
        self.assertEqual(len(plot.get_ax()[0].collections), 0)
        self.assertEqual(len(plot.get_ax()[0].containers), 0)
        plot.time_from = old_time_from
        plot.time_to = old_time_to

    def test_raster_psth_bins_error(self):
        with self.assertRaises(ValueError):
            RasterPSTHPlot(
                (15, 10),
                spiking_results=self.simulationReport.spiking_results,
                nb_bins=0,
            )

    def test_firing_rates(self):
        plot = FiringRatesPlot(
            (15, 6),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
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

    def test_firing_rates_error_kernel(self):
        with self.assertRaises(TypeError):
            FiringRatesPlot(
                (15, 10),
                spiking_results=self.simulationReport.spiking_results,
                kernel=0,
            )

    def test_plot_isis(self):
        plot = ISIPlot(
            (15, 6),
            spiking_results=self.simulationReport.spiking_results,
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

    def test_plot_isis_error_bins(self):
        with self.assertRaises(ValueError):
            ISIPlot(
                (15, 10),
                spiking_results=self.simulationReport.spiking_results,
                dict_colors=self.simulationReport.colors,
                nb_bins=0,
            )

    def test_freq_plot(self):
        plot = FrequencyPlot(
            (15, 6),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
        )
        plot.plot()
        self.assertAll(np.array(plot.firing_rates.shape) == np.array([10000, 6]))
        self.assertAll(np.array(plot.frequencies.shape) == np.array((6, 5000)))
        self.assertAll(np.array(plot.freq_powers.shape) == np.array((6, 5000)))
        self.assertEqual(
            len(plot.get_ax().lines),
            5,
            "There should be 1 line for the freq + 4 vertical lines for bands",
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

    def test_sim_table(self):
        plot = SimResultsTable(
            (5, 2.5),
            spiking_results=self.simulationReport.spiking_results,
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
        self.assertAll(np.array(list(plot.get_firing_rates().keys())) == np.array(rows))
        for expected, tested in zip([v[0] for v in plot._values], plot.get_firing_rates().values()):
            self.assertAll(np.array(tested) == np.array(expected))
        self.assertAll(np.array(list(plot.get_isis_values().keys())) == np.array(rows))
        for expected, tested in zip([v[1] for v in plot._values], plot.get_isis_values().values()):
            self.assertAll(np.array(tested) == np.array(expected))

    def test_corr_matrix(self):
        plot = SpikeCorrelationPlot(
            (10, 10.5),
            spiking_results=self.simulationReport.spiking_results,
            dict_abv=self.simulationReport.abbreviations,
        )
        plot.plot()
        self.assertEqual(plot.corrcoef.shape, (6, 6))
        self.assertAll(plot.corrcoef <= 1)
        self.assertAll(plot.corrcoef >= -1)

    def test_sorted_psth_defaults(self):
        plot = SortedPSTH(
            (10, 6),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
        )
        self.assertEqual(plot.nb_bins, 50)
        self.assertEqual(plot.sample_size, 100)

    def test_sorted_psth_bins_error(self):
        with self.assertRaises(ValueError):
            SortedPSTH(
                (10, 6),
                spiking_results=self.simulationReport.spiking_results,
                nb_bins=0,
            )

    def test_sorted_psth(self):
        nb_bins = 20
        plot = SortedPSTH(
            (10, 6),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
            nb_bins=nb_bins,
            sample_size=1,
        )
        self.assertEqual(plot.nb_bins, nb_bins)
        self.assertEqual(plot.sample_size, 1)
        plot.plot()
        ax = plot.get_ax()
        self.assertEqual(ax.get_xlabel(), "Time in ms")
        self.assertEqual(ax.get_ylabel(), "Sorted Neuron id")
        # Exactly one image should have been drawn by imshow
        self.assertEqual(len(ax.images), 1)
        # With sample_size=1 and every population spiking, the image has one
        # row per population and nb_bins columns over the RGBA channels.
        img = np.asarray(ax.images[0].get_array())
        self.assertEqual(img.shape, (len(plot.populations), nb_bins, 4))
        # alpha channel is normalized to its max along the time axis
        self.assertAll(np.nanmax(img[..., 3], axis=1) <= 1 + 1e-7)
        # y-ticks: one per population
        self.assertAll(np.array(ax.get_yticks()) == np.arange(len(plot.populations)))
        self.assertEqual(len(ax.get_yticklabels()), len(plot.populations))
        # the y-tick labels are the populations re-ordered by argsort(max_delta_times)
        self.assertEqual(
            sorted(t.get_text() for t in ax.get_yticklabels()),
            sorted(plot.populations),
        )
        # x-ticks: every other bin edge over [0, nb_bins)
        self.assertAll(np.array(ax.get_xticks()) == np.arange(0, nb_bins, 2))

    def test_sorted_psth_missing_color(self):
        plot = SortedPSTH(
            (10, 6),
            spiking_results=self.simulationReport.spiking_results,
            dict_colors=self.simulationReport.colors,
            nb_bins=10,
            sample_size=1,
        )
        # remove one population's color so it falls back to the default gray
        del plot.dict_colors["mossy_fibers"]
        plot.plot()
        img = np.asarray(plot.get_ax().images[0].get_array())
        # The first time-bin RGB values per row encode the cell-type color.
        # Exactly one row (the mossy_fibers neuron) should fall back to gray
        # [0.6, 0.6, 0.6] — every other row must still carry its cell-type color.
        rgb_per_row = img[:, 0, :3]
        is_gray = np.all(np.absolute(rgb_per_row - 0.6) <= 1e-7, axis=-1)
        self.assertEqual(int(np.count_nonzero(is_gray)), 1)

    def test_basic_simulation_report(self):
        report = BasicSimulationReport(self.scaffold, "basal_activity", "./")
        plot_keys = np.array(
            ["raster_psth", "table", "firing_rates", "isis", "freq", "corr", "legend"]
        )
        self.assertAll(np.array(list(report.plots.keys())) == plot_keys)
        filename = "test_report.pdf"
        report.print_report(filename, dpi=100)
        # should be seven cell types
        self.assertEqual(len(report.plots["table"].table_values), 6)
        # Raster PSTH plot should have two sub-plots for each population
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
