import os
import unittest

import numpy as np
from bsb import Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture

from cerebellum.analysis.plots import ScaffoldPlot
from cerebellum.analysis.spiking_results import SpikePlot, SpikeSimulationReport


class MiniCerebCircuitTest(RandomStorageFixture, unittest.TestCase, engine_name="hdf5"):
    def setUp(self):
        super().setUp()
        # one third of the canonical circuit
        self.cfg = parse_configuration_file("configurations/mouse/nest/stimulus_mossy.yaml")
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


class TestSpikePlotReport(MiniCerebCircuitTest, NumpyTestCase, engine_name="hdf5"):
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
        os.remove("test_sim_results.nio")

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
        self.assertEqual(empty_report.all_spikes.shape[0], int(self.simulation_duration / 0.1) + 1)
        self.assertEqual(empty_report.all_spikes.shape[1], 0)
        self.assertEqual(empty_report.nb_neurons.size, 0)
        self.assertEqual(empty_report.populations, [])

        with self.assertRaises(ValueError):
            SpikePlot((10, 10), self.scaffold, "bla", None, None, None, None, None)

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
