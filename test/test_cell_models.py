import pathlib
import unittest

import numpy as np
from bsb import Scaffold, parse_configuration_file
from scipy.signal import find_peaks


@unittest.skip("Needs to be updated to bsb v4.0")
class TestSingleCellModels(unittest.TestCase):
    # fixme
    @unittest.skip(reason="13Hz instead of 8")
    def test_golgi_autorythm(self):
        # Read the config and build a network with a single Golgi cell
        config = parse_configuration_file(
            pathlib.Path(__file__).parent / "test_configs" / "test_nrn_goc_autorythm.json", "json"
        )
        scaffold = Scaffold(config, self.storage)
        scaffold.compile()
        # Run a simulation without any stimulus to the Golgi cell.
        # The expected result is the autorythm of the Golgi cell at 8 pm 1 Hz
        result = scaffold.run_simulation("neurontest_goc_test")
        simulation_time = float(config.simulations.neurontest_goc_test.duration)
        resolution = float(config.simulations.neurontest_goc_test.resolution)
        time = np.arange(200, simulation_time + resolution, resolution)
        # fixme: brittle retrieval of results: filter result using device name or smth
        avg = np.mean(result.block.segments[0].analogsignals[0], axis=1)
        after_transient = int(200.0 / resolution)
        voltage = avg[after_transient:]
        peaks, _ = find_peaks(voltage, height=0)
        self.assertAlmostEqual(len(peaks), 8, msg="The expected result is 8 pm 1 Hz", delta=1)

    # fixme
    @unittest.skip(reason="0Hz instead of 40")
    def test_golgi_current_clamp_200_pa(self):
        # Build a network with a single Golgi cell
        config = parse_configuration_file(
            pathlib.Path(__file__).parent
            / "test_configs"
            / "test_nrn_goc_current_clamp_200_pa.json",
            "json",
        )
        scaffold = Scaffold(config, self.storage)
        scaffold.compile()
        # Run a simulation stimulating the Golgi cell with a 200 pA current.
        # The expected result is a burst at 40 pm 2 Hz
        result = scaffold.run_simulation("neurontest_goc_test")
        simulation_time = float(config.simulations.neurontest_goc_test.duration)
        resolution = float(config.simulations.neurontest_goc_test.resolution)
        time = np.arange(200, simulation_time + resolution, resolution)
        # fixme: brittle retrieval of results
        avg = np.mean(result.block.segments[0].analogsignals[1], axis=1)
        after_transient = int(200.0 / resolution)
        voltage = avg[after_transient:]
        peaks, _ = find_peaks(voltage, height=0)
        self.assertAlmostEqual(len(peaks), 40, msg="The expected result is 65 pm 7 Hz", delta=2)

    @unittest.skip(reason="Test too time consuming")
    def test_pc_autorythm(self):
        # Build a network with a single Purkinje cell
        config = parse_configuration_file(
            pathlib.Path(__file__).parent / "test_configs" / "test_nrn_pc_autorthythm.json", "json"
        )
        scaffold = Scaffold(config, self.storage)
        scaffold.compile()
        # todo: find more exact range of frequency
        # Run a simulation without any stimulus to the Purkinje cell.
        # The expected result is the autorythm of the Purkinje cell (17-148 Hz)
        # (Raman IM, Bean BP. Ionic currents underlying spontaneous action potentials
        # in isolated cerebellar Purkinje neurons. J Neurosci. 1999 )
        result = scaffold.run_simulation("neurontest_pc_test")
        simulation_time = float(config.simulations.neurontest_pc_test.duration)
        resolution = float(config.simulations.neurontest_pc_test.resolution)
        time = np.arange(200, simulation_time + resolution, resolution)
        # fixme: brittle retrieval of results
        avg = np.mean(result.block.segments[0].analogsignals[0], axis=1)
        after_transient = int(200.0 / resolution)
        voltage = avg[after_transient:]
        peaks, _ = find_peaks(voltage, height=0)
        self.assertGreaterEqual(len(peaks), 17)
        self.assertLessEqual(len(peaks), 148)

    @unittest.skip(reason="Test too time consuming")
    def test_granule_purkinje(self):
        # Build a network with a single Purkinje cell and ~ 70 GrCs connected to the Purkinje
        config = parse_configuration_file(
            pathlib.Path(__file__).parent / "test_configs" / "test_nrn_grc_pc.json", "json"
        )
        scaffold = Scaffold(config, self.storage)
        scaffold.compile()
        # Stimulate GrCs with a baseline of 20 pA and a 25 pA current starting at 70 ms.
        # The expected result is a spike at ~ 80 ms
        result = scaffold.run_simulation("neurontest_grc_pc_test")
        simulation_time = float(config.simulations.neurontest_grc_pc_test.duration)
        resolution = float(config.simulations.neurontest_grc_pc_test.resolution)
        time = np.arange(68, simulation_time + resolution, resolution)
        # fixme: brittle retrieval of results
        avg = np.mean(result.block.segments[0].analogsignals[-1], axis=1)
        after_transient = int(68.0 / resolution)
        voltage = avg[after_transient:]
        peaks, _ = find_peaks(voltage, height=0)
        self.assertEqual(len(peaks), 1, msg="The expected result is one peak")
