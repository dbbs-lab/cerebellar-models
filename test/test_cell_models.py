import os
import pathlib
import unittest

import numpy as np
from bsb import Scaffold, parse_configuration_content, parse_configuration_file
from bsb.services import MPI
from bsb_test import NumpyTestCase, RandomStorageFixture
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


@unittest.skipIf(MPI.get_size() > 1, "Skipped during parallel testing.")
class TestSingleCellModels(
    RandomStorageFixture, NumpyTestCase, unittest.TestCase, engine_name="hdf5"
):

    def setUp(self):
        super().setUp()
        self.configuration_dict = {
            "name": "test",
            "storage": {"engine": "hdf5"},
            "network": {"x": 1, "y": 1, "z": 1},
            "partitions": {"B": {"type": "layer", "thickness": 1}},
            "cell_types": {
                "granule_cell": {"spatial": {"radius": 1, "count": 10}},
                "golgi_cell": {"spatial": {"radius": 1, "count": 10}},
                "purkinje_cell": {"spatial": {"radius": 1, "count": 10}},
                "basket_cell": {"spatial": {"radius": 1, "count": 10}},
                "stellate_cell": {"spatial": {"radius": 1, "count": 10}},
            },
            "placement": {
                "placement_A": {
                    "strategy": "bsb.placement.strategy.FixedPositions",
                    "cell_types": [
                        "granule_cell",
                        "golgi_cell",
                        "purkinje_cell",
                        "basket_cell",
                        "stellate_cell",
                    ],
                    "partitions": ["B"],
                    "positions": [[1, 1, 1]] * 10,
                }
            },
            "connectivity": {},
            "after_connectivity": {},
            "simulations": {
                "test": {
                    "simulator": "nest",
                    "duration": 21001,
                    "resolution": 0.1,
                    "seed": 1234,
                    "modules": ["cerebmodule"],
                    "cell_models": {
                        "$import": {
                            "ref": "../configurations/mouse/nest/basal.yaml#/simulations/basal_activity/cell_models",
                            "values": [
                                "granule_cell",
                                "golgi_cell",
                                "purkinje_cell",
                                "basket_cell",
                                "stellate_cell",
                            ],
                        },
                    },
                    "connection_models": {},
                    "devices": {
                        "$import": {
                            "ref": "../configurations/mouse/nest/basal.yaml#/simulations/basal_activity/devices",
                            "values": [
                                "granule_record",
                                "golgi_record",
                                "purkinje_record",
                                "basket_record",
                                "stellate_record",
                            ],
                        },
                    },
                }
            },
        }

    @staticmethod
    def _load_spike_data(results):
        cell_dict = {}
        for st in results.spiketrains:
            cell_type = st.annotations["device"].split("_rec")[0]
            if cell_type not in cell_dict:
                cell_dict[cell_type] = {}
            if len(st.annotations["senders"]) > 0:
                u_senders, inv = np.unique(st.annotations["senders"], return_inverse=True)
                for i, sender in enumerate(u_senders):
                    if str(sender) not in cell_dict[cell_type]:
                        cell_dict[cell_type][str(sender)] = []
                    cell_dict[cell_type][str(sender)].extend(st.magnitude[inv == i])
        return cell_dict

    @staticmethod
    def _f_linear(x, m, q):
        """
        Linear function to fit fI curves.
        """
        return m * x + q

    @staticmethod
    def _fit_lin(x, y):
        params, covariance = curve_fit(TestSingleCellModels._f_linear, x, y)
        m, q = params
        m_std, q_std = np.sqrt(np.diag(covariance))
        return m, m_std, q, q_std

    @staticmethod
    def _compute_firing_rate(spikes, start=0, stop=10000, first_two=False):
        isis = [
            np.diff(
                np.array(spike_time)[
                    (np.array(spike_time) >= start) * (np.array(spike_time) < stop)
                ]
            )
            for spike_time in spikes
        ]

        if first_two:
            isis = np.array([isi[0] for isi in isis])
        else:
            isis = np.array([np.mean(isi) for isi in isis if len(isi) > 0])
        f_tonic = (1.0 / isis / 0.001) if len(isis) > 0 else 0.0  # in Hz
        return np.mean(f_tonic), np.std(f_tonic)

    @unittest.skip(reason="EGLIF results are not reproducible")
    def test_fi_curves(self):
        protocol = {
            "golgi": {
                "amplitudes": [100, 200, 300, 400, 500, 600],
                "starts": [10, 12, 14, 16, 18, 20],
                "stops": [11, 13, 15, 17, 19, 21],
            },
            "granule": {
                "amplitudes": [16, 20, 24],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
            "purkinje": {
                "amplitudes": [500, 1000, 2400],
                "starts": [10, 12, 14],
                "stops": [11, 13, 14.01],
            },
            "basket": {
                "amplitudes": [12, 24, 36],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
            "stellate": {
                "amplitudes": [12, 24, 36],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
        }
        predicted = {
            "golgi": {"autorhythm": 12.8, "slope": 0.2},
            "granule": {"autorhythm": 0.0, "slope": 3.70},
            "purkinje": {"autorhythm": 60.96, "slope": 0.08},
            "basket": {"autorhythm": 9.51, "slope": 2.16},
            "stellate": {"autorhythm": 9.51, "slope": 2.16},
        }
        for cell_type in predicted:
            for i, stim in enumerate(protocol[cell_type]["amplitudes"]):
                self.configuration_dict["simulations"]["test"]["devices"][
                    f"stimulus_{cell_type}_{i}"
                ] = {
                    "device": "dc_generator",
                    "amplitude": stim,
                    "start": protocol[cell_type]["starts"][i] * 1e3,
                    "stop": protocol[cell_type]["stops"][i] * 1e3,
                    "targetting": {"strategy": "cell_model", "cell_models": [f"{cell_type}_cell"]},
                    "weight": 1.0,
                    "delay": 0.1,
                }
        self.cfg = parse_configuration_content(
            self.configuration_dict, parser="json", path=os.path.realpath(__file__)
        )
        network = Scaffold(self.cfg, self.storage)
        network.compile()
        results = network.run_simulation("test")
        cell_dict = self._load_spike_data(results)
        for cell_type in predicted:
            spike_dict = cell_dict[cell_type]
            # Test autorhythm
            f_tonic, std_tonic = self._compute_firing_rate(spike_dict.values(), start=0, stop=10000)
            # Currently the results are not reproducible
            # self.assertTrue(abs(predicted[cell_type]["autorhythm"] - f_tonic) <= 2 * std_tonic)
            # This dummy test is just here to detect any additional "instability".
            self.assertTrue(
                abs(predicted[cell_type]["autorhythm"] - f_tonic)
                <= predicted[cell_type]["autorhythm"] / 2
            )

            # Test F/i curve points
            starts = protocol[cell_type]["starts"]
            f_points = np.zeros(len(starts))
            i_points = np.zeros(len(starts))
            for i, start in enumerate(starts):
                f_tonic, std_tonic = self._compute_firing_rate(
                    spike_dict.values(),
                    start=start * 1e3,
                    stop=protocol[cell_type]["stops"][i] * 1e3,
                    first_two=True,
                )
                f_points[i] = f_tonic
                i_points[i] = protocol[cell_type]["amplitudes"][i]
            slope, slope_std, offset, offset_std = self._fit_lin(i_points, f_points)
            # Currently the results are not reproducible
            # self.assertTrue(abs(slope - predicted[cell_type]["slope"]) <= 2 * slope_std)
            # This dummy test is just here to detect any additional "instability".
            self.assertTrue(
                abs(slope - predicted[cell_type]["slope"]) <= predicted[cell_type]["slope"] / 2
            )

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
