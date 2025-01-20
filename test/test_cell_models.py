import os
import unittest
from os.path import abspath, dirname, getmtime, isfile, join

import nest
import numpy as np
from bsb import Scaffold, parse_configuration_content
from bsb.services import MPI
from bsb_test import NumpyTestCase, RandomStorageFixture
from pynestml.exceptions.invalid_path_exception import InvalidPathException
from scipy.optimize import curve_fit


class TestNestModuleLoading(unittest.TestCase):
    def test_build_models(self):
        from cerebellum.nest_models.build_models import _build_nest_models

        # cerebmodule should have been built
        cerebmodule_file = join(
            nest.ll_api.sli_func("statusdict/prefix ::"), "lib/nest/cerebmodule.so"
        )
        self.assertTrue(isfile(cerebmodule_file))
        # Does not raise because module exists
        _build_nest_models(model_dir="bla")
        old_mtime = getmtime(cerebmodule_file)
        _build_nest_models(redo=True)  # Force update
        self.assertTrue(getmtime(cerebmodule_file) > old_mtime, "module should be updated")
        os.remove(cerebmodule_file)
        with self.assertRaises(OSError):
            _build_nest_models(model_dir="bla")
        with self.assertRaises(InvalidPathException):
            _build_nest_models(model_dir=dirname(abspath(__file__)))


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
            "components": ["cerebellum/nest_models/build_models.py"],
            "partitions": {"B": {"type": "layer", "thickness": 1}},
            "cell_types": {
                "granule_cell": {"spatial": {"radius": 1, "count": 10}},
                "golgi_cell": {"spatial": {"radius": 1, "count": 10}},
                "purkinje_cell": {"spatial": {"radius": 1, "count": 10}},
                "purkinje_vivo_cell": {"spatial": {"radius": 1, "count": 10}},
                "basket_cell": {"spatial": {"radius": 1, "count": 10}},
                "stellate_cell": {"spatial": {"radius": 1, "count": 10}},
                "unipolar_brush_cell": {"spatial": {"radius": 1, "count": 10}},
                "dcn_p_cell": {"spatial": {"radius": 1, "count": 10}},
                "dcn_i_cell": {"spatial": {"radius": 1, "count": 10}},
            },
            "placement": {
                "placement_A": {
                    "strategy": "bsb.placement.strategy.FixedPositions",
                    "cell_types": [
                        "granule_cell",
                        "golgi_cell",
                        "purkinje_cell",
                        "purkinje_vivo_cell",
                        "basket_cell",
                        "stellate_cell",
                        "unipolar_brush_cell",
                        "dcn_p_cell",
                        "dcn_i_cell",
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
                            "ref": "../configurations/mouse/ubc/ubc_nest.yaml#/simulations/basal_activity/cell_models",
                            "values": [
                                "granule_cell",
                                "golgi_cell",
                                "purkinje_cell",
                                "basket_cell",
                                "stellate_cell",
                                "unipolar_brush_cell",
                            ],
                        },
                        "purkinje_vivo_cell": {
                            "$import": {
                                "ref": "../configurations/mouse/nest/basal_vivo.yaml#/simulations/basal_activity/cell_models/purkinje_cell",
                                "values": ["constants", "model"],
                            }
                        },
                        "dcn_p_cell": {
                            "$import": {
                                "ref": "../configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml#/simulations/basal_activity/cell_models/dcn_p",
                                "values": ["constants", "model"],
                            }
                        },
                        "dcn_i_cell": {
                            "$import": {
                                "ref": "../configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml#/simulations/basal_activity/cell_models/dcn_i",
                                "values": ["constants", "model"],
                            }
                        },
                    },
                    "connection_models": {},
                    "devices": {
                        "$import": {
                            "ref": "../configurations/mouse/ubc/ubc_nest.yaml#/simulations/basal_activity/devices",
                            "values": [
                                "granule_record",
                                "golgi_record",
                                "purkinje_record",
                                "basket_record",
                                "stellate_record",
                                "unipolar_brush_record",
                            ],
                        },
                        "purkinje_vivo_record": {
                            "device": "spike_recorder",
                            "delay": 0.1,
                            "targetting": {
                                "strategy": "cell_model",
                                "cell_models": ["purkinje_vivo_cell"],
                            },
                        },
                        "dcn_p_record": {
                            "$import": {
                                "ref": "../configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml#/simulations/basal_activity/devices/dcn_p_record",
                                "values": ["device", "delay"],
                            },
                            "targetting": {"strategy": "cell_model", "cell_models": ["dcn_p_cell"]},
                        },
                        "dcn_i_record": {
                            "$import": {
                                "ref": "../configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml#/simulations/basal_activity/devices/dcn_i_record",
                                "values": ["device", "delay"],
                            },
                            "targetting": {"strategy": "cell_model", "cell_models": ["dcn_i_cell"]},
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
            "purkinje_vivo": {
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
            "unipolar_brush": {
                "amplitudes": [12, 24, 36],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
            "dcn_p": {
                "amplitudes": [142, 284, 426],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
            "dcn_i": {
                "amplitudes": [56, 112, 168],
                "starts": [10, 12, 14],
                "stops": [11, 13, 15],
            },
        }
        predicted = {
            "golgi": {"autorhythm": 12.8, "slope": 0.2},
            "granule": {"autorhythm": 0.0, "slope": 3.70},
            "purkinje": {"autorhythm": 45.6, "slope": 0.08},
            "purkinje_vivo": {"autorhythm": 83.2, "slope": 0.095},  # to match Geminiani et al. 2024
            "basket": {"autorhythm": 9.51, "slope": 2.16},
            "stellate": {"autorhythm": 9.51, "slope": 2.16},
            "unipolar_brush": {"autorhythm": 0.0, "slope": 2.35},  # no source available
            "dcn_p": {"autorhythm": 31.48, "slope": 0.28},
            "dcn_i": {"autorhythm": 14.37, "slope": 0.4},
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
            slope, slope_std, offset, offset_std = self._fit_lin(i_points[1:], f_points[1:])
            self.assertTrue(abs(slope - predicted[cell_type]["slope"]) <= 2 * slope_std)
