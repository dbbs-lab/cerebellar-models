import os
import unittest
from os.path import abspath, dirname, join

import numpy as np
from bsb import Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture
from scipy.stats import norm


class TestDuplicateSynapses(
    RandomStorageFixture, NumpyTestCase, unittest.TestCase, engine_name="hdf5"
):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        ROOT_FOLDER = abspath(dirname(dirname(__file__)))
        os.chdir(ROOT_FOLDER)

    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file(
            join(dirname(__file__), "test_configurations/canonical_mouse_awake_io_nest.json")
        )
        del self.cfg.after_connectivity["print_structure_report"]
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold.compile(
            only=[
                "io_layer_placement",
                "purkinje_layer_placement",
                "io_to_purkinje",
            ],
            skip_after_connectivity=True,
            skip_after_placement=False,
        )

    def test_duplicate_synapses(self):
        origin_conns = np.array(
            list(
                self.scaffold.get_connectivity_set("io_to_purkinje").load_connections().as_globals()
            )
        )
        u_conns = np.unique(np.concatenate(np.moveaxis(origin_conns, 0, 1), axis=-1), axis=0)
        self.scaffold.run_after_connectivity(pipelines=False)
        new_conns = np.array(
            list(
                self.scaffold.get_connectivity_set("io_to_purkinje").load_connections().as_globals()
            )
        )
        u_new_conns = np.unique(np.concatenate(np.moveaxis(new_conns, 0, 1), axis=-1), axis=0)
        self.assertAll(u_conns == u_new_conns)
        distrib_params = self.cfg.after_connectivity["duplicate_conn_io_pc"].contacts.parameters
        # Apply central limit theorem: S = sum of n i.i.d. draws X_i ~ N(loc, scale^2).
        # E[S] = n*loc, Std[S] = scale*sqrt(n)  =>  z = (S - n*loc) / (scale * sqrt(n))
        # Threshold 3.89 gives rejection probability 0.0001 (fails once in 10000 trials).
        self.assertLess(
            np.abs(new_conns.shape[0] - origin_conns.shape[0] * distrib_params["loc"])
            / distrib_params["scale"]
            / np.sqrt(origin_conns.shape[0]),
            3.89,
            "This test should fail only once in every 10000 trials",
        )
