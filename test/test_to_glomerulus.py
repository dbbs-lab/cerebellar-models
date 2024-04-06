"""
Unit tests of the `ConnectomeMossyGlomerulus`, `ConnectomeGlomerulus` strategies
"""

import unittest

import numpy as np
from bsb import Scaffold
from bsb_test import (
    MPI,
    FixedPosConfigFixture,
    NetworkFixture,
    NumpyTestCase,
    RandomStorageFixture,
)


class TestConnectomeMossyGlomerulus(
    RandomStorageFixture,
    NetworkFixture,
    FixedPosConfigFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.x_length = 60
        self.z_length = 20
        self.cfg.connectivity.add(
            "x_to_glomerulus",
            dict(
                strategy="cerebellum.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
                presynaptic=dict(cell_types=["test_cell"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                x_length=self.x_length,
                z_length=self.z_length,
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)

    def test_connect(self):
        self.network.compile()
        cs = self.network.get_connectivity_set("x_to_glomerulus")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()

        cell_targets = np.array([-1, -1])
        distances = np.full(len(cell_positions), -1.0)
        self.assertEqual(len(cs), len(cell_positions), "As many connection as postsyn cell")
        for from_, to_ in cs.load_connections().as_globals():
            self.assertClose(from_[1:], cell_targets)
            self.assertClose(to_[1:], cell_targets)
            diff = np.absolute(cell_positions[from_[0]] - cell_positions[to_[0]])
            self.assertTrue(diff[0] <= self.x_length)
            self.assertTrue(diff[2] <= self.z_length)
            distances[to_[0]] = np.linalg.norm(diff)
        self.assertAll(distances >= 0, "Each postsyn cell has a connection")

    def test_distance(self):
        # Override positions
        self.chunk_size = 100.0
        pos_1 = np.array([0.5, 0.5, 0.5]) * self.chunk_size
        pos_2 = np.array([0.9, 0.5, 0.5]) * self.chunk_size  # close enough
        pos_3 = np.array([0.5, 0.5, 1.0]) * self.chunk_size  # to far away
        self.cfg.network.chunk_size = self.chunk_size
        self.cfg.cell_types["test_cell"].spatial.count = 3
        self.cfg.placement.ch4_c25.positions = MPI.bcast(np.vstack((pos_1, pos_2, pos_3)))
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)
        nb_trials = 100
        sources = np.full(nb_trials, -1)
        for i in range(nb_trials):
            self.network.compile(clear=True)
            cs = self.network.get_connectivity_set("x_to_glomerulus")
            cell_positions = self.network.get_placement_set("test_cell").load_positions()
            for from_, to_ in cs.load_connections().as_globals():
                if np.all(cell_positions[to_[0]] == pos_1):
                    if np.all(cell_positions[from_[0]] == pos_2):
                        sources[i] = 1
                    elif np.all(cell_positions[from_[0]] == pos_1):
                        sources[i] = 0
                    else:
                        sources[i] = 2
                    break
        self.assertAll(sources != 2, "pos_2 should be unreachable")
        self.assertTrue(np.count_nonzero(sources == 0) > np.count_nonzero(sources == 1))
