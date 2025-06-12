import unittest

import numpy as np
from bsb import Configuration, ConfigurationError, Scaffold
from bsb_test import NetworkFixture, NumpyTestCase, RandomStorageFixture


class TestGlomerulus_to_UBC(
    RandomStorageFixture, NetworkFixture, NumpyTestCase, unittest.TestCase, engine_name="hdf5"
):
    def setUp(self):
        super().setUp()
        # radius is greater than 1 chunk dimensions but less than two
        self.radius = 40
        self.chunk_size = np.array([30, 30, 30])
        self.cfg = Configuration.default(
            network=dict(
                chunk_size=self.chunk_size,
                x=self.chunk_size[0] * 2,
                y=self.chunk_size[1] * 2,
                z=self.chunk_size[2] * 2,
            ),
            cell_types=dict(
                pre_cell=dict(spatial=dict(radius=2, count=80)),
                pre_cell2=dict(spatial=dict(radius=2, count=80)),
                single_cell=dict(spatial=dict(radius=2, count=8)),
                test_cell=dict(spatial=dict(radius=2, count=128)),
            ),
            partitions=dict(
                layer=dict(
                    type="rhomboid",
                    dimensions=[
                        self.chunk_size[0] * 2,
                        self.chunk_size[1] * 2,
                        self.chunk_size[2] * 2,
                    ],
                )
            ),
            placement=dict(
                random_placement=dict(
                    strategy="bsb.placement.RandomPlacement",
                    partitions=["layer"],
                    cell_types=["pre_cell", "pre_cell2", "test_cell", "single_cell"],
                ),
            ),
            connectivity=dict(),
        )
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)

    def test_negative_presyn_ratio(self):
        with self.assertRaises(ConfigurationError):
            # Test negative presyn ratios
            self.cfg.connectivity.add(
                "glom_ubc",
                dict(
                    strategy="cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
                    presynaptic=dict(cell_types=["pre_cell", "pre_cell2"]),
                    postsynaptic=dict(cell_types=["test_cell"]),
                    radius=self.radius,
                    ratios_ubc=dict(pre_cell=0.2, pre_cell2=-0.4),
                ),
            )

    def test_all_null_presyn_ratios(self):
        with self.assertRaises(ConfigurationError):
            # Test all 0 presyn ratios
            self.cfg.connectivity.add(
                "glom_ubc",
                dict(
                    strategy="cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
                    presynaptic=dict(cell_types=["pre_cell", "pre_cell2"]),
                    postsynaptic=dict(cell_types=["test_cell"]),
                    radius=self.radius,
                    ratios_ubc=dict(
                        pre_cell=0.0,
                        pre_cell2=0.0,
                        pre_cell3=-0.4,
                    ),
                ),
            )

    def test_presyn_ratios(self):
        self.cfg.connectivity.add(
            "glom_ubc",
            dict(
                strategy="cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
                presynaptic=dict(cell_types=["pre_cell", "pre_cell2"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                radius=self.radius,
                ratios_ubc=dict(
                    pre_cell=0.1,
                    pre_cell2=0.3,
                    pre_cell3=-0.4,
                ),
            ),
        )
        predicted = {"pre_cell": 0.25, "pre_cell2": 0.75}
        self.assertEqual(self.cfg.connectivity["glom_ubc"].ratios_ubc.keys(), predicted.keys())
        for k, v in predicted.items():
            self.assertClose(self.cfg.connectivity["glom_ubc"].ratios_ubc[k], v)
        self.network.configuration = self.cfg
        self.network.compile(skip_placement=True, append=True)
        cs1 = self.network.get_connectivity_set("glom_ubc_pre_cell_to_test_cell")
        cs2 = self.network.get_connectivity_set("glom_ubc_pre_cell2_to_test_cell")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()
        cell_targets = np.array([-1, -1])
        self.assertEqual(len(cs2) / len(cs1), 3)
        self.assertEqual(len(cs1) + len(cs2), len(cell_positions))
        for cs in [cs1, cs2]:
            pre_cell_position = self.network.get_placement_set(cs.pre_type.name).load_positions()
            dict_pres = {}
            for from_, to_ in cs.load_connections().as_globals():
                self.assertAll(from_[1:] == cell_targets)
                self.assertAll(to_[1:] == cell_targets)
                post_chunk = np.floor(cell_positions[to_[0]] / self.chunk_size)
                self.assertAll(
                    np.floor(pre_cell_position[from_[0]] / self.chunk_size) - post_chunk <= 1.0,
                    "Chunk size distance should be less than radius",
                )
                post_chunk = str(post_chunk)
                if post_chunk not in dict_pres:
                    dict_pres[post_chunk] = np.zeros(len(pre_cell_position), dtype=int)
                dict_pres[post_chunk][from_[0]] += 1
            for post_chunk, v in dict_pres.items():
                # 64 postsyn cells -> 8 postsyn cell per chunk
                # 80 presyn cells (for each type) -> 10 presyn cell per chunk
                # -> 4 chunks close enough -> 40 presyn cell per postsyn chunk
                # There are more than enough presyn cells so each should be used once per postsyn chunk.
                self.assertAll(
                    np.unique(v) == np.array([0, 1]),
                    f"Each presyn cell of postsyn cells in chunk {post_chunk} should be used at most once",
                )

    def test_uniqueness(self):
        self.cfg.connectivity.add(
            "glom_ubc",
            dict(
                strategy="cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
                presynaptic=dict(cell_types=["single_cell"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                radius=self.radius,
                ratios_ubc=dict(single_cell=1.0),
            ),
        )
        self.network.configuration = self.cfg
        self.network.compile(skip_placement=True, append=True)
        cs = self.network.get_connectivity_set("glom_ubc")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()
        dict_pres = {}
        for from_, to_ in cs.load_connections().as_globals():
            post_chunk = np.floor(cell_positions[to_[0]] / self.chunk_size)
            post_chunk = str(post_chunk)
            if post_chunk not in dict_pres:
                dict_pres[post_chunk] = np.zeros(8, dtype=int)
            dict_pres[post_chunk][from_[0]] += 1
        for post_chunk, v in dict_pres.items():
            # 128 postsyn cells -> 16 postsyn cell per chunk
            # 8 presyn cells -> 1 presyn cell per chunk -> 8 chunks close enough -> 8 presyn cell per postsyn chunk
            # Presyn cells are taken twice each per postsyn chunk
            self.assertAll(
                np.unique(v) == np.array([2]),
                f"Each presyn cell of postsyn cells in chunk {post_chunk} should be used exactly twice",
            )
