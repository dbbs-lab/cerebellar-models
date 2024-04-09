"""
Unit tests of the `ConnectomeGlomerulusGolgi` strategy
"""

import unittest

import numpy as np
import yaml
from bsb import Configuration, Scaffold, WorkflowError
from bsb_test import NetworkFixture, NumpyTestCase, RandomStorageFixture

from cerebellum.connectome.glomerulus_golgi import ConnectomeGlomerulusGolgi


class TestGlomerulusGolgi(
    RandomStorageFixture,
    NetworkFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.radius = 40
        self.chunk_size = np.array([30, 30, 30])
        self.cfg = Configuration.default(
            network=dict(
                chunk_size=self.chunk_size,
                x=self.chunk_size[0] * 2,
                y=self.chunk_size[1] * 2,
                z=self.chunk_size[2] * 2,
            ),
            morphologies=[
                dict(
                    file="../morphologies/GolgiCell.swc",
                    parser=dict(tags={16: ["dendrites", "basal_dendrites"]}),
                )
            ],
            cell_types=dict(
                pre_cell=dict(spatial=dict(radius=2, count=100)),
                test_cell=dict(spatial=dict(radius=2, count=5, morphologies=["GolgiCell"])),
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
                    cell_types=["pre_cell", "test_cell"],
                ),
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)

    def test_error_mf_golgi_strat(self):
        self.network.connectivity["glom_to_golgi"] = ConnectomeGlomerulusGolgi(
            presynaptic=dict(cell_types=["pre_cell"]),
            postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["wrong_label"]),
            radius=self.radius,
        )
        # Fixme: test the sub-error: should be ConnectivityError
        with self.assertRaises(WorkflowError):
            self.network.compile(append=True, skip_placement=True)

    def test_mf_golgi_strat(self):
        self.network.connectivity["glom_to_golgi"] = ConnectomeGlomerulusGolgi(
            presynaptic=dict(cell_types=["pre_cell"]),
            postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["basal_dendrites"]),
            radius=self.radius,
        )
        self.network.compile(append=True, skip_placement=True)
        cs = self.network.get_connectivity_set("glom_to_golgi")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()
        pre_cell_positions = self.network.get_placement_set("pre_cell").load_positions()
        morpho = self.network.cell_types["test_cell"].morphologies[0].load()

        for from_, to_ in cs.load_connections().as_globals():
            self.assertTrue(len(morpho.branches) > to_[1])
            branch = morpho.branches[to_[1]]
            self.assertTrue(len(branch) == to_[2] + 1, "Target should be last branch point")
            self.assertTrue(
                "basal_dendrites" in branch.labelsets[branch.labels[to_[2]]],
                "Target branch should be a dendrite",
            )
            self.assertAll(
                from_[1:] == np.array([-1, -1]), "Presyn has no morpho, locs should be -1."
            )
            self.assertTrue(
                np.linalg.norm(pre_cell_positions[from_[0]] - cell_positions[to_[0]]) <= self.radius
            )
            self.assertTrue(
                np.linalg.norm(
                    (
                        np.floor(pre_cell_positions[from_[0]] / self.chunk_size)
                        - np.floor(cell_positions[to_[0]] / self.chunk_size)
                    )
                    * self.chunk_size
                )
                <= self.radius,
                "Chunk size distance should be less than radius",
            )
