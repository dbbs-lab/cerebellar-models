import unittest

import numpy as np
from bsb import Configuration, Scaffold
from bsb_test import NumpyTestCase, RandomStorageFixture


class TestMicrozones(
    RandomStorageFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.chunk_size = np.array([30, 30, 30])
        self.labels = ["plus", "minus"]
        self.cfg = Configuration.default(
            network=dict(
                chunk_size=self.chunk_size,
                x=self.chunk_size[0] * 2,
                y=self.chunk_size[1] * 2,
                z=self.chunk_size[2] * 2,
            ),
            cell_types=dict(
                io=dict(spatial=dict(radius=2, count=16)),
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
                    cell_types=["io"],
                ),
            ),
            after_placement=dict(
                label_cell=dict(
                    strategy="cerebellar_models.placement.microzones.LabelMicrozones",
                    cell_types=["io"],
                    labels=self.labels,
                ),
            ),
        )
        self.scaffold = Scaffold(self.cfg, self.storage)

    def _check_labels(self, axis=0):
        ps = self.scaffold.get_placement_set("io")
        unique_labels = [
            list(label)[0] for label in self.scaffold.get_placement_set("io").get_unique_labels()
        ]
        self.assertAll(np.isin(self.labels, unique_labels))
        current_min = np.array([0.0, 0.0, 0.0])
        current_max = np.array(self.scaffold.partitions["layer"].dimensions)
        current_max[axis] /= 2
        for label in unique_labels:
            filt = ps.get_labelled([label])
            self.assertEqual(filt.size, 8)
            positions = ps.load_positions()[filt]
            self.assertAll(positions >= current_min)
            self.assertAll(positions < current_max)
            current_max[axis] += self.chunk_size[axis]
            current_min[axis] += self.chunk_size[axis]

    def test_microzones(self):
        self.scaffold.compile()
        self._check_labels()

    def test_microzones_axis(self):
        self.cfg.after_placement["label_cell"].axis = 2
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold.compile()
        self._check_labels(axis=2)
