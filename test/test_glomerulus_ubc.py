import unittest

import numpy as np
from bsb import Configuration, ConfigurationError, Scaffold, WorkflowError
from bsb_test import NetworkFixture, NumpyTestCase, RandomStorageFixture


class TestGlomerulus_to_UBC(
    RandomStorageFixture, NetworkFixture, NumpyTestCase, unittest.TestCase, engine_name="hdf5"
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
            cell_types=dict(
                pre_cell=dict(spatial=dict(radius=2, count=100)),
                pre_cell2=dict(spatial=dict(radius=2, count=100)),
                single_cell=dict(spatial=dict(radius=2, count=10)),
                test_cell=dict(spatial=dict(radius=2, count=100)),
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
        self.assertClose(len(cs2) / len(cs1), 3, atol=3e-1)
        self.assertClose(len(cs1) + len(cs2), 100, atol=1)
