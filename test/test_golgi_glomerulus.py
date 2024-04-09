"""
Unit tests of the `ConnectomeGlomerulusGranule` strategy
"""

import unittest
from os.path import abspath, dirname, join

import numpy as np
from bsb import Configuration, ConfigurationError, Scaffold, WorkflowError
from bsb_test import NetworkFixture, NumpyTestCase, RandomStorageFixture

from cerebellum.connectome.golgi_glomerulus import ConnectomeGolgiGlomerulus


class TestGlomerulusGranule(
    RandomStorageFixture,
    NetworkFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.radius = 40
        self.divergence = 4
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
                    file=abspath(join(dirname(dirname(__file__)), "morphologies", "GolgiCell.swc")),
                    parser=dict(tags={16: ["dendrites", "basal_dendrites"]}),
                ),
            ],
            cell_types=dict(
                glom_cell=dict(spatial=dict(radius=2, count=100)),
                error_cell=dict(spatial=dict(radius=2, count=1)),
                golgi_cell=dict(spatial=dict(radius=2, count=5, morphologies=["GolgiCell"])),
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
                    cell_types=["glom_cell", "error_cell", "golgi_cell"],
                ),
            ),
            connectivity=dict(
                glom_to_post=dict(
                    strategy="cerebellum.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi",
                    presynaptic=dict(cell_types=["glom_cell"]),
                    postsynaptic=dict(
                        cell_types=["golgi_cell"], morphology_labels=["basal_dendrites"]
                    ),
                    radius=self.radius,
                ),
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)

    def test_error_golgi_glom_strat(self):
        with self.assertRaises(ConfigurationError):
            # Test postsyn does not match dependency postsynaptic cell
            self.cfg.connectivity.add(
                "golgi_glom",
                dict(
                    strategy="cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
                    presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["axon"]),
                    postsynaptic=dict(cell_types=["error_cell"]),
                    radius=self.radius,
                    divergence=self.divergence,
                    glom_cell_type="glom_cell",
                    glom_post_strat="glom_to_post",
                ),
            )
        with self.assertRaises(ConfigurationError):
            # Test glom_cell_type is a presynaptic cell type of dependency rule
            self.cfg.connectivity.add(
                "golgi_glom2",
                dict(
                    strategy="cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
                    presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["axon"]),
                    postsynaptic=dict(cell_types=["golgi_cell"]),
                    radius=self.radius,
                    divergence=self.divergence,
                    glom_cell_type="error_cell",
                    glom_post_strat="glom_to_post",
                ),
            )
        # Test wrong morphology label
        self.cfg.connectivity.add(
            "golgi_glom3",
            dict(
                strategy="cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
                presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["wrong_label"]),
                postsynaptic=dict(cell_types=["golgi_cell"]),
                radius=self.radius,
                divergence=self.divergence,
                glom_cell_type="glom_cell",
                glom_post_strat="glom_to_post",
            ),
        )
        # Fixme: test the sub-error: should be ConnectivityError
        with self.assertRaises(WorkflowError):
            self.network.compile(append=True, skip_placement=True)

    def test_golgi_glom_strat(self):
        self.network.connectivity["golgi_glom"] = ConnectomeGolgiGlomerulus(
            presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["axon"]),
            postsynaptic=dict(cell_types=["golgi_cell"]),
            radius=self.radius,
            divergence=self.divergence,
            glom_cell_type="glom_cell",
            depends_on=[self.cfg.connectivity["glom_to_post"]],
            glom_post_strat="glom_to_post",
        )
        self.network.compile(append=True, skip_placement=True)

    def test_extra_depends_on(self):
        self.cfg.connectivity.add(
            "golgi_glom",
            dict(
                strategy="cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
                presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["axon"]),
                postsynaptic=dict(cell_types=["golgi_cell"]),
                radius=self.radius,
                divergence=self.divergence,
                glom_cell_type="glom_cell",
                depends_on=[self.cfg.connectivity["glom_to_post"]],
                glom_post_strat="glom_to_post",
            ),
        )
        self.assertEqual(len(self.cfg.connectivity["golgi_glom"].depends_on), 1)
