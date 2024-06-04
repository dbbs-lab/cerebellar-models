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
        # radius greater than one chunk
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
                # one glom and golgi per chunk
                glom_cell=dict(spatial=dict(radius=2, density=1.0 / np.prod(self.chunk_size))),
                error_cell=dict(spatial=dict(radius=2, count=1)),
                golgi_cell=dict(
                    spatial=dict(
                        radius=2, density=1.0 / np.prod(self.chunk_size), morphologies=["GolgiCell"]
                    )
                ),
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
                    # radius is less than a chunk so that glom connect to only onw golgi
                    radius=np.min(self.chunk_size) - 1.0,
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
                    glom_cell_types=["glom_cell"],
                    glom_post_strats=["glom_to_post"],
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
                    glom_cell_types=["error_cell"],
                    glom_post_strats=["glom_to_post"],
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
                glom_cell_types=["glom_cell"],
                glom_post_strats=["glom_to_post"],
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
            glom_cell_types=["glom_cell"],
            glom_post_strats=["glom_to_post"],
        )
        self.network.compile(append=True, skip_placement=True)

        cs = self.network.get_connectivity_set("golgi_glom")
        glom_locs, post_locs = (
            self.network.get_connectivity_set("glom_to_post").load_connections().as_globals().all()
        )
        cell_positions = self.network.get_placement_set("golgi_cell").load_positions()
        glom_positions = self.network.get_placement_set("glom_cell").load_positions()
        # each glom should connect to one postsyn so max nb connection per golgi is divergence
        self.assertTrue(
            len(cs) <= len(cell_positions) * self.divergence,
            "Maximum nb connection pre golgi cell should be divergence",
        )
        morpho = self.network.cell_types["golgi_cell"].morphologies[0].load()
        post_gloms = np.full((len(cell_positions), self.divergence), -1)
        last_glom = np.zeros(len(cell_positions), dtype=int)
        for from_, to_ in cs.load_connections().as_globals():
            self.assertTrue(len(morpho.branches) > from_[1])
            branch = morpho.branches[from_[1]]
            self.assertTrue(len(branch) == from_[2] + 1, "Target should be last branch point")
            self.assertTrue(
                "axon" in branch.labelsets[branch.labels[from_[2]]],
                "Source branch should be a axon",
            )
            filter_glom = np.where(post_locs[:, 0] == to_[0])[0]
            # There should be only one glom per postsyn cell due to cfg.
            self.assertTrue(len(filter_glom) == 1)
            loc_glom = glom_locs[filter_glom[0], 0]
            self.assertTrue(last_glom[from_[0]] < self.divergence)
            post_gloms[from_[0], last_glom[from_[0]]] = loc_glom
            last_glom[from_[0]] += 1
            self.assertAll(
                to_[1:] == post_locs[filter_glom[0], 1:],
                "Postsyn targets on morpho should be the same as dependency's connections.",
            )
            self.assertTrue(
                np.linalg.norm(cell_positions[from_[0]] - glom_positions[loc_glom]) <= self.radius,
                "Distance between golgi and glom should be less than radius",
            )
        # all intermediate cells should be unique
        self.assertAll(
            np.array([np.unique(x[:i]).size == i for i, x in zip(last_glom, post_gloms)])
        )

    def test_extra_depends_on(self):
        self.cfg.connectivity.add(
            "golgi_glom",
            dict(
                strategy="cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
                presynaptic=dict(cell_types=["golgi_cell"], morphology_labels=["axon"]),
                postsynaptic=dict(cell_types=["golgi_cell"]),
                radius=self.radius,
                divergence=self.divergence,
                glom_cell_types=["glom_cell"],
                depends_on=[self.cfg.connectivity["glom_to_post"]],
                glom_post_strats=["glom_to_post"],
            ),
        )
        self.assertEqual(len(self.cfg.connectivity["golgi_glom"].depends_on), 1)
