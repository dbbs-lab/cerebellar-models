"""
Unit tests of the `ConnectomeGlomerulusGranule` strategy
"""

import unittest

import numpy as np
from bsb import (
    Configuration,
    ConfigurationError,
    RequirementError,
    Scaffold,
    WorkflowError,
)
from bsb_test import NetworkFixture, NumpyTestCase, RandomStorageFixture

from cerebellum.connectome.glomerulus_granule import ConnectomeGlomerulusGranule


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
        self.convergence = 4
        self.cfg = Configuration.default(
            morphologies=[dict(file="../morphologies/GranuleCell.swc")],
            cell_types=dict(
                pre_cell=dict(spatial=dict(radius=2, count=100)),
                pre_cell2=dict(spatial=dict(radius=2, count=1)),
                test_cell=dict(spatial=dict(radius=2, count=100, morphologies=["GranuleCell"])),
            ),
            partitions=dict(layer=dict(type="rhomboid", dimensions=[100.0, 100.0, 100.0])),
            placement=dict(
                random_placement=dict(
                    strategy="bsb.placement.RandomPlacement",
                    partitions=["layer"],
                    cell_types=["pre_cell", "pre_cell2", "test_cell"],
                ),
            ),
            connectivity=dict(
                x_to_glomerulus=dict(
                    strategy="cerebellum.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
                    presynaptic=dict(cell_types=["pre_cell"]),
                    postsynaptic=dict(cell_types=["test_cell"]),
                    x_length=60,
                    z_length=20,
                ),
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)

    def test_error_mf_glom_strat(self):
        with self.assertRaises(RequirementError):
            self.cfg.connectivity.add(
                "glom_to_gran",
                dict(
                    strategy="cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
                    presynaptic=dict(cell_types=["test_cell"]),
                    postsynaptic=dict(cell_types=["test_cell"]),
                    radius=self.radius,
                    convergence=self.convergence,
                    mf_cell_type="pre_cell",
                ),
            )
        with self.assertRaises(ConfigurationError):
            self.cfg.connectivity.add(
                "glom_to_gran",
                dict(
                    strategy="cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
                    presynaptic=dict(cell_types=["pre_cell"]),
                    postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["dendrites"]),
                    radius=self.radius,
                    convergence=self.convergence,
                    mf_cell_type="pre_cell",
                    mf_glom_strat="x_to_glomerulus",
                ),
            )

        self.cfg.connectivity["x_to_glomerulus"].postsynaptic.cell_types = [
            self.cfg.cell_types["pre_cell"],
            self.cfg.cell_types["pre_cell2"],
        ]
        with self.assertRaises(ConfigurationError):
            self.cfg.connectivity.add(
                "glom_to_gran2",
                dict(
                    strategy="cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
                    presynaptic=dict(cell_types=["test_cell"]),
                    postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["dendrites"]),
                    radius=self.radius,
                    convergence=self.convergence,
                    mf_cell_type="pre_cell",
                    mf_glom_strat="x_to_glomerulus",
                ),
            )

        del (
            self.cfg.connectivity["x_to_glomerulus"],
            self.cfg.connectivity["glom_to_gran"],
            self.cfg.connectivity["glom_to_gran2"],
        )
        self.cfg.connectivity.add(
            "x_to_glomerulus2",
            dict(
                strategy="cerebellum.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
                presynaptic=dict(cell_types=["pre_cell2"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                x_length=60,
                z_length=20,
            ),
        )
        self.cfg.connectivity.add(
            "glom_to_gran_f",
            dict(
                strategy="cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
                presynaptic=dict(cell_types=["test_cell"]),
                postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["dendrites"]),
                radius=self.radius,
                convergence=self.convergence,
                mf_cell_type="pre_cell2",
                mf_glom_strat="x_to_glomerulus2",
            ),
        )
        network = Scaffold(self.cfg, self.storage)
        # Fixme: test the sub-error: should be TooFewGlomeruliClusters
        with self.assertRaises(WorkflowError):
            network.compile(clear=True)

    def test_mf_glom_strat(self):
        self.network.connectivity["glom_to_gran"] = ConnectomeGlomerulusGranule(
            presynaptic=dict(cell_types=["test_cell"]),
            postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["dendrites"]),
            radius=self.radius,
            convergence=self.convergence,
            mf_cell_type="pre_cell",
            mf_glom_strat="x_to_glomerulus",
        )
        self.network.compile(append=True, skip_placement=True)
        cs = self.network.get_connectivity_set("glom_to_gran")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()
        self.assertEqual(
            len(cs),
            len(cell_positions) * self.convergence,
            "As many connection as presyn cell times convergence",
        )
        morpho = self.network.cell_types["test_cell"].morphologies[0].load()
        for from_, to_ in cs.load_connections().as_globals():
            self.assertTrue(
                np.linalg.norm(cell_positions[from_[0]] - cell_positions[to_[0]]) <= self.radius
            )
            self.assertTrue(len(morpho.branches) > to_[1])
            branch = morpho.branches[to_[1]]
            self.assertTrue(len(branch) == to_[2] + 1)
            self.assertTrue("dendrites" in branch.labelsets[branch.labels[to_[2]]])
            self.assertAll(from_[1:] == np.array([-1, -1]))

    def test_extra_depends_on(self):
        self.cfg.connectivity.add(
            "glom_to_gran",
            dict(
                strategy="cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
                presynaptic=dict(cell_types=["test_cell"]),
                postsynaptic=dict(cell_types=["test_cell"], morphology_labels=["dendrites"]),
                radius=self.radius,
                convergence=self.convergence,
                mf_cell_type="pre_cell",
                depends_on=[self.cfg.connectivity["x_to_glomerulus"]],
                mf_glom_strat="x_to_glomerulus",
            ),
        )
        self.assertEqual(len(self.cfg.connectivity["glom_to_gran"].depends_on), 1)
