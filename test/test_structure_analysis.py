import unittest

import numpy as np
from bsb import Configuration, RandomPlacement, Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture

from cerebellum.analysis.structure_analysis import ConnectivityTable, PlacementTable


class TestPlacementTable(
    RandomStorageFixture, unittest.TestCase, NumpyTestCase, engine_name="hdf5"
):
    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file("configurations/mouse/mouse_cerebellar_cortex.yaml")
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.plot = PlacementTable(
            (5, 2.5), scaffold=self.scaffold, dict_abv={"granule_cell": "GrC"}
        )

    def test_placement_table(self):
        self.plot.plot()
        self.assertAll(
            np.array(self.plot.rows)
            == np.array(
                [
                    "glomerulus",
                    "mossy_fibers",
                    "GrC",
                    "golgi_cell",
                    "purkinje_cell",
                    "basket_cell",
                    "stellate_cell",
                ]
            )
        )
        expected = np.full((7, 2), "0.00E+00")
        self.assertAll(np.asarray(self.plot.values) == expected)
        self.scaffold.compile(only=["granular_layer_innervation"])
        # Placement counts need to be extracted because of bsb randomness
        counts = len(self.scaffold.get_placement_set("mossy_fibers"))
        expected[1][0] = "{:.2E}".format(counts)
        expected[1][1] = "{:.2E}".format(counts / (130 * 200 * 300))
        self.plot.is_updated = False
        self.plot.plot()
        self.assertAll(np.asarray(self.plot.values) == expected)

    def test_warn(self):
        scaffold = Scaffold(Configuration.default(), self.storage)
        self.plot.set_scaffold(scaffold)
        with self.assertWarns(UserWarning):
            self.plot.plot()


class TestConnectivityTable(
    RandomStorageFixture, unittest.TestCase, NumpyTestCase, engine_name="hdf5"
):
    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file("configurations/mouse/mouse_cerebellar_cortex.yaml")
        self.cfg.placement["granular_layer_placement"] = {
            "strategy": "bsb.placement.RandomPlacement",
            "partitions": ["granular_layer"],
            "cell_types": ["glomerulus"],
        }
        self.cfg.connectivity.clear()
        self.cfg.connectivity["mossy_fibers_to_glomerulus"] = {
            "strategy": "cerebellum.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
            "presynaptic": {"cell_types": ["mossy_fibers"]},
            # to test complex connectivity set names
            "postsynaptic": {"cell_types": ["glomerulus", "mossy_fibers"]},
            "x_length": 60,
            "y_length": 20,
        }
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.plot = ConnectivityTable(
            (5, 2.5), scaffold=self.scaffold, dict_abv={"mossy_fibers": "mf"}
        )

    def test_connectivity_table(self):
        self.scaffold.compile(
            only=[
                "granular_layer_innervation",
                "granular_layer_placement",
                "mossy_fibers_to_glomerulus",
            ],
            skip_after_connectivity=True,
        )
        self.plot.plot()
        self.assertAll(np.array(self.plot.rows) == np.array(["mf to glomerulus", "mf to mf"]))
        # Connectivity counts need to be extracted because of bsb randomness
        counts1 = len(
            self.scaffold.get_connectivity_set(
                "mossy_fibers_to_glomerulus_mossy_fibers_to_glomerulus"
            )
        )
        counts2 = len(
            self.scaffold.get_connectivity_set(
                "mossy_fibers_to_glomerulus_mossy_fibers_to_mossy_fibers"
            )
        )
        expected = np.array(
            [
                [
                    str(counts1),
                    "{:.2} $\pm$ {:.2}".format(1.0, 0.0),
                    "{:.2} $\pm$ {:.2}".format(1.0, 0.0),
                ],
                [
                    str(counts2),
                    "{:.2} $\pm$ {:.2}".format(1.0, 0.0),
                    "{:.2} $\pm$ {:.2}".format(1.0, 0.0),
                ],
            ]
        )
        self.assertAll(np.asarray(self.plot.values)[:, :3] == expected)

    def test_warn(self):
        # No connection sets in storage before compile
        self.plot.set_scaffold(self.scaffold)
        with self.assertWarns(UserWarning):
            self.plot.plot()
