import os
import unittest
from test.test_reports import DictTestCase

import numpy as np
from bsb import Configuration, Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture

from cerebellum.analysis.report import LIST_CT_INFO
from cerebellum.analysis.structure_analysis import (
    CellPlacement3D,
    ConnectivityTable,
    PlacementTable,
    StructureReport,
)


class TestPlacementTable(
    RandomStorageFixture, unittest.TestCase, NumpyTestCase, DictTestCase, engine_name="hdf5"
):
    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file("configurations/mouse/mouse_cerebellar_cortex.yaml")
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.plot = PlacementTable(
            (5, 2.5), scaffold=self.scaffold, dict_abv={"granule_cell": "GrC"}
        )

    def test_placement_table(self):
        keys = [
            "glomerulus",
            "mossy_fibers",
            "GrC",
            "golgi_cell",
            "purkinje_cell",
            "basket_cell",
            "stellate_cell",
        ]
        granule_volume = 130.0 * 200.0 * 300.0
        values = [
            granule_volume,
            granule_volume,
            granule_volume,
            granule_volume,
            (15.0 * 200.0 * 300.0),
            (50.0 * 200.0 * 300.0),
            (100.0 * 200.0 * 300.0),
        ]
        self.plot.plot()
        self.assertAll(np.array(self.plot.rows) == np.array(keys))
        expected = np.full((7, 2), "0.00E+00")
        self.assertAll(np.asarray(self.plot.table_values) == expected)
        self.scaffold.compile(only=["granular_layer_innervation"], skip_after_connectivity=True)
        # Placement counts need to be extracted because of bsb randomness
        counts = len(self.scaffold.get_placement_set("mossy_fibers"))
        expected[1][0] = "{:.2E}".format(counts)
        expected[1][1] = "{:.2E}".format(counts / granule_volume)
        self.plot.is_updated = False
        self.plot.plot()
        self.assertAll(np.asarray(self.plot.table_values) == expected)
        self.assertDict(self.plot.get_volumes(), {key: value for key, value in zip(keys, values)})
        dict_counts = {key: 0.0 for key in keys}
        dict_counts["mossy_fibers"] = counts
        self.assertDict(self.plot.get_counts(), dict_counts)
        dict_densities = {
            key: value / volume for (key, value), volume in zip(dict_counts.items(), values)
        }
        self.assertDict(self.plot.get_densities(), dict_densities)

    def test_warn(self):
        scaffold = Scaffold(Configuration.default(), self.storage)
        self.plot.set_scaffold(scaffold)
        with self.assertWarns(UserWarning):
            self.plot.plot()
        self.assertEqual(self.plot.get_counts(), {})
        self.assertEqual(self.plot.get_volumes(), {})
        self.assertEqual(self.plot.get_densities(), {})


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
        keys = np.array(["mf to glomerulus", "mf to mf"])
        self.assertAll(np.array(self.plot.rows) == keys)
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
        self.assertAll(np.asarray(self.plot.table_values)[:, :3] == expected)
        self.assertDictEqual(
            {k: v for k, v in zip(keys, [counts1, counts2])}, self.plot.get_synapse_counts()
        )
        for k in keys:
            self.assertAll(np.array(self.plot.get_nb_synapse_per_pair()[k]) == 1.0)
            self.assertAll(np.array(self.plot.get_convergences()[k]) == 1.0)
        divergences = self.plot.get_divergences()["mf to glomerulus"]
        self.assertEqual(len(divergences), len(self.scaffold.get_placement_set("mossy_fibers")))
        self.assertAlmostEqual(np.mean(divergences), 20.0, delta=0.5)

    def test_warn(self):
        # No connection sets in storage before compile
        self.plot.set_scaffold(self.scaffold)
        with self.assertWarns(UserWarning):
            self.plot.plot()
        self.assertEqual(self.plot.get_synapse_counts(), {})
        self.assertEqual(self.plot.get_nb_synapse_per_pair(), {})
        self.assertEqual(self.plot.get_convergences(), {})
        self.assertEqual(self.plot.get_divergences(), {})


class TestCellPlacement3D(
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
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold.compile(
            only=[
                "granular_layer_innervation",
                "granular_layer_placement",
            ],
            skip_after_connectivity=True,
        )

    def test_placement(self):
        old_color = np.array([0.847, 0, 0.451, 0.5])
        new_color = np.array([0.847, 0, 0.451, 1.0])
        plot = CellPlacement3D(
            (5, 5),
            scaffold=self.scaffold,
            ignored_ct=["mossy_fibers"],
            dict_colors={"glomerulus": old_color},
        )
        plot.plot()
        glom_ps = self.scaffold.get_placement_set("glomerulus")
        scat_data = plot.get_ax().collections[0]
        self.assertEqual(len(glom_ps), len(scat_data._facecolors))
        self.assertAll(old_color == scat_data._facecolors[0])
        self.assertEqual(2.25, scat_data._sizes)
        self.assertAll(glom_ps.load_positions() == np.array(scat_data._offsets3d).T)
        plot.set_color("glomerulus", new_color[:3])
        plot.plot()
        self.assertAll(new_color == plot.get_ax().collections[0]._facecolors[0])
        plot.dict_colors = {}
        plot.plot()
        self.assertAll(
            np.array([0.6, 0.6, 0.6, 1.0]) == plot.get_ax().collections[0]._facecolors[0]
        )

    def test_empty_cell_placement(self):
        color = np.array([0.847, 0, 0.451, 1.0])
        # all placed cells are ignored
        plot = CellPlacement3D(
            (5, 5),
            scaffold=self.scaffold,
            dict_colors={"glomerulus": color, "mossy_fibers": color},
        )
        self.assertEqual(len(plot.get_ax().collections), 0)


class TestStructureReport(
    RandomStorageFixture, NumpyTestCase, unittest.TestCase, engine_name="hdf5"
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
            "postsynaptic": {"cell_types": ["glomerulus"]},
            "x_length": 60,
            "y_length": 20,
        }
        self.auto_filename = "test_auto_report.pdf"
        self.cfg.after_connectivity["print_structure_report"].output_filename = self.auto_filename
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold.compile(
            only=[
                "granular_layer_innervation",
                "granular_layer_placement",
                "mossy_fibers_to_glomerulus",
            ],
            # to check that the automatic struct report is saved
            skip_after_connectivity=False,
        )
        self.report = StructureReport(self.scaffold, LIST_CT_INFO)

    def test_structure_report(self):
        plot_keys = np.array(["density_table", "connectivity_table", "placement_3d", "legend"])
        self.assertAll(np.array(list(self.report.plots.keys())) == plot_keys)
        filename = "test_report.pdf"
        self.report.print_report(filename, dpi=100)
        # should be seven cell types
        self.assertEqual(len(self.report.plots["density_table"].table_values), 7)
        # only one connectivity rule
        self.assertEqual(len(self.report.plots["connectivity_table"].table_values), 1)
        # no cell displayed because all ignored
        self.assertEqual(len(self.report.plots["placement_3d"].get_ax().collections), 0)
        # only 6 cell types in the legend
        self.assertEqual(len(self.report.plots["legend"].get_ax().legend_.legend_handles), 6)
        self.assertTrue(filename in os.listdir())
        os.remove(filename)

        # Test automatic report
        self.assertTrue(self.auto_filename in os.listdir())
        os.remove(self.auto_filename)
