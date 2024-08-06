import os
import unittest

import numpy as np
from bsb import Configuration, Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture

from cerebellum.analysis.plots import Legend, ScaffoldPlot
from cerebellum.analysis.report import LIST_CT_INFO, BSBReport, Report


class TestReport(unittest.TestCase, NumpyTestCase):

    def _test_dict(self, tested, expected):
        self.assertAll(np.array(list(tested.keys())) == np.array(list(expected.keys())))
        self.assertAll(np.array(list(tested.values())) == np.array(list(expected.values())))

    def test_report(self):
        report = Report()
        self.assertEqual({}, report.abbreviations)
        self.assertEqual({}, report.colors)
        with self.assertRaises(ValueError):
            report.set_color("blue", [0.0, 0.0])

        old_blue = [0.0, 0.0, 1.0]
        new_blue = [0.1, 0.1, 0.9]
        plot = Legend((7.5, 6.8), 2, None)
        plot2 = Legend((7.5, 6.8), 2, None)

        report.set_color("blue", old_blue)
        self._test_dict(report.colors, {"blue": old_blue})
        self._test_dict(report.abbreviations, {"blue": "blue"})

        report.add_plot("legend", plot)
        self._test_dict(plot.dict_colors, {"blue": old_blue})
        self._test_dict(report.plots, {"legend": plot})
        with self.assertWarns(UserWarning):
            report.add_plot("legend", plot2)
        self._test_dict(report.plots, {"legend": plot})

        report.set_color("blue", new_blue)
        self._test_dict(report.colors, {"blue": new_blue})
        self._test_dict(plot.dict_colors, {"blue": new_blue})

        report.set_plot_colors(plot2)
        self._test_dict(plot2.dict_colors, {"blue": new_blue})

    def test_print_report(self):
        plot = Legend((7.5, 6.8), 2, None)
        filename = "test_figure.png"
        self.assertFalse(plot.is_plotted)
        self.assertFalse(plot.is_updated)
        report = Report(LIST_CT_INFO)
        report.add_plot("legend", plot)
        self._test_dict({i.name: i.color for i in LIST_CT_INFO}, report.colors)
        report.save_plot("bla", filename, 100)
        self.assertFalse(filename in os.listdir())
        report.save_plot("legend", filename, 100)
        self.assertTrue(plot.is_plotted)
        self.assertTrue(plot.is_updated)
        self.assertTrue(filename in os.listdir())
        os.remove(filename)
        filename = "test_report.pdf"
        plot.is_plotted = False
        report.print_report(filename, dpi=100)
        self.assertTrue(filename in os.listdir())
        os.remove(filename)


class TestBSBReport(RandomStorageFixture, unittest.TestCase, NumpyTestCase, engine_name="hdf5"):
    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file("configurations/mouse/mouse_cerebellar_cortex.yaml")
        self.cfg2 = Configuration.default()
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.scaffold2 = Scaffold(self.cfg2, self.storage)
        self.report = BSBReport(self.scaffold, LIST_CT_INFO)

    def test_bsb_report(self):
        plot = ScaffoldPlot((10, 10), self.scaffold)
        plot.is_updated = True
        plot2 = Legend((7.5, 6.8), 2, None)
        plot2.is_updated = True

        plot3 = ScaffoldPlot((10, 10), self.scaffold2)
        plot3.is_updated = True
        self.report.add_plot("scaffold1", plot)
        self.assertTrue(plot.is_updated)
        self.assertEqual(self.scaffold, plot.scaffold)
        self.report.add_plot("legend", plot2)
        self.assertTrue(plot2.is_updated)
        self.report.add_plot("scaffold2", plot3)
        self.assertFalse(plot3.is_updated)
        self.assertEqual(self.scaffold, plot3.scaffold)

        self.assertAll(
            np.array(self.report.cell_names)
            == np.array(
                [
                    "glomerulus",
                    "mossy_fibers",
                    "granule_cell",
                    "golgi_cell",
                    "purkinje_cell",
                    "basket_cell",
                    "stellate_cell",
                ]
            )
        )

        # test that we can create a report from filename
        BSBReport(self.scaffold.storage.root, LIST_CT_INFO)
