import os
import unittest
from unittest.mock import patch

import numpy as np
from bsb import Scaffold, parse_configuration_file
from bsb_test import NumpyTestCase, RandomStorageFixture
from matplotlib import pyplot as plt

from cerebellum.analysis.plots import Legend, Plot, ScaffoldPlot


def mock_update(self):
    self.update_counter += 1
    self.is_updated = True


class TestPlot(unittest.TestCase, NumpyTestCase):

    def setUp(self):
        self.colors = {"blue": np.array([0, 0, 1]), "red": np.array([1, 0, 0])}
        self.plot = Plot((7.5, 6.8), 2, 3, self.colors)

    def tearDown(self):
        plt.close("all")

    def test_init(self):
        self.plot.init_plot()
        self.assertAll(
            np.array(self.plot.figure.bbox_inches.bounds) == np.array([0.0, 0.0, 7.5, 6.8])
        )
        plot = Plot((2, 2), 1, 6)
        plot.init_plot()
        plot2 = Plot((2, 2))
        plot2.init_plot()
        with self.assertRaises(ValueError):
            Plot((2, 2), nb_rows=1, nb_cols=0)
        self.assertEqual(self.plot.axes.shape, (2, 3))
        self.assertEqual(len(plot.axes), 6)
        self.assertAll(np.array(self.plot.figure.axes) == np.array(self.plot.get_axes()))
        self.assertAll(np.array(plot.get_axes()) == np.array(plot.figure.axes))
        self.assertAll(np.array(plot2.get_axes()) == np.array(plot2.figure.axes))
        for i in range(6):
            self.assertEqual(self.plot.get_ax(i), self.plot.figure.axes[i])
            self.assertEqual(plot.get_ax(i), plot.figure.axes[i])
        self.assertEqual(plot2.get_ax(), plot2.figure.axes[0])
        with self.assertRaises(IndexError):
            plot.get_ax(-1)
        with self.assertRaises(IndexError):
            plot.get_ax(6)
        with self.assertRaises(ValueError):
            Plot((1, 1), dict_colors={"black": [0.0, 0.0]})

    def test_set_color(self):
        with self.assertRaises(ValueError):
            self.plot.set_color(
                "black",
                [
                    0.0,
                ],
            )

        self.assertAll(np.array(list(self.plot.dict_colors.keys())) == np.array(list(self.colors)))
        new_red = [0.9, 0.1, 0.1]
        new_black = np.array([0.0, 0.0, 0.0, 1.0])
        self.plot.set_color("red", new_red)
        self.assertAll(self.plot.dict_colors["red"] == np.array(new_red))
        self.plot.set_color("black", new_black)
        self.assertAll(
            np.array(list(self.plot.dict_colors.keys())) == np.array(["blue", "red", "black"])
        )

    @patch.object(Plot, "update", mock_update)
    def test_plot(self):
        self.plot.update_counter = 0
        self.assertEqual(self.plot.is_plotted, False)
        self.assertEqual(self.plot.is_updated, False)
        self.plot.plot()
        self.assertEqual(self.plot.is_plotted, True)
        self.assertEqual(self.plot.is_updated, True)
        self.assertEqual(self.plot.update_counter, 1)

        # Create a line to test if the plots are cleared
        self.plot.get_ax().plot([0, 1], [0, 1])
        self.assertEqual(len(self.plot.get_ax().lines), 1)
        for ax in self.plot.get_axes():
            self.assertTrue(ax.axison)
        self.plot.set_axis_off()
        for ax in self.plot.get_axes():
            self.assertFalse(ax.axison)
        self.plot.plot()
        self.assertEqual(self.plot.update_counter, 1, "second plot() should not call update")
        self.assertEqual(
            len(self.plot.get_ax().lines), 0, "second plot() should have cleared the axes"
        )

        # Create a line to test if the plots are cleared after new color
        self.plot.get_ax().plot([0, 1], [0, 1])
        self.plot.set_color("black", [0.0, 0.0, 0.0])
        self.plot.plot()
        self.assertEqual(self.plot.update_counter, 1, "third plot() should not call update")
        self.assertEqual(
            len(self.plot.get_ax().lines), 0, "second plot() should have cleared the axes"
        )

        filename = "test_figure.png"
        self.plot.is_plotted = False
        self.plot.is_updated = False
        self.plot.save_figure(filename, 100)
        self.assertEqual(self.plot.update_counter, 2, "save fig should call update")
        self.assertTrue(filename in os.listdir())
        os.remove(filename)


class TestLegend(unittest.TestCase, NumpyTestCase):
    def setUp(self):
        self.colors = {
            "blue": np.array([0, 0, 1]),
            "red": np.array([1, 0, 0]),
            "green": np.array([0, 1, 0]),
            "white": np.array([1, 1, 1]),
        }
        self.plot = Legend((7.5, 6.8), 2, self.colors)

    def test_remove_ct(self):
        to_keep = ["blue", "red_cell", "black", "green"]
        to_remove = ["white", "purple", "green"]
        self.plot.remove_ct(to_keep, to_remove)
        self.assertAll(np.array(list(self.plot.dict_colors.keys())) == np.array(["blue", "red"]))

    def test_plot(self):
        self.assertEqual(self.plot.nb_cols, 1)
        self.assertEqual(self.plot.nb_rows, 1)
        self.assertTrue(self.plot.get_ax().legend_ is None)
        self.plot.plot()
        self.assertFalse(self.plot.get_ax().legend_ is None)


class TestScaffoldPlot(RandomStorageFixture, unittest.TestCase, NumpyTestCase, engine_name="hdf5"):
    def setUp(self):
        super().setUp()
        self.cfg = parse_configuration_file("configurations/mouse/mouse_cerebellar_cortex.yaml")
        self.scaffold = Scaffold(self.cfg, self.storage)
        self.plot = ScaffoldPlot((10, 10), self.scaffold)

    def test_set_scaffold(self):
        self.assertEqual(self.plot.scaffold, self.scaffold)
        self.plot.plot()
        self.assertFalse(self.plot.set_scaffold(self.scaffold))
        self.assertEqual(self.plot.scaffold, self.scaffold)
        self.assertTrue(self.plot.is_plotted)
        self.assertTrue(self.plot.is_updated)
        scaffold2 = Scaffold(self.cfg, self.storage)
        self.assertTrue(self.plot.set_scaffold(scaffold2))
        self.assertFalse(self.plot.is_plotted)
        self.assertFalse(self.plot.is_updated)
