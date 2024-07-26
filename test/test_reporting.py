import os
import unittest
from unittest.mock import patch

import numpy as np
from analysis.plots import Legend, Plot
from bsb_test import NumpyTestCase


def mock_update(self):
    self.update_counter += 1
    self.is_updated = True


class TestPlot(unittest.TestCase, NumpyTestCase):

    def setUp(self):
        self.colors = {"blue": np.array([0, 0, 1]), "red": np.array([1, 0, 0])}
        self.plot = Plot((7.5, 6.8), 2, 3, self.colors)

    def test_init(self):
        self.assertAll(
            np.array(self.plot.figure.bbox_inches.bounds) == np.array([0.0, 0.0, 7.5, 6.8])
        )
        plot = Plot((2, 2), 1, 6)
        plot2 = Plot((2, 2))
        self.assertEqual(self.plot.axes.shape, (2, 3))
        self.assertEqual(len(plot.axes), 6)
        self.assertAll(np.array(self.plot.figure.axes) == np.array(self.plot.get_axes()))
        self.assertAll(np.array(plot.get_axes()) == np.array(plot.figure.axes))
        self.assertAll(np.array(plot2.get_axes()) == np.array(plot2.figure.axes))
        for i in range(6):
            self.assertEqual(self.plot.get_ax(i), self.plot.figure.axes[i])
            self.assertEqual(self.plot.get_ax(i), self.plot.figure.axes[i])
        self.assertEqual(plot2.get_ax(), plot2.figure.axes[0])
        with self.assertRaises(AssertionError):
            Plot((1, 1), dict_colors={"black": [0.0, 0.0]})

    def test_set_color(self):
        with self.assertRaises(AssertionError):
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
