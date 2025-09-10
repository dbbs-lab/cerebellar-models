import unittest

import numpy as np
from bsb_test import NumpyTestCase

from cerebellar_models.placement import utils


class TestUtils(unittest.TestCase, NumpyTestCase):

    def test_modulo(self):
        self.assertAll(
            utils.signed_modulo(np.array([0, -5, 50, 14, -22, 13, 2]), 5)
            == np.array([0, 0, 0, -1, -2, -2, 2])
        )

    def test_line(self):
        self.assertAll(
            utils.bresenham_line(np.array([0, -5, 50]), np.array([0, -5, 50]))
            == np.array([[0, -5, 50]])
        )
        expected = np.zeros((6, 3), dtype=int)
        expected[:3, 1] = -5
        expected[3:, 1] = -4
        expected[3:, 0] = -1
        expected[:, 2] = np.arange(6)
        self.assertAll(
            utils.bresenham_line(np.array([0, -5, 0]), np.array([-1, -4, 5])) == expected
        )
        expected[:, 1] = np.arange(-5, 1, 1)
        self.assertAll(utils.bresenham_line(np.array([0, -5, 0]), np.array([-1, 0, 5])) == expected)

    def test_boundaries_index_of(self):
        self.assertEqual(
            utils.boundaries_index_of(np.array([12, -4, 5]), np.array([-1, -4, 7])), (0, 1, 2)
        )
