"""
Implementation of the BSB framework for cerebellar cortex reconstructions and simulations.
"""

import os

__version__ = "0.0.1"

import numpy as np


def templates():
    """
    :meta private:
    """
    return [os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))]


def _swap_morphology_axes(morphology, axis1: int, axis2: int):
    """
    Interchange two axes of a morphology points.

    :param bsb.morphologies.Morphology morphology: Morphology to modify
    :param int axis1: index of the first axis to exchange
    :param int axis2: index of the second axis to exchange
    :return: the modified morphology
    :rtype: bsb.morphologies.Morphology
    """
    for b in morphology.branches:
        old_column = np.copy(b.points[:, axis1])
        b.points[:, axis1] = b.points[:, axis2]
        b.points[:, axis2] = old_column

    return morphology
