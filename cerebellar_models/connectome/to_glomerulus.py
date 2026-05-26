"""
Module for the configuration node of every presynaptic cell to Glomerulus ConnectionStrategy
"""

import abc
import itertools

import numpy as np
from bsb import Chunk, ConnectionStrategy, InvertedRoI, config

from cerebellar_models.connectome.presyn_dist_strat import PresynDistStrat


def norm_exp_dist(size: int = 1, b: float = 2.0):
    """
    Normalized exponential random generator for distance based selection

    :param int size: number of random to sample
    :param float b: strength of the exponential decay.
    :return: random numbers sampled
    :rtype: numpy.ndarray
    """
    return (np.exp(-b * np.random.rand(size)) - np.exp(-b)) / (1 - np.exp(-b))


class ConnectomeGlomerulus(InvertedRoI, ConnectionStrategy):
    """
    BSB Connection strategy to connect a presynaptic cell to Glomeruli.
    """

    def connect(self, pre, post):
        # We use a truncated exponential distribution to favour the presynaptic fibers closer the
        # postsynaptic glomerulus.
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                presyn_pos = pre_ps.load_positions()
                glomeruli_pos = post_ps.load_positions()
                n_glom = len(glomeruli_pos)

                # We work assuming that there is at least 1 presynaptic fiber in the ROI.
                # Otherwise, there is something wrong in the placement phase.
                pre_locs = np.full((n_glom, 3), -1, dtype=int)
                post_locs = np.full((n_glom, 3), -1, dtype=int)

                # We connect each glomerulus to a presynaptic cell.
                for j, glomerulus in enumerate(glomeruli_pos):
                    pre_ids = self.pre_selection(presyn_pos, glomerulus)
                    if len(pre_ids) > 0:
                        pre_locs[j, 0] = pre_ids[int(np.floor(len(pre_ids) * norm_exp_dist()[0]))]
                    else:
                        # if there is no presyn within the box use the closest one
                        pre_locs[j, 0] = np.argmin(np.linalg.norm(presyn_pos - glomerulus, axis=1))
                    post_locs[j, 0] = j
                self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)

    @abc.abstractmethod
    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):  # pragma: nocover
        """
        Order presynaptic cell ids based on their respective distance to glomerulus

        :param numpy.ndarray presyn_pos: list of presynaptic cell positions
        :param numpy.ndarray glom_pos: single glomerulus position
        :return: presynaptic cell ids sorted by distance to glomerulus
        :rtype: numpy.ndarray
        """
        pass


@config.node
class ConnectomeMossyGlomerulus(ConnectomeGlomerulus):
    """
    BSB Connection strategy to connect Mossy fibers to Glomeruli.
    """

    x_length: float = config.attr(type=float, required=True)
    """Length of the box along the x axis surrounding the glomerulus cell soma in which the 
        presynaptic cell can be connected."""
    y_length: float = config.attr(type=float, required=True)
    """Length of the box along the y axis surrounding the glomerulus cell soma in which the 
        presynaptic cell can be connected."""

    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):
        diff = np.absolute(glom_pos - presyn_pos)
        ids_to_keep = np.where((diff[:, 0] <= self.x_length) * (diff[:, 1] <= self.y_length))[0]
        dist = np.linalg.norm(diff[ids_to_keep], axis=1)
        return ids_to_keep[np.argsort(dist)]


@config.node
class ConnectomeUBCGlomerulus(ConnectomeGlomerulus, PresynDistStrat):
    """
    BSB Connection strategy to connect UBC to Glomeruli.
    """

    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):
        dist = np.linalg.norm(presyn_pos - glom_pos, axis=1)
        ids_to_keep = np.where(dist <= self.radius)[0]
        return ids_to_keep[np.argsort(dist[ids_to_keep])]
