"""
    Module for the configuration node of every presynaptic cell to Glomerulus ConnectionStrategy
"""

import abc
import itertools

import numpy as np
from bsb import Chunk, ConnectionStrategy, InvertedRoI, config


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
                to_keep = 0
                for j, glomerulus in enumerate(glomeruli_pos):
                    pre_ids = self.pre_selection(presyn_pos, glomerulus)
                    if len(pre_ids) > 0:
                        roll = int(np.floor(len(pre_ids) * norm_exp_dist()))
                        pre_locs[to_keep, 0] = pre_ids[roll]
                        post_locs[to_keep, 0] = j
                        to_keep += 1

                self.connect_cells(pre_ps, post_ps, pre_locs[:to_keep], post_locs[:to_keep])

    @abc.abstractmethod
    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):  # pragma: no cover
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
    z_length: float = config.attr(type=float, required=True)
    """Length of the box along the z axis surrounding the glomerulus cell soma in which the 
        presynaptic cell can be connected."""

    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):
        diff = np.absolute(glom_pos - presyn_pos)
        ids_to_keep = np.where((diff[:, 0] <= self.x_length) * (diff[:, 2] <= self.z_length))[0]
        dist = np.linalg.norm(diff[ids_to_keep], axis=1)
        return ids_to_keep[np.argsort(dist)]

    def get_region_of_interest(self, chunk):
        # Chunk here is a postsynaptic chunk because of InvertedRoI
        # We look for chunks containing mossy fibers that are within a rectangle of size
        # x_length * z_length centered on the postsynaptic chunk containing the glomerulus.
        chunks = set(
            itertools.chain.from_iterable(
                ct.get_placement_set().get_all_chunks() for ct in self.presynaptic.cell_types
            )
        )
        selected_chunks = []
        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[2] - c[2])
            x_dist = x_dist * chunk.dimensions[0]
            z_dist = z_dist * chunk.dimensions[2]

            if (x_dist < self.x_length / 2) and (z_dist < self.z_length / 2):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks
