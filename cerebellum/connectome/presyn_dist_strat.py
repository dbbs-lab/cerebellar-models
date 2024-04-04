"""
    Module for the utility class for postsynaptically-sorted ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import Chunk, InvertedRoI, config


class PresynDistStrat(InvertedRoI):
    """
    Mixin class that id used for ConnectionStrategy that deal with the connections for a pre- and
    post-synaptic pair sorting them by the post-synaptic cell chunk.
    """

    radius = config.attr(type=int, required=True)
    """Radius of the sphere to filter the presynaptic chunks within it."""

    def get_region_of_interest(self, chunk):
        """
        Finds all the presynaptic chunks that are within a sphere of defined radius, centered on the
        postsynaptic chunk.

        :param chunk:
        :type chunk: bsb.Chunk
        :return: list of presynaptic chunks
        :rtype: list[bsb.Chunk]
        """

        # Fixme: Distance between chunk is done corner to corner. It might not detect all chunks #34
        chunks = set(
            itertools.chain.from_iterable(
                ct.get_placement_set().get_all_chunks() for ct in self.presynaptic.cell_types
            )
        )
        # Look for chunks which are less than radius away from the current one.
        selected_chunks = []
        for c in chunks:
            if np.linalg.norm(chunk * chunk.dimensions - c * c.dimensions) <= self.radius:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks
