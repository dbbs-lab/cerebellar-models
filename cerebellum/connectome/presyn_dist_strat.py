"""
Module for the utility class for postsynaptically-sorted ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import BoxTree, Chunk, InvertedRoI, config


def get_close_chunks(chunk, target_chunks, radius):
    """
    Look for target chunks which are less than radius away from the current one.

    :param chunk: Source chunk.
    :type chunk: bsb.storage._chunks.Chunk
    :param target_chunks: Target chunk.
    :type target_chunks: Set[bsb.storage._chunks.Chunk]
    :radius: Maximum distance from the source chunk.
    :return: list of presynaptic chunks
    :rtype: list[bsb.storage._chunks.Chunk]
    """
    postsyn_chunk = [
        np.concatenate((chunk * chunk.dimensions - radius, (chunk + 1) * chunk.dimensions + radius))
    ]
    chunks_mbb = [np.concatenate((c, c + 1) * c.dimensions) for c in target_chunks]
    selected = np.array([i for i in BoxTree(chunks_mbb).query(postsyn_chunk)])
    return [c for i, c in enumerate(target_chunks) if i in selected]


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

        :param chunk: Postsynaptic chunk.
        :type chunk: bsb.storage._chunks.Chunk
        :return: list of presynaptic chunks
        :rtype: list[bsb.storage._chunks.Chunk]
        """

        chunks = set(
            itertools.chain.from_iterable(
                ct.get_placement_set().get_all_chunks() for ct in self.presynaptic.cell_types
            )
        )
        return get_close_chunks(chunk, chunks, self.radius)
