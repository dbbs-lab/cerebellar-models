import itertools

import numpy as np
from bsb import Chunk, InvertedRoI, config


class PresynDistStrat(InvertedRoI):
    radius = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
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
