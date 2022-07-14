from typing_extensions import Required
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config


@config.node
class ConnectomeGranuleGolgi(ConnectionStrategy):

    radius = config.attr(type=int, required=True)
    convergence = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] - c[0], 2)
                + np.power(chunk[1] - c[1], 2)
                + np.power(chunk[2] - c[2], 2)
            )
            if dist < self.radius and c[1] < chunk[1]:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        granule_pos = pre_ps.load_positions()
        golgi_pos = post_ps.load_positions()
        n_golgi = len(golgi_pos)
        n_granule = len(granule_pos)
        n_conn = n_golgi * n_granule
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        for i, golgi in enumerate(golgi_pos):
            # Compute the distance between the somata of the current golgi cell and granule cells
            dist = np.sqrt(
                np.pow(golgi[0] - granule_pos[0], 2)
                + np.pow(golgi[1] - granule_pos[1], 2)
                + np.pow(golgi[2] - granule_pos[2], 2)
            )
            # Find the granule cells whose distence from the soma of the golgi cell is less than the radius of the dendridic tree and such that they are under the golgi cell
            radius_condition = dist < self.radius
            y_condition = granule_pos[:, 1] < golgi[1]
            final_condition = np.logical_and(radius_condition, y_condition)
            # If less than 400 granule cells can be connected, take them all
            if np.count_nonzero(final_condition) < self.convergence:
                pre_idx = np.nonzero(final_condition)[0]
                post_locs[ptr : (ptr + len(pre_locs)), 0] = i
                pre_locs[ptr : (ptr + len(pre_locs)), 0] = pre_idx
                ptr += len(pre_idx)
            # Otherwise, take the first 400 occourrence setting the subsequent ones to False
            else:
                true_indices = np.argwhere(final_condition == True)
                final_condition[true_indices[len(true_indices) - 1] :] = False
                pre_idx = np.nonzero(final_condition)[0]
                post_locs[ptr : (ptr + len(pre_locs)), 0] = i
                pre_locs[ptr : (ptr + len(pre_locs)), 0] = pre_idx
                ptr += len(pre_idx)
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
