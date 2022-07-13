import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config


@config.node
class ConnectomeAscAxonPurkinje(ConnectionStrategy):
    divergence = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        roi_xs = np.arange(chunk[0] - 1, chunk[0] + 1, 1)
        roi_ys = np.arange(chunk[1] - 1, chunk[1] + 1, 1)
        roi_zs = np.arange(chunk[2] - 1, chunk[2] + 1, 1)
        print("Chunk type?", type(roi_zs))
        my_roi_coord = np.meshgrid(roi_xs, roi_ys, roi_zs)
        print("And after meshgrid?", type(my_roi))
        my_roi = []
        for i in my_roi_coord:
            my_roi.append(
                Chunk(
                    [my_roi_coord[0], my_roi_coord[1], my_roi_coord[2]],
                    chunk.dimensions,
                )
            )
        return my_roi

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        aa_pos = pre_ps.load_positions()
        purkinje_pos = post_ps.load_positions()
        num_aa = len(aa_pos)
        num_purkinje = len(purkinje_pos)
        n_conn = num_aa * num_purkinje
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        for i, purkinje in enumerate(purkinje_pos):
            # Compute the x-distance and the z-distance between the somata of the current purkinje cell and the granule cells
            x_dist = np.abs(purkinje[0] - aa_pos[0])
            z_dist = np.abs(purkinje[2] - aa_pos[2])
            # Check if the conditions on x and z are satisfied, then compute the logical product elementwise
            x_condition = x_dist < 130 / 2
            z_condition = z_dist < 130 / 2
            final_cond = np.logical_and(x_condition, z_condition)
            # Find the indices of the aa to connect
            pre_idx = np.nonzero(final_cond)[0]
            post_locs[ptr : (ptr + len(pre_locs)), 0] = i
            pre_locs[ptr : (ptr + len(pre_locs)), 0] = pre_idx
            # Increment ptr
            ptr += len(pre_idx)
        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
