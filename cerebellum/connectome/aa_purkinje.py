import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config


@config.node
class ConnectomeAscAxonPurkinje(ConnectionStrategy):
    """
    Legacy implementation for the connections between ascending axons and purkinje cells.
    """

    divergence = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        roi_xs = np.arange(chunk[0] - 1, chunk[0] + 1, 1)
        roi_ys = np.arange(chunk[1] - 1, chunk[1] + 1, 1)
        roi_zs = np.arange(chunk[2] - 1, chunk[2] + 1, 1)
        print("Chunk type?", type(roi_zs))
        my_roi = np.meshgrid(roi_xs, roi_ys, roi_zs)
        print("And after meshgrid?", type(my_roi))
        return my_roi

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        # N x 3 (xyz)
        pre_pos = pre_ps.load_positions()
        post_pos = post_ps.load_positions()

        # numpy magic

        # cell id, branch id, point id
        # cell id, -1, -1
        pre_locs = np.zeros((N, 3))
        post_locs = np.zeros((N, 3))

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
