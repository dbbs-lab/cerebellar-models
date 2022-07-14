import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config


@config.node
class ConnectomePFPurkinje(ConnectionStrategy):

    x_extension = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        # print("CT", ct.name, "N", len(ct.get_placement_set()), "Placed in", len(ct.get_placement_set().get_all_chunks()), "chunks")
        selected_chunks = []
        purkinje_x_max = self.x_extension / 2
        for c in chunks:
            x_pos = np.abs(chunk[0] - c[0]) * chunk.dimensions[0]
            if x_pos < purkinje_x_max:
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
        purkinje_pos = post_ps.load_positions()
        n_golgi = len(granule_pos)
        n_purkinje = len(purkinje_pos)
        n_conn = n_golgi * n_purkinje
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        for i, purkinje in enumerate(purkinje_pos):
            # Compute the distance on the x-axis between the somata of purkinje and granule cells
            x_dist = np.abs(purkinje[0] - granule_pos[0])
            # Find the granule cells whose x-distance from the current purkinje cell
            pre_gr = x_dist < self.x_extension
            pre_idx = np.nonzero(pre_gr)[0]
            post_locs[ptr : (ptr + len(pre_locs)), 0] = i
            pre_locs[ptr : (ptr + len(pre_locs)), 0] = pre_idx
            ptr += len(pre_idx)
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
