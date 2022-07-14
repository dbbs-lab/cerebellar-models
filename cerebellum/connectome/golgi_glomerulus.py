from zoneinfo import available_timezones
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config

@config.node
class ConnectomeGolgiGlomerulus(ConnectionStrategy):

    divergence = config.attr(type=int, required=True)
    radius = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        # print("CT", ct.name, "N", len(ct.get_placement_set()), "Placed in", len(ct.get_placement_set().get_all_chunks()), "chunks")
        selected_chunks = []
        for c in chunks:
            dist = np.sqrt(np.power(chunk[0] - c[0], 2) + np.power(chunk[2] - c[2], 2))
            if dist < self.radius and chunk[1] > c[1]:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        golgi_pos = pre_ps.load_positions()
        glom_pos = post_ps.load_positions()
        div = self.divergence
        n_golgi = len(golgi_pos)
        n_conn = n_golgi * div

        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        # Check the divergence
        num_glom = len(glom_pos < div)
        n_conn_goc = div
        if num_glom < div:
            print("The number of glomeruli in ROI is less than", div)
            n_conn_goc = num_glom

        avail = np.ones((n_conn_goc), dtype=bool)
        ptr = 0
        for i, gpos in enumerate(golgi_pos):
            dist = np.sqrt(
                np.power(gpos[0] - glom_pos[:, 0], 2)
                + np.power(gpos[1] - glom_pos[:, 1], 2)
                + np.power(gpos[2] - glom_pos[:, 2], 2)
            )
            cand = dist < self.radius
            cand = cand & avail
            avail = avail & ~cand
            pre_locs[i, 0] = i
            cand_idx = np.nonzero(cand)[0]
            post_locs[ptr : (ptr + len(cand_idx)), 0] = cand_idx
            pre_locs[ptr : (ptr + len(cand_idx)), 0] = i
            ptr += len(cand_idx)

        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
