import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config

@config.node
class ConnectomeGlomerulusGranule(ConnectionStrategy):
    radius = config.attr(type=int, required=True)
    divergence = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        print(
            "CT",
            ct.name,
            "N",
            len(ct.get_placement_set()),
            "Placed in",
            len(ct.get_placement_set().get_all_chunks()),
            "chunks",
        )
        selected_chunks = []
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] - c[0], 2)
                + np.power(chunk[1] - c[1], 2)
                + np.power(chunk[2] - c[2], 2)
            )
            if dist < self.radius:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        glom_pos = pre_ps.load_positions()
        gran_pos = post_ps.load_positions()
        n_glom = len(glom_pos)
        n_gran = len(gran_pos)
        max_connections = self.divergence
        n_conn = n_glom * max_connections
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        ptr = 0
        for i, grpos in enumerate(gran_pos):
            gr_connection = 0
            dist = np.sqrt(
                np.power(grpos[0] - glom_pos[:, 0], 2)
                + np.power(grpos[1] - glom_pos[:, 1], 2)
                + np.power(grpos[2] - glom_pos[:, 2], 2)
            )
            avail = np.ones((n_glom), dtype=bool)
            sorted_indices = np.argsort(dist)
            sorted_dist = np.sort(dist)
            # cand = sorted_dist < 40
            # Sort the distances array
            for j, gdist in enumerate(sorted_dist):
                if gr_connection < self.divergence:
                    if gdist < self.radius:
                        post_locs[ptr + j, 0] = i
                        pre_locs[ptr + j, 0] = sorted_indices[j]
                        gr_connection += 1
                else:
                    break
            ptr += gr_connection
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
