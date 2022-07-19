import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from scipy.stats.distributions import truncexpon


@config.node
class ConnectomeGlomerulusGolgi(ConnectionStrategy):
    radius = config.attr(type=int, required=True)
    detailed = config.attr(type=bool, required=True)
    # Is it the correct type?
    compartments = config.attr(type=list, required=True)

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
        # Look for chunks which are less than radius away
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] * chunk.dimensions[0] - c[0] * c.dimensions[0], 2)
                + np.power(chunk[1] * chunk.dimensions[1] - c[1] * c.dimensions[0], 2)
                + np.power(chunk[2] * chunk.dimensions[2] - c[2] * c.dimensions[0], 2)
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

        # If synaptic contacts need to be made we use this exponential distribution
        # to pick the closer by compartments.
        exp_dist = truncexpon(b=5, scale=0.03)

        glom_pos = pre_ps.load_positions()
        golgi_pos = post_ps.load_positions()
        n_glom = len(glom_pos)
        n_golgi = len(golgi_pos)
        n_conn = n_glom * n_golgi
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        compartments = np.zeros((0, 1))

        ptr = 0
        for i, golgi in enumerate(golgi_pos):
            dist = np.sqrt(
                np.power(golgi[0] - glom_pos[:, 0], 2)
                + np.power(golgi[1] - glom_pos[:, 1], 2)
                + np.power(golgi[2] - glom_pos[:, 2], 2)
            )
            to_connect_bool = dist < self.radius
            to_connect_idx = np.nonzero(to_connect_bool)[0]
            connected_gloms = len(to_connect_idx)
            pre_locs[ptr : (prt + connected_gloms), 0] = to_connect_idx
            post_locs[ptr : (prt + connected_gloms), 0] = i

            # If GOC is detailed, find the closest compartments
            if self.detailed:
                rolls = exp_dist.rvs(size=connected_gloms)
                gg_comps = np.zeros((connected_gloms, 1))
                comps = selfs.compartments
                total_compartments = len(comps)

                # Does c.start[0] give the x position relative to the Golgi soma position?
                for c in comps:
                    c_dist = np.sqrt(
                        np.power(c.start[0] + golgi[0] - glom_pos[:, 0], 2)
                        + np.power(c.start[1] + golgi[1] - glom_pos[:, 1], 2)
                        + np.power(c.start[2] + golgi[2] - glom_pos[:, 2], 2)
                    )
                    d_comps = np.argsort(c_dist)

                    # Pick compartments according to a exponential distribution mapped
                    # through the distance indices: high chance to pick closeby comps.
                    gg_comps = [
                        comps[d_comps[int(k * total_compartments)]].id for k in rolls
                    ]
                    compartments = np.vstack((compartments, gg_comps))

            prt += connected_gloms

        if self.detailed:
            self.scaffold.connect_cells(
                self,
                pre_ps,
                post_ps,
                pre_locs[:ptr],
                post_locs[:ptr],
                compartments=compartments,
            )
        else:
            self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
