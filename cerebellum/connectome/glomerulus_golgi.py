import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from scipy.stats.distributions import truncexpon
from bsb.morphologies import Morphology


@config.node
class ConnectomeGlomerulusGolgi(ConnectionStrategy):
    radius = config.attr(type=int, required=True)

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
        
        glomeruli_pos = pre_ps.load_positions()
        golgi_pos = post_ps.load_positions()
        n_glom = len(glomeruli_pos)
        n_golgi = len(golgi_pos)
        n_conn = n_glom * n_golgi
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        # Find cells to connect
        ptr = 0
        for i, golgi in enumerate(golgi_pos):
            dist = np.sqrt(
                np.power(golgi[0] - glomeruli_pos[:, 0], 2)
                + np.power(golgi[1] - glomeruli_pos[:, 1], 2)
                + np.power(golgi[2] - glomeruli_pos[:, 2], 2)
            )
            to_connect_bool = dist < self.radius
            to_connect_idx = np.nonzero(to_connect_bool)[0]
            connected_gloms = len(to_connect_idx)
            pre_locs[ptr : (ptr + connected_gloms), 0] = to_connect_idx
            post_locs[ptr : (ptr + connected_gloms), 0] = i

            # Find which dendrite to connect
            basal_dendrides_branches = post_ct.morphology.get_branches(["basal_dendrites"])
            # We keep track of the index of a branch and of the index of a point on a branch
            id_pairs = np.full((len(post_ct.morphology.points), 2), -1, dtype=int)
            basal_points_coord = np.full((len(post_ct.morphology.points), 3), -1, dtype=float)
            ptr_idx = 0
            for idx_b,branch in enumerate(post_ct.morphology.get_branches(basal_dendrides_branches)):
                for idx_p,coordinates in enumerate(branch.points):
                    id_pairs[ptr_idx] = [idx_b,idx_p]
                    basal_points_coord[ptr_idx] = coordinates
                    ptr_idx += 1
            # Draw rolls from the exponential distribution equal to the total amount
            # of synaptic contacts to be made between this Golgi cell and all its
            # glomeruli.
            num_basal_points = len(basal_points_coord)
            rolls = exp_dist.rvs(size=num_basal_points)
            # Compute the distance between terminal points of basal dendrites 
            # and the soma of the avaiable glomeruli
            for id_g,glom_p in enumerate(glomeruli_pos):
                for pt_coord in basal_points_coord:
                    pts_dist = np.sqrt(
                        np.power(basal_points_coord[:,0] + golgi[0] - glom_p[0], 2)
                        + np.power(basal_points_coord[:,1] + golgi[1] - glom_p[1], 2)
                        + np.power(basal_points_coord[:,2] + golgi[2] - glom_p[2], 2)
                    )
                sorted_pts_ids = np.argsort(pts_dist)
                # Pick the point in which we form a synapse according to a exponential distribution mapped
                # through the distance indices: high chance to pick closeby points.
                pt_idx = [max(int(k * num_basal_points)) for k in rolls]
                post_locs[ptr+i,1] = id_pairs[pt_idx,0]
                post_locs[ptr+i,2] = id_pairs[pt_idx,1]
                
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
