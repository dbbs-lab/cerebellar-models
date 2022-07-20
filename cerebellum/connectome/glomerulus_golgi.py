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
        # Find Golgi cells to connect
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
            '''
            # Find which dendrite to connect
            basal_dendrides = post_ct.morphology.get_branches(["basal_dendrites"])
            basal_points = []
            for basal_branch in basal_dendrides:
                if basal_branch.is_terminal:
                    basal_points.append(
                        basal_branch.get_points_labelled("basal_dendrites")
                    )
            num_basal_points = len(basal_points)
            rolls = exp_dist.rvs(size=connected_gloms)
            # Compute the distance between terminal points of basal dendrites 
            # and the soma of the avaiable glomeruli
            for b_point in basal_points:
                bp_dist = np.sqrt(
                    np.power(b_point[0] + golgi[0] - glomeruli_pos[:, 0], 2)
                    + np.power(b_point[1] + golgi[1] - glomeruli_pos[:, 1], 2)
                    + np.power(b_point[2] + golgi[2] - glomeruli_pos[:, 2], 2)
                )
            # To employ vectorized/parallelized operations in numpy, we select which dendride connects
            # with each one of the glomeruli in post, then we take only the connections satisfying
            # the geometric constraints using an elementwise multiplication and taking only
            # the non-zero elements
            synapse_points = [basal_points[int(k * num_basal_points)].id for k in rolls]
            synapse_points = int(to_connect_bool)*synapse_points
            post_locs[ptr + connected_gloms), 1] = index
            '''
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])