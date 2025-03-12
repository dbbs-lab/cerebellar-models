"""
Module for the configuration node of the Glomerulus to Golgi ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import ConnectionStrategy, ConnectivityError, config
from scipy.stats.distributions import truncexpon

from cerebellum.connectome.presyn_dist_strat import PresynDistStrat


@config.node
class ConnectomeGlomerulusGolgi(PresynDistStrat, ConnectionStrategy):
    """
    BSB Connection strategy to connect Glomerulus to Golgi cells.
    """

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _connect_type(self, pre_ps, post_ps):
        glomeruli_pos = pre_ps.load_positions()
        golgi_pos = post_ps.load_positions()

        # If synaptic contacts need to be made we use this exponential distribution
        # to pick the closer by subcell_labels.
        exp_dist = truncexpon(b=5, scale=0.03)
        golgi_morphos = post_ps.load_morphologies().iter_morphologies(cache=True, hard_cache=True)

        n_conn = len(glomeruli_pos) * len(golgi_pos)
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        for i, golgi, morpho in zip(itertools.count(), golgi_pos, golgi_morphos):
            to_connect_idx = np.nonzero(
                np.linalg.norm(golgi - glomeruli_pos, axis=1) < self.radius
            )[0]
            connected_gloms = len(to_connect_idx)
            pre_locs[ptr : ptr + connected_gloms, 0] = to_connect_idx
            post_locs[ptr : ptr + connected_gloms, 0] = i

            # Find terminal points on branches
            basal_dendrides_branches = morpho.get_branches()
            terminal_branches_ids = np.nonzero([b.is_terminal for b in basal_dendrides_branches])[0]
            basal_dendrides_branches = np.take(
                basal_dendrides_branches, terminal_branches_ids, axis=0
            )
            if basal_dendrides_branches.size == 0:
                raise ConnectivityError(
                    "The golgi morphology provided has no terminal branches.\n"
                    "Check the morphology_labels."
                )

            # Find the point-on-branch ids of the tips
            tips_coordinates = np.array([b.points[-1] for b in basal_dendrides_branches])

            # Connect each close by glom to a tip of the golgi
            for id_g, glom_p in enumerate(to_connect_idx):
                sorted_pts_ids = np.argsort(
                    np.linalg.norm(tips_coordinates + golgi - glom_p, axis=1)
                )
                # Pick the golgi tip according to an exponential distribution mapped
                # through the distance to each glom: high chance to pick close by tips.
                pt_idx = sorted_pts_ids[tips_coordinates.size * exp_dist.rvs(size=1).astype(int)[0]]

                post_locs[ptr + id_g, 1] = morpho.branches.index(basal_dendrides_branches[pt_idx])
                post_locs[ptr + id_g, 2] = len(basal_dendrides_branches[pt_idx]) - 1
            ptr += connected_gloms

        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
