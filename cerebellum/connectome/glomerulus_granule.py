"""
    Module for the configuration node of the Glomerulus to Granule ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import (
    CellType,
    ConfigurationError,
    ConnectionStrategy,
    ConnectivityError,
    config,
    refs,
)

from cerebellum.connectome.presyn_dist_strat import PresynDistStrat


class TooFewGlomeruliClusters(ConnectivityError):
    """
    Error raised when too few glomerulus clusters are available for a postsynaptic cell.
    """

    pass


@config.node
class ConnectomeGlomerulusGranule(PresynDistStrat, ConnectionStrategy):
    """
    BSB Connection strategy to connect Glomerulus to Granule cells.
    With a convergence value set to `n`, this connection guarantees that each Granule cell connects
    to `n` unique Glomerulus clusters, where each Glomerulus cluster is connected to a different
    Mossy fiber.
    """

    convergence = config.attr(type=float, required=True)
    """Convergence value between Glomeruli and Granule cells. 
        Corresponds to the mean number of Glomeruli that has a single Granule cell as target"""
    mf_glom_strat: ConnectionStrategy = config.ref(refs.connectivity_ref, required=True)
    """Connection Strategy that links Mossy fibers to Glomeruli."""
    mf_cell_type: CellType = config.ref(refs.cell_type_ref, required=True)
    """Celltype used for the Mossy fibers."""

    @config.property
    def depends_on(self):
        # Get the possibly missing `_depends_on` list.
        # Fixme: does not work with str.
        deps = getattr(self, "_depends_on", None) or []
        # Strat is required, but depends on a reference that isn't available when the config loads.
        strat = getattr(self, "mf_glom_strat", None)
        if strat is None:
            return deps
        else:
            return [*{*deps, strat}]

    @depends_on.setter
    def depends_on(self, value):
        self._depends_on = value

    def _assert_dependencies(self):
        # assert dependency rule corresponds to mossy to glom
        post_ct = self.mf_glom_strat.postsynaptic.cell_types
        if len(post_ct) != 1 or post_ct[0] not in self.presynaptic.cell_types:
            raise ConfigurationError(
                "Postsynaptic cell of dependency rule does not match this rule's"
                " presynaptic cell."
            )

    def boot(self):
        self._assert_dependencies()

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _get_mf_clusters(self, pre_ps):
        # Find the glomeruli clusters

        cs = self.mf_glom_strat.get_output_names(self.mf_cell_type, pre_ps.cell_type)
        assert (
            len(cs) == 1
        ), f"Only one connection set should be given from {self.mf_glom_strat.name}."
        cs = self.scaffold.get_connectivity_set(cs[0])
        # find mf-glom connections where the postsyn chunk corresponds to the
        # glom-grc presyn chunk
        iter = cs.load_connections().to(pre_ps.get_loaded_chunks())
        mf_locs, glom_locs = iter.all()

        unique_mossy = np.unique(mf_locs[:, 0])
        if unique_mossy.size < self.convergence:
            raise TooFewGlomeruliClusters(
                "Less than 4 unique mossy fibers have been found. "
                "Check the densities of mossy fibers and glomeruli in "
                "the configuration file."
            )

        clusters = []
        for current in unique_mossy:
            glom_idx = np.where(mf_locs[:, 0] == current)[0]
            clusters.append(glom_locs[glom_idx, 0])

        return unique_mossy, clusters

    def _connect_type(self, pre_ps, post_ps):
        glom_pos = pre_ps.load_positions()
        gran_pos = post_ps.load_positions()

        # Find the glomeruli clusters
        unique_mossy, clusters = self._get_mf_clusters(pre_ps)

        gran_morphos = post_ps.load_morphologies().iter_morphologies(cache=True, hard_cache=True)

        n_conn = int(np.round(len(gran_pos) * self.convergence))
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        for i, gr_pos, morpho in zip(itertools.count(), gran_pos, gran_morphos):
            # morpho should have enough dendrites to match convergence
            dendrites = morpho.get_branches()
            assert len(dendrites) >= self.convergence

            # Randomize the order of the clusters and dendrites
            cluster_idx = np.arange(0, unique_mossy.size)
            np.random.shuffle(cluster_idx)
            dendrites_idx = np.arange(0, len(dendrites))
            np.random.shuffle(dendrites_idx)

            # The following loop connects a glomerulus from each cluster to a grc dendrite
            # until the convergence rule is reached.
            # First it checks for glomerulus that are close enough,
            # otherwise, it chooses the closest glomerulus from each remaining cluster.
            gr_connections = 0
            current_cluster = 0
            check_dist = True
            while gr_connections < self.convergence:
                if current_cluster >= len(cluster_idx):
                    # Not enough glom were found close enough to the GrC.
                    # Select from the remaining (more distant) gloms.
                    current_cluster = 0
                    check_dist = False
                nc = cluster_idx[current_cluster]
                dist = np.linalg.norm(gr_pos - glom_pos[clusters[nc]], axis=1)
                if check_dist:
                    # Try to select a cell from 4 clusters satisfying the conditions
                    close_indices = np.nonzero(dist < self.radius)[0]
                    if len(close_indices) == 0:
                        current_cluster += 1
                        continue
                    # Id of the glomerulus, randomly selected between the available ones
                    rnd = np.random.randint(low=0, high=len(close_indices))
                    id_glom = clusters[nc][close_indices[rnd]]
                else:
                    # If there are some free dendrites, connect them to the closest glomeruli,
                    # even if they do not satisfy the geometric conditions.
                    # Id of the glomerulus, randomly selected between the available ones
                    id_glom = clusters[nc][np.argmin(dist)]
                pre_locs[ptr + gr_connections, 0] = id_glom
                # Id of the granule cell
                post_locs[ptr + gr_connections, 0] = i
                # Select one of the 4 dendrites
                dendrite = dendrites[dendrites_idx[gr_connections]]
                post_locs[ptr + gr_connections, 1] = morpho.branches.index(dendrite)
                # Select the terminal point of the branch
                post_locs[ptr + gr_connections, 2] = len(dendrite) - 1

                gr_connections += 1
                # remove cluster used already
                cluster_idx = np.delete(cluster_idx, current_cluster)
            ptr += gr_connections

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
