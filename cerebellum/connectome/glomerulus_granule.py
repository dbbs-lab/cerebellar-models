"""
    Module for the configuration node of the Glomerulus to Granule ConnectionStrategy
"""

import functools
import itertools

import numpy as np
from bsb import (
    ConfigurationError,
    ConnectionStrategy,
    ConnectivityError,
    InvertedRoI,
    config,
    refs,
)


class TooFewGlomeruliClusters(ConnectivityError):
    """
    Error raised when too few glomerulus clusters are available for a postsynaptic cell.
    """

    pass


@config.node
class ConnectomeGlomerulusGranule(InvertedRoI, ConnectionStrategy):
    """
    BSB Connection strategy to connect Glomerulus to Granule cells.
    With a convergence value set to `n`, this connection guarantees that each Granule cell connects
    to `n` unique Glomerulus clusters, where each Glomerulus cluster is connected to a different
    Mossy fiber.
    """

    radius = config.attr(type=int, required=True)
    """Radius of the sphere to filter the presynaptic chunks within it."""
    convergence: float = config.attr(type=float, required=True)
    """Convergence value between Glomeruli and Granule cells. 
        Corresponds to the mean number of Glomeruli that has a single Granule cell as target"""
    pre_glom_strats = config.reflist(refs.connectivity_ref, required=True)
    """Connection Strategies that links Pre-presyn cell to Glomeruli."""
    pre_cell_types = config.reflist(refs.cell_type_ref, required=True)
    """Celltype used for the pre-presyn cell."""

    @config.property
    def depends_on(self):
        # Get the possibly missing `_depends_on` list.
        # Fixme: does not work with str.
        deps = getattr(self, "_depends_on", None) or []
        # Strat is required, but depends on a reference that isn't available when the config loads.
        strats = getattr(self, "pre_glom_strats", None)
        if strats is None:
            return deps
        else:
            return [*{*deps, *strats}]

    @depends_on.setter
    def depends_on(self, value):
        self._depends_on = value

    def _assert_dependencies(self):
        # assert dependency rule corresponds to mossy to glom
        for strat in self.pre_glom_strats:
            post_ct = strat.postsynaptic.cell_types
            if len(post_ct) != 1 or post_ct[0] not in self.presynaptic.cell_types:
                raise ConfigurationError(
                    "Postsynaptic cell of dependency rule does not match this rule's"
                    " presynaptic cell."
                )
        for strat in self.pre_glom_strats:
            pre_ct = strat.presynaptic.cell_types
            found = False
            for pre_pre_ct in self.pre_cell_types:
                if pre_pre_ct in pre_ct:
                    found = True
                    break
            if not found:
                raise ConfigurationError(
                    f"Presynaptic cells of dependency rule {strat} does not contain any of the provided presynaptic "
                    f"cell types: {[ct.name for ct in pre_ct]}."
                )

    def boot(self):
        self._assert_dependencies()

    def connect(self, pre, post):
        for post_ps in post.placement:
            self._connect_type(pre, post_ps)

    @functools.cache
    def load_connections(self):
        dict_cs = {}
        for pre_ct in self.presynaptic.cell_types:
            for strat in self.pre_glom_strats:
                for pre_pre_ct in self.pre_cell_types:
                    try:
                        cs = strat.get_output_names(pre_pre_ct, pre_ct)
                    except ValueError:
                        continue
                    assert (
                        len(cs) == 1
                    ), f"Only one connection set should be given from {strat.name} with type {pre_pre_ct.name}."
                    dict_cs[cs[0]] = list(
                        self.scaffold.get_connectivity_set(cs[0]).load_connections().all()
                    )
        return dict_cs

    def clean_cache(self, scaffold):
        self.load_connections.cache_clear()

    def queue(self, pool):
        super().queue(pool)
        pool.queue(self.clean_cache, deps=pool.get_submissions_of(self))

    def _get_pre_clusters(self, pre_ps):
        # Find the glomeruli clusters

        clusters = []
        unique_pres = []
        for strat in self.pre_glom_strats:
            for pre_ct in self.pre_cell_types:
                try:
                    cs = strat.get_output_names(pre_ct, pre_ps.cell_type)
                except ValueError:
                    continue
                assert (
                    len(cs) == 1
                ), f"Only one connection set should be given from {strat.name} with type {pre_ct.name}."
                # find pre-glom connections where the postsyn chunk corresponds to the
                # glom-grc presyn chunk
                pre_locs, glom_locs = self.load_connections()[cs[0]]
                ct_uniques = np.unique(pre_locs[:, 0])
                unique_pres.extend(ct_uniques)

                for current in ct_uniques:
                    glom_idx = np.where(pre_locs[:, 0] == current)[0]
                    clusters.append(glom_locs[glom_idx, 0])

        return unique_pres, clusters

    def _connect_type(self, pre, post_ps):
        gran_pos = post_ps.load_positions()
        gran_morphos = post_ps.load_morphologies().iter_morphologies(cache=True, hard_cache=True)

        # Find the glomeruli clusters
        class presyn_dict:
            def __init__(cls, glom_pos, unique_pre, clusters):
                cls.glom_pos = glom_pos
                cls.unique_pre = unique_pre
                cls.clusters = clusters

        presyn_dicts = [
            presyn_dict(pre_ps.load_positions(), *self._get_pre_clusters(pre_ps))
            for pre_ps in pre.placement
        ]
        unique_pre = np.concatenate([np.arange(len(d.unique_pre)) for d in presyn_dicts])
        dict_ids = np.repeat(
            np.arange(len(presyn_dicts)), [len(d.unique_pre) for d in presyn_dicts]
        )
        if len(unique_pre) < self.convergence:
            raise TooFewGlomeruliClusters(
                "Less than 4 unique pre-presynaptic cells have been found. "
                "Check the densities of pre-presynaptic cells and glomeruli in "
                "the configuration file."
            )

        # TODO: implement random rounding and adapt tests.
        n_conn = int(np.round(len(gran_pos) * self.convergence))
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        selected_pre = np.full(n_conn, -1, dtype=int)
        ptr = 0
        for i, gr_pos, morpho in zip(itertools.count(), gran_pos, gran_morphos):
            # morpho should have enough dendrites to match convergence
            dendrites = morpho.get_branches()
            assert len(dendrites) >= self.convergence

            # Randomize the order of the clusters and dendrites
            cluster_idx = np.arange(0, len(unique_pre))
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
                loc_dict = presyn_dicts[dict_ids[nc]]
                clusters = loc_dict.clusters[unique_pre[nc]]
                dist = np.linalg.norm(gr_pos - loc_dict.glom_pos[clusters], axis=1)
                if check_dist:
                    # Try to select a cell from 4 clusters satisfying the conditions
                    close_indices = np.nonzero(dist < self.radius)[0]
                    if len(close_indices) == 0:
                        current_cluster += 1
                        continue
                    # Id of the glomerulus, randomly selected between the available ones
                    rnd = np.random.randint(low=0, high=len(close_indices))
                    id_glom = clusters[close_indices[rnd]]
                else:
                    # If there are some free dendrites, connect them to the closest glomeruli,
                    # even if they do not satisfy the geometric conditions.
                    # Id of the glomerulus, randomly selected between the available ones
                    id_glom = clusters[np.argmin(dist)]
                selected_pre[ptr + gr_connections] = dict_ids[nc]
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
        for i, pre_ps in enumerate(pre.placement):
            self.connect_cells(
                pre_ps, post_ps, pre_locs[selected_pre == i], post_locs[selected_pre == i]
            )
