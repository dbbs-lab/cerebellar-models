"""
    Module for the configuration node of the Golgi to Glomerulus ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import ConfigurationError, ConnectionStrategy, config, refs
from connectome.presyn_dist_strat import get_close_chunks


@config.node
class ConnectomeGolgiGlomerulus(ConnectionStrategy):
    """
    BSB Connection strategy to connect Golgi cells to postsynaptic cells through Glomeruli.
    With a divergence value set to `n`, this connection guarantees that each golgi cell
    connects to all postsynaptic cells that are themselves connected to `n` unique Glomerulus.
    """

    divergence: float = config.attr(type=float, required=True)
    """Divergence value between Golgi cells and Glomeruli. 
        Corresponds to the mean number of Glomeruli targeted by a single Golgi cell"""
    radius: float = config.attr(type=float, required=True)
    """Radius of the sphere surrounding the Golgi cell soma in which glomeruli can be connected."""
    glom_post_strats = config.reflist(refs.connectivity_ref, required=True)
    """Connection Strategies that links Glomeruli to the postsynaptic cells."""
    glom_cell_types = config.reflist(refs.cell_type_ref, required=True)
    """Cell types used for the Glomeruli."""

    @config.property
    def depends_on(self):
        # Get the possibly missing `_depends_on` list.
        deps = getattr(self, "_depends_on", None) or []
        # Strat is required, but depends on a reference that isn't available when the config loads.
        strat = getattr(self, "glom_post_strats", None)
        if strat is None:
            return deps
        else:
            return [*{*deps, *strat}]

    @depends_on.setter
    def depends_on(self, value):
        self._depends_on = value

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            self._connect_type(pre_ps, post)

    def _assert_dependencies(self):
        # assert dependency rule corresponds to glom to post
        for glom_post_strat in self.glom_post_strats:
            # At least one glomerulus type should match for each strat.
            found = False
            for glom_cell_type in self.glom_cell_types:
                if found:
                    break
                for pre_ct in glom_post_strat.presynaptic.cell_types:
                    if pre_ct == glom_cell_type:
                        found = True
                        break
            if not found:
                raise ConfigurationError(
                    f"Presynaptic cell of dependency rule {glom_post_strat.name} should match "
                    f"at least one of the provided glom_cell_types: "
                    f"{[glom_cell_type.name for glom_cell_type in self.glom_cell_types]}."
                )

            post_ct = glom_post_strat.postsynaptic.cell_types
            found = False
            for ct in self.postsynaptic.cell_types:
                if ct in post_ct:
                    found = True
                    break
            if not found:
                raise ConfigurationError(
                    f"The dependency rule {glom_post_strat.name} does not connect glomeruli to this connection's "
                    f"postsynaptic cell: {post_ct.name}."
                )

    def boot(self):
        self._assert_dependencies()

    def _combine_unique_gloms(self, glom_ids, loc_glom_pos, loc_post_locs):
        # Extract unique glomerulus groups
        unique_gloms = np.unique(glom_ids)
        loc_post_locs = np.asarray(loc_post_locs, dtype=object)
        post_locs = []
        glom_pos = np.zeros((len(unique_gloms), 3), dtype=float)
        for c, u_glom in enumerate(unique_gloms):
            ids = np.where(glom_ids == u_glom)[0]
            post_locs.append(
                [locs for locs in np.asarray(np.concatenate(loc_post_locs[ids]), dtype=int)]
            )
            glom_pos[c] = loc_glom_pos[ids[0]]

        return (unique_gloms, glom_pos, post_locs)

    def _get_glom_cluster(self, chunk, post_ps, glom_type):
        # Get the glom_to_post connections
        glom_ids = []
        post_locs = []
        loc_glom_pos = np.empty((0, 3))
        for glom_post_strat in self.glom_post_strats:
            try:
                # Extract connection sets that link glom_type to the postsynaptic cell type.
                cs = glom_post_strat.get_output_names(glom_type, post_ps.cell_type)
            except ValueError:
                continue
            assert (
                len(cs) == 1
            ), f"Only one connection set should be given from {glom_post_strat.name}."
            cs = self.scaffold.get_connectivity_set(cs[0])

            # We extract global ids to filter the postsynaptic neuron that match the ones from the dependency
            chunks = cs.pre_type.get_placement_set().get_all_chunks()
            # We filter glomeruli in chunks which are less than radius away from the current one.
            loc_glom_locs, loc_post_locs = (
                cs.load_connections()
                .from_(get_close_chunks(chunk, chunks, self.radius))
                .as_globals()
                .all()
            )

            # Filter unique glomeruli connections and grouping postsynaptic targets
            loc_glom_locs, ids = np.unique(loc_glom_locs[:, 0], return_inverse=True)
            glom_ids.extend(loc_glom_locs)
            # Sorting connections by presynaptic glomeruli id
            sorting = np.argsort(ids)
            # Returns the index positions of each glomerulus group in the sorted connection array
            _, ids = np.unique(ids[sorting], return_index=True)
            # Sorts the postsynaptic info and then splitting it according to the presynaptic ids
            post_locs.extend(np.split(loc_post_locs[sorting], ids[1:]))

            # Since we use global ids there is no need to filter chunks of the presynaptic cell placement set
            loc_glom_pos = np.concatenate(
                [cs.pre_type.get_placement_set().load_positions()[loc_glom_locs], loc_glom_pos]
            )

        return self._combine_unique_gloms(glom_ids, loc_glom_pos, post_locs)

    def _connect_type(self, pre_ps, post):
        # Chunks are sorted presynaptically so there should be only one chunk
        chunk = pre_ps.get_loaded_chunks()
        assert len(chunk) == 1, "There should be exactly one presynaptic chunk"
        chunk = chunk[0]

        # Extract all unique glomeruli connections from connection strategy dependencies grouped by glomeruli
        glom_pos = np.empty([0, 3])
        postsyn_locs = []
        ps_ids = []
        for glom_cell_type in self.glom_cell_types:
            ct_glom_ids = []
            ct_glom_pos = np.empty([0, 3])
            ct_postsyn_locs = []
            ct_ps_ids = []
            for i, post_ps in enumerate(post.placement):
                (loc_glom_ids, loc_glom_pos, loc_post_locs) = self._get_glom_cluster(
                    chunk, post_ps, glom_cell_type
                )
                ct_glom_ids.extend(loc_glom_ids)
                ct_glom_pos = np.concatenate([loc_glom_pos, ct_glom_pos])
                ct_postsyn_locs.extend(loc_post_locs)
                ct_ps_ids.extend([np.repeat(i, len(ids)) for ids in loc_post_locs])
                i += 1
            # group gloms from different postsynaptic cell type together (e.g. UBC, GrC)
            u_glom_ids, ct_glom_pos, ct_postsyn_locs = self._combine_unique_gloms(
                ct_glom_ids, ct_glom_pos, ct_postsyn_locs
            )
            glom_pos = np.concatenate([ct_glom_pos, glom_pos])
            postsyn_locs.extend(ct_postsyn_locs)
            ct_ps_ids = np.asarray(ct_ps_ids, dtype=object)
            for u_glom in u_glom_ids:
                ids = np.where(ct_glom_ids == u_glom)[0]
                ps_ids.append(np.concatenate(ct_ps_ids[ids].tolist(), dtype=int))

        golgi_pos = pre_ps.load_positions()
        # Cache morphologies and generate the morphologies iterator
        golgi_morphos = pre_ps.load_morphologies().iter_morphologies(cache=True, hard_cache=True)

        # TODO: implement random rounding and adapt tests.
        num_glom_to_connect = np.min([int(self.divergence), len(postsyn_locs)])
        n_conn = (
            len(golgi_pos)
            * num_glom_to_connect
            * np.max([len(post_conn) for post_conn in postsyn_locs])
        )
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        selected_ps = np.full(n_conn, -1, dtype=int)
        ptr = 0
        for i, golgi, morpho in zip(itertools.count(), golgi_pos, golgi_morphos):
            # Find terminal branches
            axon_branches = morpho.get_branches()
            terminal_branches_ids = np.nonzero([b.is_terminal for b in axon_branches])[0]
            axon_branches = np.take(axon_branches, terminal_branches_ids, axis=0)
            terminal_branches_ids = np.array([morpho.branches.index(b) for b in axon_branches])

            # Find the point-on-branch ids of the tips
            tips_coordinates = np.array([len(b.points) - 1 for b in axon_branches])

            # Compute and sort the distances between the golgi soma and the glomeruli
            to_connect = np.linalg.norm(golgi - glom_pos, axis=1)
            sorting = np.argsort(to_connect)
            to_connect = sorting[to_connect[sorting] <= self.radius]
            # Keep the closest glomeruli
            to_connect = to_connect[:num_glom_to_connect]
            # For each glomerulus, connect the corresponding postsyn cells directly to the current Golgi
            for post_conn in to_connect:
                take_post = np.array(postsyn_locs[post_conn])
                post_to_connect = len(take_post)
                # Select postsyn cells ids
                post_locs[ptr : ptr + post_to_connect] = take_post
                selected_ps[ptr : ptr + post_to_connect] = ps_ids[post_conn]
                # Select Golgi axon branch
                pre_locs[ptr : ptr + post_to_connect, 0] = i
                ids_branches = np.random.randint(
                    low=0, high=len(axon_branches), size=post_to_connect
                )
                pre_locs[ptr : ptr + post_to_connect, 1] = terminal_branches_ids[ids_branches]
                pre_locs[ptr : ptr + post_to_connect, 2] = tips_coordinates[ids_branches]
                ptr += post_to_connect

        for i, ps in enumerate(post.placement):
            selected = selected_ps[:ptr] == i
            # We use the scaffold placement set for each post cell type so that the global postsynaptic ids are used
            self.connect_cells(
                pre_ps,
                self.scaffold.get_placement_set(ps.cell_type),
                pre_locs[:ptr][selected],
                post_locs[:ptr][selected],
            )
