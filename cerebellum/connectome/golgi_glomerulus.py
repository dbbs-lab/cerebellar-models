"""
    Module for the configuration node of the Golgi to Glomerulus ConnectionStrategy
"""

import itertools

import numpy as np
from bsb import Chunk, ConfigurationError, ConnectionStrategy, config, refs


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

    def _pre_chunks_filter(self, chunk, cs):
        # Chunk distance from cs pre_types to provided chunk
        chunks = cs.pre_type.get_placement_set().get_all_chunks()
        # Fixme: Distance between chunk is done corner to corner. It might not detect all chunks #34
        pre_chunks = []
        for c in chunks:
            if np.linalg.norm(chunk * chunk.dimensions - c * c.dimensions) <= self.radius:
                pre_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return pre_chunks

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
            # We filter glomeruli in chunks which are less than radius away from the current one.
            loc_glom_locs, loc_post_locs = (
                cs.load_connections().from_(self._pre_chunks_filter(chunk, cs)).as_globals().all()
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

        # Extract unique glomerulus groups
        unique_gloms = np.unique(glom_ids)
        post_locs = np.asarray(post_locs, dtype=object)
        postsyn_connections = []
        postsyn_connections_branch_point = []
        glom_pos = np.zeros((len(unique_gloms), 3), dtype=float)
        for c, u_glom in enumerate(unique_gloms):
            ids = np.where(glom_ids == u_glom)[0]
            postsyn_connections.append(
                [loc_post_locs[0] for loc_post_locs in np.concatenate(post_locs[ids])]
            )
            postsyn_connections_branch_point.append(
                [loc_post_locs[1:] for loc_post_locs in np.concatenate(post_locs[ids])]
            )
            glom_pos[c] = loc_glom_pos[ids[0]]

        return (
            glom_pos,
            postsyn_connections,
            postsyn_connections_branch_point,
        )

    def _connect_type(self, pre_ps, post):
        # Chunks are sorted presynaptically so there should be only one chunk
        chunk = pre_ps.get_loaded_chunks()
        assert len(chunk) == 1, "There should be exactly one presynaptic chunk"
        chunk = chunk[0]

        # Extract all unique glomeruli connections from connection strategy dependencies grouped by glomeruli
        glom_pos = np.empty([0, 3])
        postsyn_connections = []
        postsyn_connections_branch_point = []
        post_ps_list = []
        ps_ids = []
        i = 0
        for post_ps in post.placement:
            for glom_cell_type in self.glom_cell_types:
                (
                    loc_glom_pos,
                    loc_postsyn_connections,
                    loc_postsyn_connections_branch_point,
                ) = self._get_glom_cluster(chunk, post_ps, glom_cell_type)
                glom_pos = np.concatenate([loc_glom_pos, glom_pos])
                postsyn_connections.extend(loc_postsyn_connections)
                postsyn_connections_branch_point.extend(loc_postsyn_connections_branch_point)
                ps_ids.extend(np.repeat(i, len(loc_postsyn_connections)))
                post_ps_list.append(self.scaffold.get_placement_set(post_ps.cell_type))
                i += 1

        golgi_pos = pre_ps.load_positions()
        # Cache morphologies and generate the morphologies iterator
        golgi_morphos = pre_ps.load_morphologies().iter_morphologies(cache=True, hard_cache=True)

        # TODO: implement random rounding and adapt tests.
        num_glom_to_connect = np.min([int(self.divergence), len(postsyn_connections)])
        n_conn = (
            len(golgi_pos)
            * num_glom_to_connect
            * np.max([len(post_conn) for post_conn in postsyn_connections])
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
            for j, post_conn in enumerate(to_connect):
                take_post = postsyn_connections[post_conn]
                post_to_connect = len(take_post)
                # Select postsyn cells ids
                post_locs[ptr : ptr + post_to_connect, 0] = take_post
                # Store branch-ids and points-on-branch-ids of the postsyn cells
                post_locs[ptr : ptr + post_to_connect, 1:] = postsyn_connections_branch_point[
                    post_conn
                ]
                selected_ps[ptr : ptr + post_to_connect] = ps_ids[post_conn]
                # Select Golgi axon branch
                pre_locs[ptr : ptr + post_to_connect, 0] = i
                ids_branches = np.random.randint(
                    low=0, high=len(axon_branches), size=post_to_connect
                )
                pre_locs[ptr : ptr + post_to_connect, 1] = terminal_branches_ids[ids_branches]
                pre_locs[ptr : ptr + post_to_connect, 2] = tips_coordinates[ids_branches]
                ptr += post_to_connect

        # So that the global postsynaptic ids are used
        for i, ps in enumerate(post_ps_list):
            selected = selected_ps[:ptr] == i
            self.connect_cells(
                pre_ps,
                ps,
                pre_locs[:ptr][selected],
                post_locs[:ptr][selected],
            )
