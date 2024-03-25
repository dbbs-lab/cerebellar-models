import itertools
import numpy as np
from bsb import config, refs, ConnectionStrategy, ConfigurationError, Chunk


@config.node
class ConnectomeGolgiGlomerulusGranule(ConnectionStrategy):
    # Read vars from the configuration file
    # The radius is the maximum length of Golgi axon plus the radius of the granule
    # dendritic tree, needed in get_region_of_interest
    # to find the chunks containing the candidate cells to connect
    convergence = config.attr(type=int, required=True)
    glom_grc_tag = config.attr(required=True)
    radius = config.attr(type=int, required=True)

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _get_glom_cluster(self, pre_ps, post_ps):
        # Get the glom_to_granule connections
        cs = self.scaffold.get_connectivity_set(self.glom_grc_tag)
        # Chunks are sorted pre-synaptically so there should be only one chunk
        chunk = pre_ps.get_loaded_chunks()
        assert len(chunk) == 1, "There should be exactly one presynaptic chunk"
        chunk = chunk[0]
        # Distance from glom to golgis
        chunks = cs.pre_type.get_placement_set().get_all_chunks()
        # Look for chunks containing glom which are less than radius away from the current one.
        # Fixme: Distance between chunk is done corner to corner. It might not detect all chunks
        pre_chunks = []
        for c in chunks:
            if (
                np.linalg.norm(chunk * chunk.dimensions - c * c.dimensions)
                <= self.radius
            ):
                pre_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        # We need global ids to filter grc that match the ones from the dependency
        iter = cs.load_connections().from_(pre_chunks).as_globals()
        _, grc_locs = iter.all()
        iter = cs.load_connections().from_(pre_chunks)
        glom_locs, _ = iter.all()
        unique_gloms = np.unique(glom_locs[:, 0])
        granule_connections = []
        granule_connections_branch_point = []
        for u_glom in unique_gloms:
            ids = np.where(glom_locs[:, 0] == u_glom)[0]
            granule_connections.append(grc_locs[ids, 0])
            post_ps.load_ids()
            granule_connections_branch_point.append(grc_locs[ids, 1:])

        glom_pos = cs.pre_type.get_placement_set(chunks=pre_chunks).load_positions()[
            unique_gloms
        ]

        return (
            unique_gloms,
            glom_pos,
            granule_connections,
            granule_connections_branch_point,
        )

    def _connect_type(self, pre_ps, post_ps):

        # Consider only the glomeruli which are connected to at least a granule cell in the ROI.
        # Select only the gloms for which a connection is found
        golgi_pos = pre_ps.load_positions()
        (
            unique_gloms,
            glom_pos,
            granule_connections,
            granule_connections_branch_point,
        ) = self._get_glom_cluster(pre_ps, post_ps)

        # Cache morphologies and generate the morphologies iterator
        golgi_morphos = pre_ps.load_morphologies().iter_morphologies(
            cache=True, hard_cache=True
        )

        num_glom_to_connect = np.min([self.convergence, len(granule_connections)])
        n_conn = (
            len(golgi_pos)
            * num_glom_to_connect
            * np.max([len(gr_conn) for gr_conn in granule_connections])
        )
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        ptr = 0
        # For each Golgi cell...
        for i, golgi, morpho in zip(itertools.count(), golgi_pos, golgi_morphos):
            # Find terminal branches
            axon_branches = morpho.get_branches()
            terminal_branches_ids = np.nonzero([b.is_terminal for b in axon_branches])[
                0
            ]
            axon_branches = np.take(axon_branches, terminal_branches_ids, axis=0)
            terminal_branches_ids = np.array(
                [morpho.branches.index(b) for b in axon_branches]
            )

            # Find the point-on-branch ids of the tips
            tips_coordinates = np.array([len(b.points) - 1 for b in axon_branches])

            # Compute and sort the distances between the golgi soma and the glomeruli
            to_connect = np.argsort(np.linalg.norm(golgi - glom_pos, axis=1))
            # Keep the closest glomeruli
            to_connect = to_connect[:num_glom_to_connect]
            # For each glomerulus, connect the corresponding granule cells directly to the current
            # Golgi
            for j in range(num_glom_to_connect):
                take_granule = granule_connections[to_connect[j]]
                granule_to_connect = len(take_granule)
                # Select granule cells ids
                post_locs[ptr : ptr + granule_to_connect, 0] = take_granule
                # Store branch-ids and points-on-branch-ids of the granule cells
                post_locs[ptr : ptr + granule_to_connect, 1:] = (
                    granule_connections_branch_point[to_connect[j]]
                )
                # Select Golgi axon branch
                pre_locs[ptr : ptr + granule_to_connect, 0] = i
                ids_branches = np.random.randint(
                    low=0, high=len(axon_branches), size=granule_to_connect
                )
                pre_locs[ptr : ptr + granule_to_connect, 1] = terminal_branches_ids[
                    ids_branches
                ]
                pre_locs[ptr : ptr + granule_to_connect, 2] = tips_coordinates[
                    ids_branches
                ]
                ptr += granule_to_connect

        # So that the global ids are used
        self.connect_cells(
            pre_ps,
            self.scaffold.get_placement_set(post_ps.cell_type),
            pre_locs[:ptr],
            post_locs[:ptr],
        )
