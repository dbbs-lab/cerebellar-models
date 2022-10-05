import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from bsb.morphologies import Morphology
from bsb.connectivity.strategy import Hemitype
from bsb.connectivity.strategy import HemitypeCollection


@config.node
class ConnectomeGlomerulusGranule(ConnectionStrategy):
    x_length = config.attr(type=int, required=True)
    z_length = config.attr(type=int, required=True)
    max_radius = config.attr(type=int, required=True)
    convergence = config.attr(type=int, required=True)
    # We need acces to mossy fibers to know the glomeruli clusters
    prepresynaptic = config.attr(type=Hemitype, required=True)

    """def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, [chunk])
        return pre, post"""

    def get_region_of_interest(self, chunk):
        selected_chunks = []

        # We look for chunks containing the granule cells to connect
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()

        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] * chunk.dimensions[0] - c[0] * c.dimensions[0], 2)
                + np.power(chunk[1] * chunk.dimensions[1] - c[1] * c.dimensions[0], 2)
                + np.power(chunk[2] * chunk.dimensions[2] - c[2] * c.dimensions[0], 2)
            )
            if dist < self.max_radius:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        # We look for chunks containing the mossy fibers
        ct = self.prepresynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()

        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[2] - c[2])
            x_dist = x_dist * chunk.dimensions[0]
            z_dist = z_dist * chunk.dimensions[2]

            if x_dist < self.x_length / 2 and z_dist < self.z_length / 2:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        return selected_chunks

    def connect(self, pre, post):
        # pre_type = pre.cell_types[0]
        # post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):

        glom_pos = pre_ps.load_positions()
        gran_pos = post_ps.load_positions()

        print("Pre", pre_ps.get_loaded_chunks())
        print("Post", post_ps.get_loaded_chunks())

        print(
            "Connecting",
            len(gran_pos),
            "granule cells in chunk: ",
            post_ps.get_loaded_chunks(),
        )
        print("Number of glomeruli in ROI:", len(glom_pos))
        n_glom = len(glom_pos)
        n_gran = len(gran_pos)
        max_connections = self.convergence
        n_conn = n_gran * max_connections
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        # Find the glomeruli clusters
        cs = self.scaffold.get_connectivity_set("mossy_fibers_to_glomerulus")
        clusters = []
        total_count = 0
        """for lchunk, itr in cs.nested_iter_connections("inc",local_=pre_ps.get_loaded_chunks(),global_=pre_ps.get_loaded_chunks()):
            for gchunk, (global_locs,local_locs) in itr:
                #local_n_glom = pre_ct.get_placement_set(chunks=[lchunk])
                unique_mossy = np.unique(local_locs,axis=0)
                print(unique_mossy)
                for i,current in enumerate(unique_mossy):
                    glom_ids = np.where(local_locs[:,0]==current[0])
                    #Update the global id
                    temp = local_locs[:,0] + total_count 
                    clusters.append(temp[glom_ids[0]])
            #Update the total count to get the correct indices
            total_count += len(pre_ct.get_placement_set(chunks=[lchunk]))"""

        for lchunk, itr in cs.nested_iter_connections(
            "out",
            local_=post_ps.get_loaded_chunks(),
            global_=pre_ps.get_loaded_chunks(),
        ):
            print("Local chunk:", lchunk)
            for gchunk, (local_locs, global_locs) in itr:
                # print(gchunk)
                print(local_locs)
                # local_n_glom = pre_ct.get_placement_set(chunks=[lchunk])
                unique_mossy = np.unique(local_locs, axis=0)
                # print(unique_mossy)
                for i, current in enumerate(unique_mossy):
                    glom_ids = np.where(local_locs[:, 0] == current[0])
                    # Update the global id
                    print(global_locs)
                    temp = np.unique(global_locs[:, 0], axis=0) + total_count
                    print(temp)
                    clusters.append(temp[glom_ids[0]])
            # Update the total count to get the correct indices
            # total_count += len(post_ct.get_placement_set(chunks=[lchunk]))
            print(total_count)
        print(total_count, "vs", n_glom)
        num_clusters = len(clusters)
        print("Number of clusters:", num_clusters)
        for i, cl in enumerate(clusters):
            print("Glomeruli in cluster", i, ":", len(cl))

        print(clusters)

        # Cache morphologies
        morphologies = []
        for morpho in post_ct.get_morphologies():
            morphologies.append(morpho.load())
        # Get the granule-morphology correspondence
        morpho_id = post_ps.load_morphologies().get_indices()

        # Get the name of the branches to connect
        branches_names = []
        for string in self.postsynaptic.subcell_labels:
            branches_names.append(string)

        # We keep track of the entries of pre_locs and post_locs we actually used.
        ptr = 0

        for i, grpos in enumerate(gran_pos):

            # Get the morphology of the current granule cell
            morpho = morphologies[morpho_id[i]]
            dendrites = morpho.get_branches(labels=branches_names)

            # Get the starting branch id of the dendritic branches
            first_dendride_id = morpho.branches.index(dendrites[0])

            gr_connections = 0

            # Select (randomly) the order of the clusters
            cluster_idx = np.arange(0, num_clusters)
            np.random.shuffle(cluster_idx)

            # Select (randomly) the order of the dendrites
            # We generate it now because it is faster than
            # generating a random id in the cluster loop
            dendrites_idx = np.arange(0, len(dendrites))

            np.random.shuffle(dendrites_idx)

            # Try to select a cell from 4 clusters satisfying the conditions
            for nc in cluster_idx:
                dist = np.sqrt(
                    np.power(grpos[0] - glom_pos[clusters[nc]][:, 0], 2)
                    + np.power(grpos[1] - glom_pos[clusters[nc]][:, 1], 2)
                    + np.power(grpos[2] - glom_pos[clusters[nc]][:, 2], 2)
                )
                bool_arr = dist < self.max_radius
                sorted_indices = np.nonzero(bool_arr)[0]
                if len(sorted_indices) > 0 and gr_connections < 4:
                    rnd = np.random.randint(low=0, high=len(sorted_indices))
                    # Id of the granule cell
                    post_locs[ptr + gr_connections, 0] = i
                    # Id of the glomerulus, randomly selected between the avaiable ones
                    pre_locs[ptr + gr_connections, 0] = clusters[nc][
                        sorted_indices[rnd]
                    ]

                    # Select one of the 4 dendrites
                    dendrite = dendrites[dendrites_idx[gr_connections]]
                    # Select the terminal point of the branch
                    tip = len(dendrite) - 1
                    post_locs[ptr + gr_connections, 1] = (
                        first_dendride_id + dendrites_idx[gr_connections]
                    )
                    post_locs[ptr + gr_connections, 2] = tip
                    gr_connections += 1

                # When 4 connection are formed, exit the loop
                if gr_connections >= 4:
                    break

            # If there are some free dendrites, connect them to the closest glomeruli,
            # even if they do not satisfy the geometric condtitions.
            if gr_connections < self.convergence:
                # Connect the clostest glomeruli from four random clusters
                for nc in cluster_idx:
                    dist = np.sqrt(
                        np.power(grpos[0] - glom_pos[clusters[nc]][:, 0], 2)
                        + np.power(grpos[1] - glom_pos[clusters[nc]][:, 1], 2)
                        + np.power(grpos[2] - glom_pos[clusters[nc]][:, 2], 2)
                    )
                    min_dist_idx = np.argmin(dist)
                    post_locs[ptr + gr_connections, 0] = i
                    # Id of the glomerulus, randomly selected between the avaiable ones
                    pre_locs[ptr + gr_connections, 0] = clusters[nc][min_dist_idx]
                    # Select one of the 4 dendrites
                    dendrite = dendrites[dendrites_idx[gr_connections]]
                    # Select the terminal point of the branch
                    tip = len(dendrite) - 1
                    post_locs[ptr + gr_connections, 1] = (
                        first_dendride_id + dendrites_idx[gr_connections]
                    )
                    post_locs[ptr + gr_connections, 2] = tip
                    gr_connections += 1
                    # When 4 connection are formed, exit the loop
                    if gr_connections >= 4:
                        break

            ptr += gr_connections

        print("Connected", n_gran, "granular cells.\n")
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
