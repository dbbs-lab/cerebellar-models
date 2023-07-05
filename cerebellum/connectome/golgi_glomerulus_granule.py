import itertools
import numpy as np
from bsb.connectivity import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from scipy.stats.distributions import truncexpon
from bsb.morphologies import Morphology
from bsb.connectivity.strategy import Hemitype
from bsb.connectivity.strategy import HemitypeCollection


@config.node
class ConnectomeGolgiGlomerulusGranule(ConnectionStrategy):
    # Read vars from the configuration file
    # The radius is the maximum length of Golgi axon plus the radius of the granule
    # dendritic tree, needed in get_region_of_interest
    # to find the chunks containing the candidate cells to connect
    radius = config.attr(type=int, required=True)
    convergence = config.attr(type=int, required=True)
    intermediate = config.attr(type=Hemitype, required=True)
    
    #We override this method because it is easier to connect the Golgi cells
    #in a single post chunks to the granule cells in many pre-chunks
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, chunk)
        return pre, post

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        #print(type(ct))
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        # Look for chunks which are less than radius away from the current one in the xy plane
        # and for neighbouring chunks along z direction.
        # Note: All the chunks have the same dimensions
        for c in chunks:    
            dist = np.sqrt(
                np.power((chunk[0] - c[0]) * chunk.dimensions[0], 2)
                + np.power((chunk[1]  - c[1]) * chunk.dimensions[1], 2)
                + np.power((chunk[2]  - c[2]) * chunk.dimensions[2], 2)
            )
            if (dist < self.radius):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        golgi_pos = pre_ps.load_positions()
        granule_pos = post_ps.load_positions()

        print("Number of granule", len(post_ps))
        print("Number of Golgi", len(pre_ps))
        print("Number of gloms", len(self.intermediate.cell_types[0].get_placement_set(chunks=pre_ps.get_loaded_chunks())))

        #Look for connected glomeruli
        total_count = 0
        granule_connections = []
        granule_connections_branch_point = []
        glom_ids_list = []
        n_local_glom = 0

        #Sort the chunks in increasing chunk_id order
        pre_chunks = list(pre_ps.get_loaded_chunks())
        chunks_ids = [c.id for c in pre_chunks]
        sorted_ids = np.argsort(chunks_ids)
        sorted_chunks = []
        for i in sorted_ids:
            sorted_chunks.append(pre_chunks[i])

        #Get the glom_to_granule connections
        cs = self.scaffold.get_connectivity_set("glomerulus_to_granule")
        iter = cs.load_connections().from_(sorted_chunks)
        gid = iter.all()
        gid_local = gid[0]
        gid_global = gid[1]
        
        #Local are granule cells; global are the gloms.
        local_locs = gid_global[:,:]
        global_locs = gid_local[:,:]
        
        unique_gloms, unique_gloms_id = np.unique(global_locs,axis=0,return_index=True)
        unique_gloms_id = global_locs[unique_gloms_id,0]
        granule_connections = []
        granule_connections_branch_point = []

        for ugloms in unique_gloms:
            ids = np.where(global_locs[:,0] == ugloms[0])
            starting = np.min(local_locs[:,0])
            granule_connections.append(local_locs[ids][:,0]-starting)
            granule_connections_branch_point.append(local_locs[ids][:,1:2])
        
        #Consider only the glomeruli which are connected to at least a granule cell in the ROI.
        #Select only the gloms for which a connection is found
        glom_pos = self.intermediate.cell_types[0].get_placement_set(chunks=sorted_chunks).load_positions()[unique_gloms_id]
        n_glomeruli = len(glom_pos)
        n_golgi = len(golgi_pos)

        post_locs = np.full((n_glomeruli*n_golgi*len(post_ps)*40,3),-1,dtype=int)
        pre_locs = np.full((n_glomeruli*n_golgi*len(post_ps)*40,3),-1,dtype=int)
        
        #Consider only the glomeruli which are connected to at least a granule cell in the ROI.
        #Select only the gloms for which a connection is found

        #Cache morphologies and generate the morphologies iterator
        morpho_set = post_ps.load_morphologies()
        golgi_morphos = morpho_set.iter_morphologies(cache=True, hard_cache=True)

        #Cache Golgi morphologies 
        morphologies = []
        for morpho in pre_ct.get_morphologies():
            morphologies.append(morpho.load())

        ptr = 0

        #For each Golgi cell...
        for id_golgi, golgi, morpho in zip(itertools.count(), golgi_pos, golgi_morphos):

            axon_branches = morpho.get_branches()
            first_axon_branch_id = axon_branches.index(axon_branches[0])

            #Find terminal branches
            terminal_ids = np.full(len(axon_branches), 0, dtype=int)
            for i,b in enumerate(axon_branches):
                if b.is_terminal:
                    terminal_ids[i] = 1
            terminal_branches_ids = np.nonzero(terminal_ids)[0]

            #Keep only terminal branches
            axon_branches = np.take(axon_branches, terminal_branches_ids, axis=0)
            terminal_branches_ids = terminal_branches_ids + first_axon_branch_id

            #Find the point-on-branch ids of the tips
            tips_coordinates = np.full((len(axon_branches)), 0, dtype=float)
            for i,branch in enumerate(axon_branches):
                tips_coordinates[i] = len(branch.points[-1])-1

            #Compute the distances between the golgi baricenter and the glomeruli
            dist = np.sqrt(
                    np.power(golgi[0] - glom_pos[:, 0], 2)
                    + np.power(golgi[1] - glom_pos[:, 1], 2)
                    + np.power(golgi[2] - glom_pos[:, 2], 2)
                )
            
            #Select the 40 closer gloms 
            to_connect = np.argsort(dist)
            num_glom_to_connect = np.min([self.convergence,len(to_connect),len(granule_connections)])
            to_connect = to_connect[0:num_glom_to_connect]
            
            #For each glomerulus, connect the corresponding granule cells directly to the current Golgi
            for i in range(num_glom_to_connect):
                take_granule = granule_connections[to_connect[i]]
                #print(take_granule)
                #print(len(take_granule))
                #print(len(post_locs))
                granule_to_connect = len(take_granule)
                #Select cells ids
                post_locs[ptr:ptr+granule_to_connect,0] = take_granule
                pre_locs[ptr:ptr+granule_to_connect,0] = id_golgi
                #Store branch-ids and points-on-branch-ids of the granule cells
                take_branch_point = granule_connections_branch_point[to_connect[i]]
                post_locs[ptr:ptr+granule_to_connect,1:] = take_branch_point
                #Select Golgi axon branch
                ids_branches = np.random.randint(low=0,high=len(axon_branches),size=granule_to_connect)
                pre_locs[ptr:ptr+granule_to_connect,1] = terminal_branches_ids[ids_branches]
                pre_locs[ptr:ptr+granule_to_connect,2] = tips_coordinates[ids_branches]
                ptr += granule_to_connect
        
        print("Connected", n_golgi, "Golgi cells to", ptr, "granule cells.")
        #print(pre_ps.get_loaded_chunks())
        """
        for chunk in pre_ps.get_loaded_chunks():
            with pre_ps.chunk_context(chunk):
                ln = len(pre_ps)
                print("Golgi cells in:", chunk, ln)
        """
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
