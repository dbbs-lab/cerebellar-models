from asyncio import unix_events
import chunk
import itertools
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from bsb.morphologies import Morphology
from bsb.connectivity.strategy import Hemitype
from bsb.connectivity.strategy import HemitypeCollection

class TooFewGlomeruliClusters(Exception):
    pass

@config.node
class ConnectomeGlomerulusGranule(ConnectionStrategy):
    max_radius = config.attr(type=int, required=True)
    convergence = config.attr(type=int, required=True)
    #We need acces to mossy fibers to know the glomeruli clusters
    prepresynaptic = config.attr(type=Hemitype, required=True)

    #We need to override the default _get_connect_args_from_job function because 
    #we need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    #In this way, if for each glomerulus in the ROI we select randomly the mossy fiber
    #to which is connected, the divergence ratio is automatically ~20.
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, chunk)
        return pre, post

    def get_region_of_interest(self, chunk):
        selected_chunks = []
        
        #We look for chunks containing the granule cells to connect
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] * chunk.dimensions[0] - c[0] * c.dimensions[0], 2)
                + np.power(chunk[1] * chunk.dimensions[1] - c[1] * c.dimensions[0], 2)
                + np.power(chunk[2] * chunk.dimensions[2] - c[2] * c.dimensions[0], 2)
            )
            #If c is intersecting the granular layer partition and if it is close enough to chunk, then it is in the ROI
            if (dist <= self.max_radius) and c[1]*c.dimensions[1] < self.scaffold.partitions.granular_layer.thickness:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        return selected_chunks

    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        glom_pos = pre_ps.load_positions()
        gran_pos = post_ps.load_positions()
        n_granular = len(gran_pos)
        max_connections = self.convergence
        n_conn = n_granular * max_connections
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        #Find the glomeruli clusters
        cs = self.scaffold.get_connectivity_set("mossy_fibers_to_glomerulus")
        iter = cs.load_connections().to(pre_ps.get_loaded_chunks()).as_globals()
        clusters = []
        
        gid_global,gid_local = iter.all()
        unique_mossy = np.unique(gid_global,axis=0)

        for current in unique_mossy:
            glom_ids = np.where((gid_global[:,0]==current[0]))
            starting = np.min(gid_local[:,0])
            clusters.append(gid_local[glom_ids[0],0]-starting)

        num_clusters = len(clusters)

        if (num_clusters < 4):
            raise TooFewGlomeruliClusters("Less then 4 clusters of glomeruli have been found. Check the densities of mossy fibers and glomeruli in the configuration file.")
        
        #Cache morphologies 
        morpho_set = post_ps.load_morphologies()

        #We keep track of the entries of pre_locs and post_locs we actually used.
        ptr = 0 
        gran_morphos = morpho_set.iter_morphologies(cache=True, hard_cache=True)
        for i, grpos, morpho in zip(itertools.count(), gran_pos, gran_morphos):
            
            dendrites = morpho.get_branches()

            #Get the starting branch id of the dendritic branches
            first_dendride_id = morpho.branches.index(dendrites[0])

            gr_connections = 0

            #Select (randomly) the order of the clusters
            cluster_idx = np.arange(0,num_clusters)
            np.random.shuffle(cluster_idx)

            #Select (randomly) the order of the dendrites
            #We generate it now because it is faster than
            #generating a random id in the cluster loop
            dendrites_idx = np.arange(0,len(dendrites))

            np.random.shuffle(dendrites_idx)  

            #Try to select a cell from 4 clusters satisfying the conditions
            for nc in cluster_idx:
                dist = np.sqrt(
                    np.power(grpos[0] - glom_pos[clusters[nc]][:, 0], 2)
                    + np.power(grpos[1] - glom_pos[clusters[nc]][:, 1], 2)
                    + np.power(grpos[2] - glom_pos[clusters[nc]][:, 2], 2)
                )
                bool_arr = dist < self.max_radius
                sorted_indices = np.nonzero(bool_arr)[0]
                if (len(sorted_indices)>0 and gr_connections<4):
                    rnd = np.random.randint(low=0,high=len(sorted_indices))
                    #Id of the granule cell
                    post_locs[ptr + gr_connections, 0] = i
                    #Id of the glomerulus, randomly selected between the avaiable ones
                    pre_locs[ptr + gr_connections, 0] = clusters[nc][sorted_indices[rnd]]
                    #Select one of the 4 dendrites
                    dendrite = dendrites[dendrites_idx[gr_connections]]
                    #Select the terminal point of the branch
                    tip = len(dendrite)-1
                    post_locs[ptr + gr_connections, 1] = first_dendride_id+dendrites_idx[gr_connections]
                    post_locs[ptr + gr_connections, 2] = tip
                    gr_connections += 1
                    
                #When 4 connection are formed, exit the loop
                if (gr_connections >= 4):
                    break
            
            #If there are some free dendrites, connect them to the closest glomeruli,
            #even if they do not satisfy the geometric condtitions.
            if (gr_connections<self.convergence):
                #Connect the clostest glomeruli from four random clusters
                for nc in cluster_idx:
                    dist = np.sqrt(
                        np.power(grpos[0] - glom_pos[clusters[nc]][:, 0], 2)
                        + np.power(grpos[1] - glom_pos[clusters[nc]][:, 1], 2)
                        + np.power(grpos[2] - glom_pos[clusters[nc]][:, 2], 2)
                    )
                    min_dist_idx = np.argmin(dist)
                    #Check if the chosen glom is already connected to the granule
                    if clusters[nc][min_dist_idx] not in pre_locs[ptr : ptr + gr_connections, 0]:

                        post_locs[ptr + gr_connections, 0] = i
                        #Id of the glomerulus, randomly selected between the avaiable ones
                        pre_locs[ptr + gr_connections, 0] = clusters[nc][min_dist_idx]
                        #Select one of the 4 dendrites
                        dendrite = dendrites[dendrites_idx[gr_connections]]
                        #Select the terminal point of the branch
                        tip = len(dendrite)-1
                        post_locs[ptr + gr_connections, 1] = first_dendride_id+dendrites_idx[gr_connections]
                        post_locs[ptr + gr_connections, 2] = tip
                        gr_connections += 1
                        #When 4 connection are formed, exit the loop
                        if (gr_connections >= 4):
                            break
                
            ptr += gr_connections     

        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
