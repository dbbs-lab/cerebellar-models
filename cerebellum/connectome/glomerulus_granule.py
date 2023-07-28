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
    max_convergence = config.attr(type=int, required=True)
    #We need acces to mossy fibers to know the glomeruli clusters
    prepresynaptic = config.attr(type=Hemitype, required=True)

    #We need to override the default _get_connect_args_from_job function because 
    #we need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    #In this way, given  if for each glomerulus in the ROI we select randomly the mossy fiber
    #to connect, the divergence ratio is automatically around 20.
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
            if (dist <= self.max_radius) and c[1]*c.dimensions[0] < self.scaffold.partitions.granular_layer.thickness:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        return selected_chunks

    def connect(self, pre, post):
        for post_ct, post_ps in post.placement.items():
            pre_ps_mf = None
            pre_ps_ubc = None
            for pre_ct, pre_ps in pre.placement.items():        
                if (pre_ct.name == "glomerulus"):
                    pre_ps_mf = pre_ps
                if (pre_ct.name == "ubc_glomerulus"):
                    pre_ps_ubc = pre_ps
            self._connect_type(pre_ps_mf, pre_ps_ubc, post_ps)

    def _connect_type(self, pre_ps_mf, pre_ps_ubc, post_ps):
        
        glom_mf_pos = pre_ps_mf.load_positions()
        glom_ubc_pos = pre_ps_ubc.load_positions()
        gran_pos = post_ps.load_positions()

        print("Pre",pre_ps_mf.get_loaded_chunks())
        print("Post",post_ps.get_loaded_chunks())
        
        print("Connecting", len(gran_pos), "granule cells in chunk: ",post_ps.get_loaded_chunks())
        print("To", len(glom_mf_pos), "glomeruli in chunk: ",pre_ps_mf.get_loaded_chunks())
        print("Number of glomeruli in ROI:",len(glom_mf_pos))
        print("Number of granule in ROI:",len(gran_pos))
        
        n_glom_mf = len(glom_mf_pos)
        n_glom_ubc = len(glom_ubc_pos)
        n_gran = len(gran_pos)
        max_connections = self.max_convergence
        n_conn = n_gran * max_connections
        
        #For the sake of speed we save the connectivity data separately for glom_mf and glom_ubc
        
        pre_locs_mf = np.full((n_conn, 3), -1, dtype=int)
        post_locs_mf = np.full((n_conn, 3), -1, dtype=int)
        pre_locs_ubc = np.full((n_conn, 3), -1, dtype=int)
        post_locs_ubc = np.full((n_conn, 3), -1, dtype=int)
        
        #Find the glomeruli clusters
        cs = self.scaffold.get_connectivity_set("mossy_fibers_to_glomerulus")
        iter = cs.load_connections().to(pre_ps_mf.get_loaded_chunks()).as_globals()
        clusters = []
        
        gid_global,gid_local = iter.all()
        unique_mossy = np.unique(gid_global,axis=0)

        for current in unique_mossy:
            glom_ids = np.where((gid_global[:,0]==current[0]))
            starting = np.min(gid_local[:,0])
            clusters.append(gid_local[glom_ids[0],0]-starting)

        num_clusters = len(clusters)
        print("Number of clusters:", num_clusters)
        for i,cl in enumerate(clusters):
            print("Glomeruli in cluster", i, ":", len(cl))
        
        if (num_clusters < 4):
            raise TooFewGlomeruliClusters("Less then 4 clusters of glomeruli have been found. Check the densities of mossy fibers and glomeruli in the configuration file.")
        
        #Cache morphologies 
        morpho_set = post_ps.load_morphologies()

        #We keep track of the entries of pre_locs and post_locs we actually used.
        ptr_mf = 0 
        ptr_ubc = 0

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
            
            selected_cluster = []
            selected_indices = []
            selected_distances = []

            for nc in cluster_idx:
                dist_mf = np.sqrt(
                    np.power(grpos[0] - glom_mf_pos[clusters[nc]][:, 0], 2)
                    + np.power(grpos[1] - glom_mf_pos[clusters[nc]][:, 1], 2)
                    + np.power(grpos[2] - glom_mf_pos[clusters[nc]][:, 2], 2)
                )
                
                sorted_by_distance = np.sort(dist_mf)
                ids_sorted_by_distance = np.argsort(dist_mf)
                bool_arr_mf = sorted_by_distance < self.max_radius
                sorted_by_distance = sorted_by_distance[bool_arr_mf]
                ids_sorted_by_distance = ids_sorted_by_distance[bool_arr_mf]
                ln = np.count_nonzero(bool_arr_mf)
                if ln > 0 and ln < 5:
                    for id in range(0,ln):
                        selected_indices.append(ids_sorted_by_distance[id])
                        selected_distances.append(sorted_by_distance[id])
                        selected_cluster.append(nc)
            
            #Select first four gloms (according to distance) in each cluster

            if (len(glom_ubc_pos)) > 0:
                dist_ubc = np.linalg.norm(grpos-glom_ubc_pos,axis=1)
                sorted_by_distance = np.sort(dist_ubc)
                ids_sorted_by_distance = np.argsort(dist_ubc)
                mask_ubc = sorted_by_distance < self.max_radius
                sorted_by_distance = sorted_by_distance[mask_ubc]
                ids_sorted_by_distance = ids_sorted_by_distance[mask_ubc]
                ln = np.count_nonzero(mask_ubc)
                if ln > 0:
                    for idx in range(0,min(4,ln)):
                        selected_indices.append(ids_sorted_by_distance[idx])
                        selected_distances.append(sorted_by_distance[idx])
                        selected_cluster.append(-1)

            selected_indices = np.array(selected_indices)
            selected_distances = np.array(selected_distances)
            selected_cluster = np.array(selected_cluster)

            global_sorted_idx = np.argsort(selected_distances)
            selected_cluster = selected_cluster[global_sorted_idx]
            selected_indices = selected_indices[global_sorted_idx]

            if len(selected_indices) >= 4:
                for idx in range(min(4,len(selected_indices))):                  
                    
                    #If the id of the cluster is -1, we are connecting a ubc_glom
                    if (selected_cluster[idx] == -1 and gr_connections < 4):
                    #Select one of the 4 dendrites
                        dendrite = dendrites[dendrites_idx[gr_connections]]
                        #Select the terminal point of the branch
                        tip = len(dendrite)-1
                        post_locs_ubc[ptr_ubc, 0] = i
                        post_locs_ubc[ptr_ubc, 1] = first_dendride_id+dendrites_idx[gr_connections]
                        post_locs_ubc[ptr_ubc, 2] = tip
                        pre_locs_ubc[ptr_ubc, 0] = selected_indices[idx]
                        gr_connections += 1
                        ptr_ubc = ptr_ubc + 1

                    #If the id of the cluster is not -1, we are connecting a mf_glom
                    if (selected_cluster[idx] != -1 and gr_connections < 4):
                        
                        cluster_id = selected_cluster[idx]
                        #Id of the granule cell
                        post_locs_mf[ptr_mf, 0] = i
                        #Id of the glomerulus, randomly selected between the avaiable ones
                        pre_locs_mf[ptr_mf, 0] = clusters[cluster_id][selected_indices[idx]]
                        #Select one of the 4 dendrites
                        dendrite = dendrites[dendrites_idx[gr_connections]]
                        #Select the terminal point of the branch
                        tip = len(dendrite)-1
                        post_locs_mf[ptr_mf, 1] = first_dendride_id+dendrites_idx[gr_connections]
                        post_locs_mf[ptr_mf, 2] = tip
                        gr_connections += 1
                        ptr_mf = ptr_mf + 1 
                        
            #If there are some free dendrites, connect them to the closest glomeruli,
            #even if they do not satisfy the geometric condtitions.
            if (gr_connections<self.max_convergence):
                
                #Connect the clostest glomeruli from four random clusters
                for nc in cluster_idx:
                    dist = np.sqrt(
                        np.power(grpos[0] - glom_mf_pos[clusters[nc]][:, 0], 2)
                        + np.power(grpos[1] - glom_mf_pos[clusters[nc]][:, 1], 2)
                        + np.power(grpos[2] - glom_mf_pos[clusters[nc]][:, 2], 2)
                    )
                    min_dist_idx = np.argmin(dist)
                    #Check if the chosen glom is already connected to the granule
                    if clusters[nc][min_dist_idx] not in pre_locs_mf[ptr_mf : ptr_mf + gr_connections, 0]:

                        post_locs_mf[ptr_mf + gr_connections, 0] = i
                        #Id of the glomerulus, randomly selected between the avaiable ones
                        pre_locs_mf[ptr_mf + gr_connections, 0] = clusters[nc][min_dist_idx]
                        #Select one of the 4 dendrites
                        dendrite = dendrites[dendrites_idx[gr_connections]]
                        #Select the terminal point of the branch
                        tip = len(dendrite)-1
                        post_locs_mf[ptr_mf + gr_connections, 1] = first_dendride_id+dendrites_idx[gr_connections]
                        post_locs_mf[ptr_mf + gr_connections, 2] = tip
                        gr_connections += 1
                        #When 4 connection are formed, exit the loop
                        if (gr_connections >= 4):
                            break
                
                ptr_mf += gr_connections 

        #print("Connected", n_gran, "granular cells.\n")
        #print(pre_locs[:ptr], post_locs[:ptr])
        
        if ptr_mf > 0:
            self.connect_cells(pre_ps_mf, post_ps, pre_locs_mf[:ptr_mf], post_locs_mf[:ptr_mf])
        if ptr_ubc:
            self.connect_cells(pre_ps_ubc, post_ps, pre_locs_ubc[:ptr_ubc], post_locs_ubc[:ptr_ubc])


