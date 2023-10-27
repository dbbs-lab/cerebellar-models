from types import CellType
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from bsb.cell_types import CellType
from bsb.connectivity.strategy import HemitypeNode

@config.node
class ConnectomeMossyGlomerulus(ConnectionStrategy):
    x_length = config.attr(type=int, required=True)
    z_length = config.attr(type=int, required=True)
    x_sigma = config.attr(type=int, required=True)
    z_sigma = config.attr(type=int, required=True)
    mean = config.attr(type=int, required=True)
    sigma = config.attr(type=int, required=True)
    glom_granule_max_radius = config.attr(type=int, required=True)
    convergence = config.attr(type=int, required=True)
    intermediate = config.attr(type=HemitypeNode, required=True)

    def get_region_of_interest(self, chunk):
        #TODO rewrite this function to employ numpy vectorized operations

        #Find the chunks for the connection mossy fibers - glomeruli
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[2]- c[2]) 
            x_dist = x_dist * chunk.dimensions[0]
            z_dist = z_dist * chunk.dimensions[2]

            if (x_dist < self.x_length/2 and  z_dist < self.z_length/2):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        


        """
        #Find the chunks for the connection glomeruli - granule cells
        ct = self.intermediate.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] * chunk.dimensions[0] - c[0] * c.dimensions[0], 2)
                + np.power(chunk[1] * chunk.dimensions[1] - c[1] * c.dimensions[0], 2)
                + np.power(chunk[2] * chunk.dimensions[2] - c[2] * c.dimensions[0], 2)
            )
            print(dist)

            if (dist < self.glom_granule_max_radius):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        """
        selected_chunks = [chunk]
        return selected_chunks

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps.cell_type, pre_ps, post_ps.cell_type, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        print("##########################################################")
        print(pre_ps.get_loaded_chunks())
        print(post_ps.get_loaded_chunks())
        print("Number of MF:")
        print(len(pre_ps))
        print("Number of gloms:")
        print(len(post_ps))
        
        print("GRANULE")
        print(len(self.intermediate.cell_types[0].get_placement_set()))
        print(len(self.intermediate.cell_types[0].get_morphologies()))
        granule_pos = self.intermediate.cell_types[0].get_placement_set().load_positions()
        print(self.intermediate.cell_types[0].get_placement_set().get_loaded_chunks())
        mossy_pos = pre_ps.load_positions()
        glomeruli_pos = post_ps.load_positions()
        n_mossy = len(mossy_pos)
        n_glom = len(glomeruli_pos)
        print("Numero mossy:", n_mossy)
        print("Numero glomeruli:", n_glom)
        n_conn = n_glom * n_mossy
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        # Since a glomerulus can be connected only to one mossy fiber, we need
        # to keep track of the glomeruli which have already been connected.
        print("Initialize free gloms ids")
        free_glomeruli = np.full(n_glom, True, dtype=bool)
        print("FREE GLOMS")
        print(free_glomeruli)
        #mossy_pos = np.random.rand(0,n_mossy, len(id_non_connected[0]), dtype=int)

        # Find glomeruli to connect         
        #    
        #To keep track of the clusters we use a list of numpy arrays containing 
        #the ids of the glomeruli associated to the i-th mossy fiber.
        cluster = []
        for i, mossy in enumerate(mossy_pos):
            cluster.append([])
        ptr = 0
        print(pre_ps.get_loaded_chunks())
        for i, mossy in enumerate(mossy_pos):
            print("*******************************")
            print("MOSSY NUMBER:", i)
            
            #Shuffle free_glomeruli and glomeruli_pos the same way
            #p = np.random.permutation(len(glomeruli_pos))
            #free_glomeruli = free_glomeruli[p]
            #s_glomeruli_pos = glomeruli_pos[p,:]

            x_dist = np.fabs(mossy[0] - glomeruli_pos[:,0])
            z_dist = np.fabs(mossy[2]- glomeruli_pos[:,2])

            #Select the candidates glomeruli using a distance-based probability rule
            roll_x = np.random.normal(self.x_length/2, self.x_sigma, n_glom)
            roll_z = np.random.normal(self.z_length/2, self.z_sigma, n_glom)
            candidates_x =  x_dist < roll_x
            #candidates_x = np.nonzero(candidates_x)
            candidates_z =  z_dist < roll_z
            #candidates_z = np.nonzero(candidates_z)


            final_candidates = (candidates_x & candidates_z) & free_glomeruli
            final_candidates = np.nonzero(final_candidates)

            
            
            print(final_candidates)
            
            #Draw the number of glomeruli to connect to the current mossy fiber
            # following a Gaussian distribution
            glomeruli_to_connect = int(np.random.normal(self.mean, self.sigma))
            print("N_TO_CONNECT:",glomeruli_to_connect)
            #If the good glomeruli are less than glomeruli_to_connect, connect them all;
            # Otherwise, select the first good glomeruli_to_connect glomeruli
            num_candidates = len(final_candidates[0])
            print("Num candidates:",num_candidates)
            if (num_candidates <= glomeruli_to_connect):
                pre_locs[ptr:ptr+num_candidates,0] = i
                print("CAND")
                print(final_candidates)
                post_locs[ptr:ptr+num_candidates,0] = final_candidates[0]
                #free_glomeruli = free_glomeruli[num_candidates:]
                free_glomeruli[final_candidates[0]] = False
                print(num_candidates)
                ptr += num_candidates
                cluster[i] = final_candidates[0][0:num_candidates]
            else:
                pre_locs[ptr:ptr+glomeruli_to_connect,0] = i
                print(final_candidates[0][0:glomeruli_to_connect])
                post_locs[ptr:ptr+glomeruli_to_connect,0] = final_candidates[0][0:glomeruli_to_connect]
                #free_glomeruli = free_glomeruli[glomeruli_to_connect:]
                free_glomeruli[final_candidates[0][:glomeruli_to_connect]] = False
                print(glomeruli_to_connect)
                ptr += glomeruli_to_connect
                cluster[i] = final_candidates[0][0:glomeruli_to_connect]
            #Add to a cluster
        print(free_glomeruli)
        
        #If there are some glomeruli which are not connected, connect them to a random near mossy fiber
        id_non_connected = np.nonzero(free_glomeruli)[0]
        print("NON CONNESSI:",id_non_connected)
        print("LUNGHEZZA NON CONNESSI:", len(id_non_connected))
        #Draw some random numbers to be used as the indices of mossy fibers to connect to
        # each non connected glomerulus.
        random_ids = np.random.randint(low=0,high=n_mossy, size=len(id_non_connected), dtype=int)
        print("random_ids")
        #print(np.random.randint(0,n_mossy, len(id_non_connected)))
        print(random_ids)
        if (ptr>0):
            for j,id_glom in enumerate(id_non_connected):
                x_dist = np.fabs(glomeruli_pos[id_glom,0] - mossy_pos[:,0])
                z_dist = np.fabs(glomeruli_pos[id_glom,2]- mossy_pos[:,2])
                #Select the candidates mossy fibers
                roll_x = np.random.normal(self.x_length/2, self.x_sigma, 1)
                roll_z = np.random.normal(self.z_length/2, self.z_sigma, 1)
                candidates_x =  x_dist < roll_x
                #candidates_x = np.nonzero(candidates_x)
                candidates_z =  z_dist < roll_z
                #candidates_z = np.nonzero(candidates_z)
                final_candidates = (candidates_x & candidates_z)
                final_candidates = np.nonzero(final_candidates)
                #print(final_candidates)
                #print(random_ids[j])
                pre_locs[ptr,0] = random_ids[j]
                #print(id_glom)
                post_locs[ptr,0] = id_glom
                free_glomeruli[id_glom] = False
                ptr += 1
                #print(cluster[random_ids[j]])
                cluster[random_ids[j]] = np.append(cluster[random_ids[j]],id_glom)

        #Connect mossy fibers to glomeruli
        #self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
        print("Numero clusters:",len(cluster))
        for cl in cluster:
            print(len(cl))
        #print(cluster[0])
       
        #NUOVA SOLUZIONE -> FARE IL CONTRARIO: ASSOCIARE AD OGNI GLOM UNA MOSSY A CASO TRA QUELLE CHE VANNO BENE
        rng_id = np.random.randint(low=0,high=len(mossy_pos),size=len(glomeruli_pos))
        print(rng_id)
        pre_locs[0:len(glomeruli_pos),0] = rng_id
        post_locs[0:len(glomeruli_pos),0] = np.arange(0,len(glomeruli_pos),dtype=int)
        print(np.arange(0,len(glomeruli_pos),dtype=int))
        self.connect_cells(pre_ps, post_ps, pre_locs[:len(glomeruli_pos)], post_locs[:len(glomeruli_pos)])

        #for i, mossy in enumerate(glomeruli_pos):
            




        """ptr = 0
        n_granule = len(granule_pos)
        n_conn = self.convergence * n_granule
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        for granule_id,granule in enumerate(granule_pos):
            #TODO USE  iter_morphologies() 
            morpho = self.intermediate.cell_types[0].get_morphologies()[0].load().branches
            print("Number of branches in the morphology:",len(morpho))
            #morpho = self.intermediate.get_morphologies()[i].load()
            #In the morphology file the dednrites correspond to branches 3-7.
            dendrites = morpho[3:7]
            print("Number of dendritic branches",len(dendrites))
            clusters_id = np.arange(len(cluster), dtype=int)
            #print("Numero clusters:",len(cluster))
            #print(clusters_id)
            np.random.shuffle(clusters_id)
            #print(clusters_id)
            connected_glomeruli = 0
            #For each cluster find the candidate gloms to connect
            for id in clusters_id:
                #print("GLOM IN CLUSTER ", i ," : ",len(cluster[id]))
                cluster_dist = np.sqrt(
                np.power(granule[0] - glomeruli_pos[cluster[id],0], 2)
                + np.power(granule[1] - glomeruli_pos[cluster[id],1], 2)
                + np.power(granule[2] - glomeruli_pos[cluster[id],2], 2))
                
                #print(cluster_dist)
                cluster_dist_bool = cluster_dist < self.glom_granule_max_radius
                #print(cluster_dist_bool)
                cluster_dist_ids = np.nonzero(cluster_dist_bool)[0]
                if (len(cluster_dist_ids)>0):
                    print(cluster_dist_ids)
                    #Select a random granule between the good ones
                    choice = np.random.randint(low=0,high=len(cluster_dist_ids), dtype=int)
                    glomerulus_id = cluster_dist_ids[choice]
                    #print("Chosen id:", glomerulus_id)
                    #branch_id = np.where(morpho == dendrites[connected_glomeruli])
                    branch_id = connected_glomeruli
                    point_on_branch_id = len(dendrites[connected_glomeruli])
                    pre_locs[ptr,0] = glomerulus_id
                    post_locs[ptr,0] = granule_id
                    post_locs[ptr,1] = branch_id
                    post_locs[ptr,2] = point_on_branch_id
                    ptr += 1
                    connected_glomeruli += 1
            
                #If 4 gloms are connected to a granule cell, exit the cluster loop.
                if (connected_glomeruli == 4):
                    break
            print(connected_glomeruli)
        self.connect_cells(post_ps, self.intermediate.cell_types[0], pre_locs[:ptr], post_locs[:ptr])"""
