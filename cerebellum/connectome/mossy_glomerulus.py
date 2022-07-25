import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config

@config.node
class ConnectomeMossyGlomerulus(ConnectionStrategy):
    x_length = config.attr(type=int, required=True)
    z_length = config.attr(type=int, required=True)
    x_sigma = config.attr(type=int, required=True)
    z_sigma = config.attr(type=int, required=True)
    mean = config.attr(type=int, required=True)
    sigma = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        # Look for chunks which are less than radius away
        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[2]- c[2])

        x_dist = x_dist * chunk.dimensions[0]
        z_dist = z_dist * chunk.dimensions[3]

        #TODO vectorize this part, removing the loop
        for c in chunks:
            if x_dist < self.x_length/2 and  z_dist < self.z_length/2:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        pre_type = pre.cell_types[0]
        post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
     
        mossy_pos = pre_ps.load_positions()
        glomeruli_pos = post_ps.load_positions()
        n_mossy = len(mossy_pos)
        n_glom = len(glomeruli_pos)
        n_conn = n_glom * n_mossy
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)
        # Since a glomerulus can be connected only to one mossy fiber, we
        # need to keep track of the connected glomeruli
        free_glomeruli = np.full(n_glom, True)

        # Find glomeruli to connect            
        ptr = 0
        for i, mossy in enumerate(mossy_pos):
            x_dist = np.fabs(mossy[0] - glomeruli_pos[:,0])
            z_dist = np.fabs(mossy[2]- glomeruli_pos[:,2])

            #Select the candidates glomeruli using a distance-based probability rule
            roll_x = np.random.normal(self.x_length/2, self.x_sigma, n_glom)
            roll_z = np.random.normal(self.z_length/2, self.z_sigma, n_glom)
            candidates_x =  x_dist < roll_x
            candidates_x = np.nonzero(candidates_x)
            candidates_z =  z_dist < roll_z
            candidates_z = np.nonzero(candidates_z)
            final_candidates = (candidates_x & candidates_z) & free_glomeruli
            final_candidates = np.nonzero(final_candidates, dtype = int)
            
            #Draw the number of glomeruli to connect to the current mossy fiber
            # following a Gaussian distribution
            glomeruli_to_connect = np.random.normal(self.mean, self.sigma)

            #If the good glomeruli are less than glomeruli_to_connect, connect them all;
            # Otherwise, select the first good glomeruli_to_connect glomeruli
            num_candidates = len(final_candidates)
            if (num_candidates < glomeruli_to_connect):
                pre_locs[ptr:ptr+num_candidates] = i
                post_locs[ptr:ptr+num_candidates] = final_candidates
                free_glomeruli[final_candidates] = False
                ptr += num_candidates
            else:
                pre_locs[ptr:ptr+glomeruli_to_connect] = i
                post_locs[ptr:ptr+glomeruli_to_connect] = final_candidates[:glomeruli_to_connect]
                free_glomeruli[final_candidates[:glomeruli_to_connect]] = False
                ptr += glomeruli_to_connect
        
        #If there are some glomeruli which are not connected, connect them to a random near mossy fiber
        id_non_connected = np.nonzero(~free_glomeruli)
        #Draw some random numbers to be used as the indices of mossy fibers to connect to
        # each non connected glomerulus.
        random_ids = np.random.randint(0, len(id_non_connected), len(id_non_connected), dtype=int)
        if (len(id_non_connected)>0):
            for id_glom in id_non_connected:
                x_dist = np.fabs(mossy[id_glom][0] - mossy_pos[:,0])
                z_dist = np.fabs(mossy[id_glom][2]- mossy_pos[:,2])
                #Select the candidates mossy fibers
                roll_x = np.random.normal(self.x_length/2, self.x_sigma, n_glom)
                roll_z = np.random.normal(self.z_length/2, self.z_sigma, n_glom)
                candidates_x =  x_dist < roll_x
                candidates_x = np.nonzero(candidates_x)
                candidates_z =  z_dist < roll_z
                candidates_z = np.nonzero(candidates_z)
                final_candidates = (candidates_x & candidates_z)
                final_candidates = np.nonzero(final_candidates, dtype = int)
                pre_locs[ptr:ptr+1] = final_candidates[random_ids[id_glom]]
                post_locs[ptr:ptr+1] = id_glom
                free_glomeruli[id_glom] = False
                ptr += 1
                
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])