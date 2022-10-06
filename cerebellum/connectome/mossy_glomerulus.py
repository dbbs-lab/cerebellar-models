from types import CellType
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.connectivity.strategy import HemitypeCollection
from bsb.storage import Chunk
from bsb import config
from bsb.cell_types import CellType
from itertools import chain
from bsb.reporting import report
from scipy.stats.distributions import truncexpon

@config.node
class ConnectomeMossyGlomerulus(ConnectionStrategy):
    x_length = config.attr(type=int, required=True)
    z_length = config.attr(type=int, required=True)

    #

    def get_region_of_interest(self, chunk):

        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        
        #We look for chunks satisfying the geometrical conditions:
        #a mossy fiber can be connected to a glomuerulus
        #inside a self.x_length x self.z_length rectangle
        #centered in the mossy
        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[2]- c[2]) 
            x_dist = x_dist * chunk.dimensions[0]
            z_dist = z_dist * chunk.dimensions[2]

            if (x_dist < self.x_length/2 and  z_dist < self.z_length/2):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        #selected_chunks = [chunk]
        return selected_chunks

    #We need to override the default _get_connect_args_from_job function because 
    #we make implicit use of 1:20 ratio between mossy fibers and glomeruli 
    #in the custom connection strategy.
    #We need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    #In this way if for each glomerulus in the ROI we select randomly the mossy fiber
    #to connect, the divergence ratio is automatically around 20.
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, [chunk])
        return pre, post

    def connect(self, pre, post):
        #pre_type = pre.cell_types[0]
        #post_type = post.cell_types[0]
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
              
        mossy_pos = pre_ps.load_positions()
        glomeruli_pos = post_ps.load_positions()
        n_mossy = len(mossy_pos)
        n_glom = len(glomeruli_pos)

        print("Pre",pre_ps.get_loaded_chunks())
        print("Post",post_ps.get_loaded_chunks())

        print("N mossy:", n_mossy)
        print("N gloms:", n_glom)

        #We work assuming that there is at least 1 mossy fiber in the ROI.
        #Otherwise there is something wrong in the placement phase.
        pre_locs = np.full((n_glom, 3), -1, dtype=int)
        post_locs = np.full((n_glom, 3), -1, dtype=int)

        #We generate the ids of the MF tp be connected before
        #to loop among the granule to speed up the computations.
        #We use a truncated exponential distribution to favour
        #the MF closer o the glomerulus.
        exp_dist = truncexpon(b=0.99)
        rolls = np.floor(n_mossy*exp_dist.rvs(size=n_glom))
        rolls = rolls.astype(int)

        #We connect each glomerulus to a mossy fiber.
        for j,glomerulus in enumerate(glomeruli_pos):
            
            dist = np.linalg.norm(glomerulus-mossy_pos,axis=1)
            id_sorted_dist = np.argsort(dist)
            #MF closer to the granule have an higher chance to be picked 
            pre_locs[j,0] = id_sorted_dist[rolls[j]]
            post_locs[j,0] = j

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
        print("Connected", n_glom, "glomeruli to", n_mossy, "mossy fibers.")
        #print(pre_locs)
        #print(post_locs)


