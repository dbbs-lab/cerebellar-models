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

@config.node
class ConnectomeGlomerulusUBC_to_UBC(ConnectionStrategy):
    max_radius = config.attr(type=float, required=True)

    #We need to override the default _get_connect_args_from_job function because 
    #we make implicit use of 1:4 ratio between ubc and ubc_gloms
    #in the custom connection strategy.
    #We need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    #In this way if for each glomerulus in the ROI we select randomly the ubc
    #to connect, the divergence ratio is automatically around 4.
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, chunk)
        return pre, post

    def get_region_of_interest(self, chunk):
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []
        
        #We look for chunks satisfying the geometrical conditions:
        #a ubc can be connected to a mf/ubc glom if it is less than 50 um away
        current_chunk_coord = chunk*chunk.dimensions

        for c in chunks:
            c_coord = c*c.dimensions
            distance = np.linalg.norm(current_chunk_coord-c_coord)

            if (distance < self.max_radius):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks
    
    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)



    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        glom_pos = pre_ps.load_positions()
        ubc_pos = post_ps.load_positions()
        id_ubc = np.arange(len(ubc_pos))
        np.random.permutation(id_ubc)
        one_third = len(ubc_pos) // 3
        selected_ubc = ubc_pos[one_third:]

        n_glom = len(glom_pos)
        n_ubc = len(selected_ubc)
        max_connections = 1
        n_conn = n_ubc * max_connections
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        ptr = 0

        connected_gloms = np.full((n_glom), False, dtype=bool)
        
        for id,ubc in enumerate(selected_ubc):
            distances = np.linalg.norm(ubc-glom_pos[~connected_gloms],axis=1)
            sorted_distances = np.sort(distances)
            id_sorted_distances = np.argsort(distances)
            mask = sorted_distances < self.max_radius
            post_locs[id] = id_ubc[id]
            pre_locs[id] = id_sorted_distances[mask][0]
            connected_gloms[id_sorted_distances[mask][0]] = True
            ptr = ptr +1

        
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
