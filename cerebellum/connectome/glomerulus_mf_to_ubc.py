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
class ConnectomeGlomerulusMF_to_UBC(ConnectionStrategy):
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

            if (distance <self.max_radius):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks
    
    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        glom_pos = pre_ps.load_positions()
        ubc_pos = post_ps.load_positions()
        print(pre_ps.get_loaded_chunks())

        n_glom = len(glom_pos)
        n_ubc = len(ubc_pos)

        #Sort the chunks in increasing chunk_id order
        post_chunks = list(post_ps.get_loaded_chunks())
        chunks_ids = [c.id for c in post_chunks]
        sorted_ids = np.argsort(chunks_ids)
        sorted_chunks = []
        for i in sorted_ids:
            sorted_chunks.append(post_chunks[i])

        #Get the ubc_glom to ubc connectivity
        cs = self.scaffold.get_connectivity_set("glomerulus_ubc_to_ubc")
        iter = cs.load_connections().to(sorted_chunks)
        gid = iter.all()
        gid_global = gid[1]
        connected_ubc_ids = gid_global[:,0]
        avaiable_ubc = np.full((n_ubc), True, dtype=bool)
        avaiable_ubc[connected_ubc_ids] = False

        n_avaiable_ubc = np.count_nonzero(avaiable_ubc)
        max_connections = 1
        n_conn = n_avaiable_ubc * max_connections
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        ptr = 0

        connected_gloms = np.full((n_glom), False, dtype=bool)

        ubc_ids = np.arange(len(ubc_pos))   
        gloms_ids = np.arange(len(glom_pos))  
        ubc_ids = ubc_ids[avaiable_ubc]       
        selected_ubc = ubc_pos[avaiable_ubc]
        print(gloms_ids)
        
        for id,ubc in enumerate(selected_ubc):
            distances = np.linalg.norm(ubc-glom_pos[~connected_gloms],axis=1)
            avaiable_gloms_ids = np.where(~connected_gloms)[0]
            sorted_distances = np.sort(distances)
            id_sorted_distances = np.argsort(distances)
            mask = sorted_distances < self.max_radius
            if (np.count_nonzero(mask) > 0):
                post_locs[id] = ubc_ids[id]
                pre_locs[id] = avaiable_gloms_ids[id_sorted_distances[mask][0]]
                connected_gloms[avaiable_gloms_ids[id_sorted_distances[mask][0]]] = True
                ptr = ptr +1

        
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
