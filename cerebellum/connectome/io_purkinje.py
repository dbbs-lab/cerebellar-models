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
class ConnectomeIOPurkinje(ConnectionStrategy):

    def get_region_of_interest(self, chunk):

        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        return chunks

    #We need to override the default _get_connect_args_from_job function because 
    #We need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, chunk)
        return pre, post

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps.cell_type, pre_ps, post_ps.cell_type, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
              
        io_pos = pre_ps.load_positions()
        pc_pos = post_ps.load_positions()
        n_io = len(io_pos)
        n_pc = len(pc_pos)

        print("N io cells:", n_io)
        print("N pcs:", n_pc)

        #Each PC is connect to a single IO cell
        pre_locs = np.full((n_pc, 3), -1, dtype=int)
        post_locs = np.full((n_pc, 3), -1, dtype=int)+

        #We connect each pc to an IO cell
        
        selected_io_ids = np.random.randint(0, high=n_io, size=n_pc, dtype=int)
        
        for j,_ in enumerate(pc_pos):
            
            pre_locs[j,0] = selected_io_ids[j]
            post_locs[j,0] = j

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
        print("Connected", n_io, "IO cells to", n_pc, "PCs.")

