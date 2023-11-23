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
class ConnectomeMossyDCNP(ConnectionStrategy):
    convergence = config.attr(type=int, required=True)

    #We need to connect a DCNp to #convergence MFs from the whole region,
    #so we just return all the chunks
    def get_region_of_interest(self, chunk):
        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        return chunks

    #We need to override the default _get_connect_args_from_job function because 
    #We need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    #In this way if for each DCNP in the ROI we select randomly the mossy fiber
    #to connect according to the convergence parameter given by the user
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, roi)
        post = HemitypeCollection(self.postsynaptic, chunk)
        return pre, post

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps.cell_type, pre_ps, post_ps.cell_type, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
              
        mossy_pos = pre_ps.load_positions()
        dcnp_pos = post_ps.load_positions()
        n_mossy = len(mossy_pos)
        n_dcnp = len(dcnp_pos)

        print("N mossy:", n_mossy)
        print("N dcnp:", n_dcnp)

        max_synapses = n_dcnp * self.convergence
        pre_locs = np.full((max_synapses, 3), -1, dtype=int)
        post_locs = np.full((max_synapses, 3), -1, dtype=int)

        ptr = 0

        #We connect each DCNp to #convergence MFs
        for j,_ in enumerate(dcnp_pos):

            #Select randomly #convergence MFs from all the MFs
            selected_mf_ids = np.random.choice(n_mossy, self.convergence, replace=False)
            pre_locs[ptr:ptr+self.convergence,0] = selected_mf_ids
            post_locs[ptr:ptr+self.convergence,0] = j
            ptr = ptr + self.convergence

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
        print("Connected", n_dcnp, "dcnp to", n_mossy, "mossy fibers.")

