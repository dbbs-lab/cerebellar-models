import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.connectivity.strategy import HemitypeCollection
from bsb.storage import Chunk
from bsb import config
from bsb.cell_types import CellType
from itertools import chain
from bsb.reporting import report
from scipy.stats.distributions import truncexpon
from bsb.connectivity.strategy import Hemitype

@config.node
class ConnectomeIO_MLI(ConnectionStrategy):
    mli_pc_connectivity = config.attr(type=str, required=True)
    io_pc_connectivity = config.attr(type=str, required=True)

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
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)
    
    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        io_pos = pre_ps.load_positions()
        mli_pos = post_ps.load_positions()


        #All the MLIs connected to the same PC must be connected to the IO cell which
        #is connected to that PC

        #We retrieve the connectivity data for the mli-pc connectivity

        cs = self.scaffold.get_connectivity_set(self.mli_pc_connectivity)
        iter = cs.load_connections().as_globals()
        mli_per_pc_list = []
        
        gid_global,gid_local = iter.all()
        unique_pc = np.unique(gid_local,axis=0)

        for current in unique_pc:
            pc_ids = np.where((gid_local[:,0]==current[0]))
            starting = np.min(gid_global[:,0])
            mli_per_pc_list.append(gid_global[pc_ids[0],0]-starting)


        #We retrieve the connectivity data for the io-pc connectivity

        cs = self.scaffold.get_connectivity_set(self.io_pc_connectivity)
        iter = cs.load_connections().as_globals()
        io_pc_list = []
        
        gid_global,gid_local = iter.all()
        unique_pc = np.unique(gid_local,axis=0)

        for current in unique_pc:
            pc_ids = np.where((gid_local[:,0]==current[0]))
            starting = np.min(gid_global[:,0])
            io_pc_list.append(gid_global[pc_ids[0],0]-starting)

        n_purkinje = len(io_pc_list)
        max_connections = mli_pos * n_purkinje
        pre_locs = np.full((max_connections, 3), -1, dtype=int)
        post_locs = np.full((max_connections, 3), -1, dtype=int)


        ptr = 0

        for j,id_io in enumerate(io_pc_list):
            for k,id_mli in enumerate(mli_per_pc_list[j]):

                pre_locs[ptr,0] = id_io
                post_locs[ptr,0] = id_mli
                ptr = ptr + 1
        
        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)

        print("Connected", len(io_pos), "io cells to", len(mli_pos), "mlis.")

        
