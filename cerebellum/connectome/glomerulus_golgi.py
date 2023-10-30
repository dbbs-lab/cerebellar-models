import numpy as np
import itertools
from bsb.connectivity import ConnectionStrategy
from bsb.storage import Chunk
from bsb import config
from scipy.stats.distributions import truncexpon
from bsb.morphologies import Morphology
from bsb.storage.interfaces import ConnectivitySet as IConnectivitySet
from bsb.connectivity.strategy import HemitypeCollection


@config.node
class ConnectomeGlomerulusGolgi(ConnectionStrategy):
    radius = config.attr(type=int, required=True)
    affinity = config.attr(type=float, default=1.0)

    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, chunk)
        post = HemitypeCollection(self.postsynaptic, roi)
        return pre, post

    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []

        # Look for chunks which are less than radius away from the current one.
        for c in chunks:
            dist = np.sqrt(
                np.power(chunk[0] * chunk.dimensions[0] - c[0] * c.dimensions[0], 2)
                + np.power(chunk[1] * chunk.dimensions[1] - c[1] * c.dimensions[1], 2)
                + np.power(chunk[2] * chunk.dimensions[2] - c[2] * c.dimensions[2], 2)
            )

            #If c is intersecting the granular layer partition and if it is close enough to chunk, then it is in the ROI
            if dist < self.radius and c[1]*c.dimensions[1] < self.scaffold.partitions.granular_layer.thickness:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks

    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
        
        # If synaptic contacts need to be made we use this exponential distribution
        # to pick the closer by subcell_labels.
        exp_dist = truncexpon(b=5, scale=0.03)
        
        glomeruli_pos = pre_ps.load_positions()
        golgi_pos = post_ps.load_positions()
        n_glom = len(glomeruli_pos)
        n_golgi = len(golgi_pos)
        n_conn = n_glom * n_golgi
        pre_locs = np.full((n_conn, 3), -1, dtype=int)
        post_locs = np.full((n_conn, 3), -1, dtype=int)

        #Cache morphologies and generate the morphologies iterator
        morpho_set = post_ps.load_morphologies()
        golgi_morphos_itr = morpho_set.iter_morphologies(cache=True, hard_cache=True)
        golgi_morphos = []

        for m in golgi_morphos_itr:
            golgi_morphos.append(m)

        #Find cells to connect
        ptr = 0
        
        for id_glom,glom_p in enumerate(glomeruli_pos):
        
            dist = np.sqrt(
                np.power(glom_p[0] - golgi_pos[:, 0], 2)
                + np.power(glom_p[1] - golgi_pos[:, 1], 2)
                + np.power(glom_p[2] - golgi_pos[:, 2], 2)
            )
                
            to_connect_bool = dist < self.radius
            to_connect_idx = np.nonzero(to_connect_bool)[0]
            connected_golgi = len(to_connect_idx)
            pre_locs[ptr : (ptr + connected_golgi), 0] = id_glom
            post_locs[ptr : (ptr + connected_golgi), 0] = to_connect_idx
            
            ptr_golgi = 0
            for counter, id in enumerate(to_connect_idx):
                basal_dendrides_branches = golgi_morphos[id].get_branches()
                            
                #Get the starting branch id of the dendritic branches
                first_dendride_id = golgi_morphos[id].branches.index(basal_dendrides_branches[0])
                
                #Find terminal points on branches
                terminal_ids = np.full(len(basal_dendrides_branches), 0, dtype=int)
                for i,b in enumerate(basal_dendrides_branches):
                    if b.is_terminal:
                        terminal_ids[i] = 1
                terminal_branches_ids = np.nonzero(terminal_ids)[0]

                #Keep only terminal branches
                basal_dendrides_branches = np.take(basal_dendrides_branches, terminal_branches_ids, axis=0)
                terminal_branches_ids = terminal_branches_ids + first_dendride_id

                #Find the point-on-branch ids of the tips
                tips_coordinates = np.full((len(basal_dendrides_branches),3), 0, dtype=float)
                for i,branch in enumerate(basal_dendrides_branches):
                    tips_coordinates[i] = branch.points[-1]

                #Choose randomly the branch where the synapse is made
                #favouring the branches closer to the glomerulus.
                rolls = exp_dist.rvs(size=len(basal_dendrides_branches))
                
                # Compute the distance between terminal points of basal dendrites 
                # and the soma of the avaiable glomeruli
                pts_dist = np.sqrt(np.power(tips_coordinates[:,0] + golgi_pos[id][0] - glom_p[0], 2)
                            + np.power(tips_coordinates[:,1] + golgi_pos[id][1] - glom_p[1], 2)
                            + np.power(tips_coordinates[:,2] + golgi_pos[id][2] - glom_p[2], 2)
                        )

                sorted_pts_ids = np.argsort(pts_dist)
                # Pick the point in which we form a synapse according to a exponential distribution mapped
                # through the distance indices: high chance to pick closeby points.
                pt_idx = sorted_pts_ids[int(len(basal_dendrides_branches)*rolls[np.random.randint(0,len(rolls))])]

                #The id of the branch is the id of the terminal_branches plus the id of the first dendritic branch
                post_locs[ptr+counter,1] = terminal_branches_ids[pt_idx]
                #We connect the tip of the branch
                post_locs[ptr+counter,2] = len(basal_dendrides_branches[pt_idx].points)-1
 
            ptr += connected_golgi

        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
