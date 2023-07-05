from bsb.core import from_storage
import numpy as np

network = from_storage("HC.hdf5")

cs= network.get_connectivity_set('ascending_axon_to_golgi').load_connections().as_globals()
data=cs.all()
data_pre = data[0]
data_post = data[1]


pre_locs = data_pre
post_locs = data_post

pre_cells, pre_counts = np.unique(pre_locs[:, 0], axis=0, return_counts=True)
print("tot synapses:", np.sum(pre_counts))
combos, combo_counts = np.unique(np.column_stack((pre_locs[:, 0], post_locs[:, 0])), axis=0, return_counts=True)
uniquePre= np.unique(combos[:,:1], axis=0, return_counts=True)
print("mean divergence:", len(combos)/len(uniquePre[0]))
uniquePost= np.unique(combos[:,1:], axis=0, return_counts=True)
print("mean convergence:", len(combos)/len(uniquePost[0]))
print("mean synapse/pair:", np.sum(pre_counts)/len(combos))
