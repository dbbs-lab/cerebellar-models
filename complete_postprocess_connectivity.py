from bsb.core import from_storage
import numpy as np

network = from_storage("MiniHuman.hdf5")

for ps in network.get_connectivity_sets():
    print("----------------------------")
    print(ps.tag)
    cs= network.get_connectivity_set(ps.tag).load_connections().as_globals()
    data=cs.all()
    data_pre = data[0]
    data_post = data[1]

    pre_locs = data_pre
    post_locs = data_post

    """pre_cells, pre_counts = np.unique(pre_locs[:, 0], axis=0, return_counts=True)
    post_cells, post_counts = np.unique(post_locs[:, 0], axis=0, return_counts=True)
    print("tot synapses:", np.sum(pre_counts))
    combos, combo_counts = np.unique(np.column_stack((pre_locs[:, 0], post_locs[:, 0])), axis=0, return_counts=True)

    uniquePost= np.unique(combos[:,:], axis=1, return_counts=True)
    uniquePre= np.unique(combos[:,:], axis=0, return_counts=True)
    print("mean divergence:", len(combos)/len(uniquePre[0]))
    print("mean std:", np.std(pre_counts))

    print("mean convergence:", len(combos)/len(uniquePost[0]))
    print("mean std:", np.std(post_counts))

    print("mean synapse/pair:", np.sum(pre_counts)/len(combos))
    print("------------------------")
    """

    pre_cells, pre_counts = np.unique(pre_locs[:, 0], axis=0, return_counts=True)
    post_cells, post_counts = np.unique(post_locs[:, 0], axis=0, return_counts=True)
    print("tot synapses:", np.sum(pre_counts))
    combos, combo_counts = np.unique(np.column_stack((pre_locs[:, 0], post_locs[:, 0])), axis=0, return_counts=True)
    print("Numero di sinapsi:",np.sum(combo_counts))

    list_pre = []
    for i,c in enumerate(pre_cells):
        sel = combos[combos[:,0] == c]
        list_pre.append(len(sel))
        #print(len(sel))

    #print(max(list_pre))
    #print(min(list_pre))
    print("Divergence:", np.sum(list_pre)/len(list_pre)," pm ", np.sqrt((np.max(list_pre)-np.min(list_pre))))

    list_post = []
    for i,c in enumerate(post_cells):
        sel = combos[combos[:,1] == c]
        list_post.append(len(sel))

    #print(max(list_post))
    #print(min(list_post))
    #print(list_post)
    print("Convergence:", np.sum(list_post)/len(list_post)," pm ", np.sqrt((np.max(list_post)-np.min(list_post))))


    """uniquePost= np.unique(combos[:,1:], axis=0, return_counts=True)
    uniquePre= np.unique(combos[:,:1], axis=0, return_counts=True)
    print("mean divergence:", len(combos)/len(uniquePre[0]))
    print("mean std:", np.std(pre_counts))

    print("mean convergence:", len(combos)/len(uniquePost[0]))
    print("mean std:", np.std(post_counts))"""

    print("mean synapse/pair:", np.sum(pre_counts)/len(combos))

