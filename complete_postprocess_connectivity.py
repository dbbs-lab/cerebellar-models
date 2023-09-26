from bsb.core import from_storage
import numpy as np

#Bult the network from file
network = from_storage("200_200_mouse.hdf5")

#Loop through tge entire connectivity set
for ps in network.get_connectivity_sets():
    
    #Print the name of the connectivity strategy
    print(ps.tag)
    
    #Get the ConnectivityIterator for the current connectivity strategy
    cs = network.get_connectivity_set(ps.tag).load_connections().as_globals()
    
    #Get the arrays containing the connecivity data
    data = cs.all()

    #data[0] contains pre data; data[1] post data
    pre_locs = data[0]
    post_locs = data[1]

    #The total number of synapses is the number of entries in pre_locs or post locs
    print("Tot synapses:", len(pre_locs))
    
    #Find the pairs of pre-post neurons (combos) and count how many synapses there are between each pair (combo_counts)
    combos, combo_counts = np.unique(np.column_stack((pre_locs[:, 0], post_locs[:, 0])), axis=0, return_counts=True)
    print("Synapses per pair:", np.mean(combo_counts), "pm", np.std(combo_counts))

    #Find the unique pre neurons and compute how many post are connected to each of them
    uniquePre, uniquePre_count = np.unique(combos[:,0], axis=0, return_counts=True)
    print("Divergence:", np.mean(uniquePre_count), "pm", np.std(uniquePre_count))

    #Find the unique post neurons and compute how many pre are connected to each of them
    niquePost, uniquePost_count = np.unique(combos[:,1], axis=0, return_counts=True)
    print("mean convergence:", np.mean(uniquePost_count), "pm", np.std(uniquePost_count))

    print("-------------------------------------------------------------")
