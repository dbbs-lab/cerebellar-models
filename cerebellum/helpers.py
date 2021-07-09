import numpy as np 
import os 
import math 

cell_properties = {'stellate_cell': {'radius': 4.0, 'id': 5, 'type': 'inh'}, 'basket_cell': {'radius': 6.0, 'id': 6, 'type': 'inh'}, 'granule_cell': {'radius': 2.5,'id': 2, 'type': 'exc'}, 'golgi_cell': {'radius': 8.0, 'id': 1, 'type': 'inh'}, 'glomerulus': {'radius': 1.5, 'density': 3e-4},'purkinje_cell': {'radius': 7.5}}
cell_per_layer = {'Stellate layer': ['stellate_cell'],'Basket layer': ['basket_cell'], 'granular layer': ['granule_cell','golgi_cell', 'glomerulus']}

def relative_folder(path):
    return os.path.join(os.path.dirname("__file__"), path)

def extract_voxels(indexes, number, dim, threshold):
    """
        It prepares voxels and cell count to be fed to the particle placement.
        It extracts the voxels (with 3D origin and size), the neuron count in each voxel.

        :param indexes: 3D indexes of voxels belonging to a certain layer from the full voxelized brain/cerebellum/cerebellum part
        :param number: number of neurons in the selected voxels  and a possible dimension along which we want to cut to
        place only in a sub-part of a certain region

    """

    voxel_sizes = np.array([25.,25.,25.])
    voxels = []
    neu_in_voxels = []

    indexes = np.transpose(indexes)
    #print("transpose size",indexes.shape)
    indexes = np.asarray(indexes*int(voxel_sizes[0]))

    for ind in range(len(indexes)):
        if indexes[ind][dim]<threshold:
            voxels.append([np.array(indexes[ind]), voxel_sizes])
            neu_in_voxels.append(number[ind])

    return voxels, neu_in_voxels
