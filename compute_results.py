# flake8: noqa
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import chart_studio.plotly as py
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
import statistics
import random
import nrrd
import json

from cerebellum.brain_regions_tree import BrainRegionsTree

PLOTTING = True
data_path = "../Flocculus3.0_Lingula/"


reference_UBC_density = {
    "Lingula (I)": 6000,
    "Uvula": 25600,
    "Nodulus": 98400,
    "Flocculus": 78000,
    "Paraflocculus": 18272.36,
}

ann, h = nrrd.read(data_path + "data/annotations.nrrd")
print(h)
print(ann.shape)
print(ann[0, 0, 0])
print(ann[395, 191, 121])
dens_cell, h = nrrd.read(data_path + "data/cell_density.nrrd")
print("dens cell ", dens_cell.shape)
dens_neuron, h = nrrd.read(data_path + "data/neu_density.nrrd")
dens_inh, h = nrrd.read(data_path + "data/inh_density.nrrd")
orientations, h = nrrd.read(data_path + "data/orientations_cereb.nrrd")
print("orientation: ", orientations.shape)
print(orientations[0].shape)

# starting to create
brain_regions_tree = BrainRegionsTree()
id_to_region_dictionary = brain_regions_tree.id_to_region_dictionary_ALLNAME()
id_region = brain_regions_tree.id_region
id_gr, id_pc, id_mol = brain_regions_tree.get_id_gr_pc_mol()

# mask_current_region = annis in id_region
VOXEL_SIZE = 25.0  # um
region_names = []
number_cells = []
number_neurons = []
number_inhibitory = []
volumes = []
cell_densities = []
neuron_densities = []
print(id_region)
print(id_mol)
mask = {}
vox_in_layer = {}
for id_ in id_region:
    region_names.append(id_to_region_dictionary[id_])
    region, layer = id_to_region_dictionary[id_].split(", ")
    print(layer)
    mask[layer] = ann == id_
    vox_in_layer[layer] = len(np.where(mask[layer])[0])
    number_cells.append(np.round(np.sum(dens_cell[mask[layer]])))
    number_neurons.append(np.round(np.sum(dens_neuron[mask[layer]])))
    number_inhibitory.append(np.round(np.sum(dens_inh[mask[layer]])))
    volumes.append(len(np.where(mask[layer])[0]) / (4.0 ** 3) / 1000.0)  # in mm3
    cell_densities.append(number_cells[-1] / volumes[-1])
    neuron_densities.append(number_neurons[-1] / volumes[-1])
