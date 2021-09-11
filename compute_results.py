# flake8: noqa
import nrrd
import numpy as np
from JSONread import *
import plotly.graph_objects as go
import plotly.express as px
import chart_studio.plotly as py
from plotly.subplots import make_subplots
import random
import matplotlib.pyplot as plt
import statistics

PLOTTING = True

data_path = "../Flocculus3.0_Lingula/"

# Literature info on GR layer inhibitory neurons
# Number of Lugaro cells according to ratio from Dieudonne et al (2011)
# 1.0/15.0 * number of Purkinje cells
fac_Lugaro = 15.0

# Density of UBC according to Sekerkova et al (2014)
# 1.5*4*1000 = 6000 UBC/mm^3 for Lingula, Central Lobule, Culmen, Declive, Follum-tuber vermis, Pyramus, Simple lobule, Ansiform lobule
# 6.4*4*1000 = 25600 UBC/mm^3 for Uvula, Paramedian lobule, Copula pyramidis
# 24.6*4*1000 = 98400 UBC/ mm^3 for Nodulus
# 19.5*4*1000 = 78000 UBC/mm^3 for Flocculus
# 18272.36 UBC/mm^3 for Paraflocculus TODO: Redo computation after annotation corrections
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

jsontextfile = open(data_path + "data/brain_regions.json", "r")
jsoncontent = json.loads(jsontextfile.read())
search_children(jsoncontent["msg"][0], 0)

region_name = "Flocculus"
# region_name = "Lingula (I)"

id_region = []
id_mol = -1
id_pc = -1
id_gr = -1

for id_reg, name in id_to_region_dictionary_ALLNAME.items():
    if region_name in name:
        if region_name + "," in name:
            id_region.append(id_reg)
            if "molecular layer" in name:
                id_mol = id_reg
            elif "granular layer" in name:
                id_gr = id_reg
            elif "Purkinje layer" in name:
                id_pc = id_reg
        else:
            id_current_region = id_reg

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

mask[region_name] = ann == id_current_region

region_names.append(region_name)
number_cells.append(np.sum(number_cells))
number_neurons.append(np.sum(number_neurons))
number_inhibitory.append(np.sum(number_inhibitory))
volumes.append(np.sum(volumes))
cell_densities.append(number_cells[-1] / volumes[-1])
neuron_densities.append(number_neurons[-1] / volumes[-1])

layers_per_cell = {
    "granule": "granular layer",
    "golgi": "granular layer",
    "purkinje": "Purkinje layer",
    "stellate": "Stellate layer",
    "basket": "Basket layer",
}


i_granular = -1
i_molecular = -1
i_purkinje = -1
print(region_names)
for i in range(len(region_names)):
    if "granular layer" in region_names[i]:
        i_granular = i
    if "molecular layer" in region_names[i]:
        i_molecular = i
    if "Purkinje layer" in region_names[i]:
        i_purkinje = i

print("i_granular: ", i_granular)
print("i_molecular: ", i_molecular)
print("Purkinje layer: ", i_purkinje)
# Volumes print
print("Volume molecular layer: " + str(volumes[i_molecular]))
print("Volume purkinje layer: " + str(volumes[i_purkinje]))
print("Volume granular layer: " + str(volumes[i_granular]))


# Number of basket and stellate cells according to relative layer distance
bounds, h = nrrd.read(data_path + "data/boundaries_mo.nrrd")
thickness_ratio = 2.0 / 3.0  # ratio of molecular layer space for stellate cells
up_layer_distance = np.abs(bounds[0])
down_layer_distance = np.abs(bounds[1])
print(np.amax(up_layer_distance))
print(np.amin(up_layer_distance))
print(up_layer_distance[up_layer_distance != 0])
relative_layer_distance = np.zeros(up_layer_distance.shape)

mask_mol = ann == id_mol

up_layer_distance = up_layer_distance[mask_mol]
down_layer_distance = down_layer_distance[mask_mol]
print("Mean down layer: ", statistics.mean(down_layer_distance))
print("Mean up layer: ", statistics.mean(up_layer_distance))
relative_layer_distance[mask_mol] = up_layer_distance / (up_layer_distance + down_layer_distance)


mask_of_stellate = relative_layer_distance * mask_mol
mask_of_stellate = (mask_of_stellate > 0) * (mask_of_stellate < thickness_ratio)
mask["Stellate layer"] = mask_of_stellate
mask_of_basket = ~mask_of_stellate * mask_mol
mask["Basket layer"] = mask_of_basket
vox_in_layer["Stellate layer"] = len(np.where(mask["Stellate layer"])[0])
vox_in_layer["Basket layer"] = len(np.where(mask["Basket layer"])[0])

volumes.append(len(np.where(mask["Stellate layer"])[0]) / (4.0 ** 3) / 1000.0)  # in mm3
volumes.append(len(np.where(mask["Basket layer"])[0]) / (4.0 ** 3) / 1000.0)  # in mm3
print("Volume Stellate layer: " + str(volumes[4]))
print("Volume Basket layer: " + str(volumes[5]))
print("Total volume: " + str(volumes[3]))

print("Voxels Granular layer: " + str(vox_in_layer["granular layer"]))
print("Voxels Purkinje layer: " + str(vox_in_layer["Purkinje layer"]))
print("Voxels Stellate layer: " + str(vox_in_layer["Stellate layer"]))
print("Voxels Basket layer: " + str(vox_in_layer["Basket layer"]))


sliding_dir = 2

# PC layer
maskPC = ann == id_pc
maskPC = maskPC * 1
PCsurface = np.where(maskPC)

# ML
maskMLI = ann == id_mol
maskMLI = maskMLI * 1
MLIsurface = np.where(maskMLI)
SCsurface = np.where(mask_of_stellate)
BCsurface = np.where(mask_of_basket)

# GL
maskGL = ann == id_gr
maskGL = maskGL * 1
GLsurface = np.where(maskGL)

# Extract all region
mask_all = np.isin(ann, id_region)
region = (
    ann - id_region[2] + 1
)  # Scale to have granular layer = 1, PC layer = 2, molecular layer = 3
region[~mask_all] = 0  # Outside of the region = 0
region[mask_of_stellate] = 4  # To differentiate SC and BC in the ML

# Cut around the region
region_index = np.nonzero(region)
region = region[
    np.amin(region_index[0]) - 10 : np.amax(region_index[0]) + 10,
    np.amin(region_index[1]) - 10 : np.amax(region_index[1]) + 10,
    np.amin(region_index[2]) - 10 : np.amax(region_index[2]) + 10,
]
dim = region.shape
example_slice = np.take(region, 0, axis=sliding_dir)
r, c = example_slice.shape

# Define colorscale
color_region = [
    # External to region: gray
    [0, "rgb(220, 220, 220)"],
    [0.2, "rgb(220, 220, 220)"],
    # Granular layer: red
    [0.2, "rgb(256, 0, 0)"],
    [0.4, "rgb(256, 0, 0)"],
    # PC layer: green
    [0.4, "rgb(58, 146, 94)"],
    [0.6, "rgb(58, 146, 94)"],
    # Molecular layer - BC: arancione
    [0.6, "rgb(249, 90, 8)"],
    [0.8, "rgb(249, 90, 8)"],
    # Molecular layer - SC: giallo
    [0.8, "rgb(245, 177, 0)"],
    [1.0, "rgb(245, 177, 0)"],
]


# for sl in range(10, 20):    #ciclo for da 10 a 19
#    print(sl)


# 69 – 386 ( range 317 )
low = np.amin(region_index[2]) - 10  # 69
high = np.amax(region_index[2]) + 10  # 387
rang = high - low  # 317

reason = False
if reason:
    for j in range(rang):  # per ogni slice z
        matrix = np.take(region, j, axis=sliding_dir)
        count = 0
        num_zero = 0
        sum = 0
        sum_PC = 0

        for i in matrix:
            # per ogni riga
            num_zero = 0
            num_uno = 0
            num_due = 0
            num_tre = 0
            num_quattro = 0
            for k in i:  # per ogni elemento della riga
                if k == 0:
                    num_zero += 1
                if k == 1:
                    num_uno += 1
                if k == 2:
                    num_due += 1
                if k == 3:
                    num_tre += 1
                if k == 4:
                    num_quattro += 1
            sum += num_uno + num_due + num_tre + num_quattro
            # print(" Voxels outside:", num_zero, " | GL:", num_uno, " | PCL:", num_due, " BL:", num_tre, " | SL:", num_quattro, " --- TOT = ", num_uno+num_due+num_tre+num_quattro)
            sum_PC += num_due
        print("slice num:", j, "(", j + 70, ") - tot voxels", sum, "|| PCL: ", sum_PC)

fig_scatter6 = go.Figure(
    data=go.Heatmap(
        z=np.take(region, 39, axis=sliding_dir),
        colorscale=color_region
        # cmin=0, cmax=200,
        # colorbar=dict(thickness=20, ticklen=4)
    )
)
fig_scatter6.show()
