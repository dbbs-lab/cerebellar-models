# flake8: noqa
# Importing scaffold
from scaffold.particles import (
    Particle,
    ParticleSystem,
    plot_particle_system,
    plot_detailed_system,
)
from scaffold.placement import ParallelArrayPlacement
from scaffold.config import JSONConfig
from scaffold.core import Scaffold, from_hdf5
from scaffold.plotting import *
from scaffold.models import MorphologySet, PlacementSet
from scaffold.output import MorphologyRepository, MorphologyCache
from scaffold.reporting import set_verbosity

# Importing compute results
from compute_results import (
    region_name,
    VOXEL_SIZE,
    mask,
    dens_neuron,
    dens_inh,
    layers_per_cell,
    mask_of_stellate,
    orientations,
    fac_Lugaro,
    reference_UBC_density,
)

docs = """
# Per Silvia
import pickle

# obj0, obj1, obj2 are created here...

# Saving the objects:
with open('data_mapping.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([region_name, VOXEL_SIZE, mask, dens_neuron, dens_inh, layers_per_cell, mask_of_stellate, orientations, fac_Lugaro, reference_UBC_density], f)

# Getting back the objects:
# with open('data_mapping.pkl') as f:  # Python 3: open(..., 'rb')
#     region_name, VOXEL_SIZE, mask, dens_neuron, dens_inh, layers_per_cell, mask_of_stellate, orientations, fac_Lugaro, reference_UBC_density = pickle.load(f)
"""
# Importing other packages
# Extract PC layer voxels - indexes
import numpy as np
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname("__file__"), ".."))
from time import time
import plotly.graph_objects as go
import random
import h5py
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import asearch
import scipy
import cProfile


PLOTTING_FOR_CHECK = False
PC_SHUFFLED = False

cell_properties = {
    "stellate_cell": {"radius": 4.0, "id": 5, "type": "inh"},
    "basket_cell": {"radius": 6.0, "id": 6, "type": "inh"},
    "granule_cell": {"radius": 2.5, "id": 2, "type": "exc"},
    "golgi_cell": {"radius": 8.0, "id": 1, "type": "inh"},
    "glomerulus": {"radius": 1.5, "density": 3e-4},
    "purkinje_cell": {"radius": 7.5},
}
cell_per_layer = {
    "Stellate layer": ["stellate_cell"],
    "Basket layer": ["basket_cell"],
    "granular layer": ["granule_cell", "golgi_cell", "glomerulus"],
}

# cell_per_layer = {'Stellate layer': ['stellate_cell']}
# cell_per_layer = {'Stellate layer': ['stellate'],'Basket layer': ['basket']}
# cell_per_layer = {'granular layer': ['granule_cell']}


def extract_voxels(indexes, number, dim, threshold):
    """
    It prepares voxels and cell count to be fed to the particle placement.
    It extracts the voxels (with 3D origin and size), the neuron count in each voxel.

    :param indexes: 3D indexes of voxels belonging to a certain layer from the full voxelized brain/cerebellum/cerebellum part
    :param number: number of neurons in the selected voxels  and a possible dimension along which we want to cut to
    place only in a sub-part of a certain region

    """

    voxel_sizes = np.array([25.0, 25.0, 25.0])
    voxels = []
    neu_in_voxels = []

    indexes = np.transpose(indexes)
    # print("transpose size",indexes.shape)
    indexes = np.asarray(indexes * int(voxel_sizes[0]))

    for ind in range(len(indexes)):
        if indexes[ind][dim] < threshold:
            voxels.append([np.array(indexes[ind]), voxel_sizes])
            neu_in_voxels.append(number[ind])

    return voxels, neu_in_voxels


# Initialize Scaffold
def relative_folder(path):
    return os.path.join(os.path.dirname("__file__"), path)


morpho_file = relative_folder("morphologies.hdf5")

config = JSONConfig(file="mouse_lingula_neuron.json")
config.connection_types["parallel_fiber_to_golgi"].transformation.quivers = orientations
config.connection_types["parallel_fiber_to_purkinje"].transformation.quivers = orientations
config.connection_types["parallel_fiber_to_basket"].transformation.quivers = orientations
config.connection_types["parallel_fiber_to_stellate"].transformation.quivers = orientations
set_verbosity(3)
scaffold = Scaffold(config)
scaffold.morphology_repository = MorphologyRepository(morpho_file)
morphology_cache = MorphologyCache(scaffold.morphology_repository)
morphology_cache.rotate_all_morphologies(30, 30)

num_to_place = {}
positions_all = {}
rotations_all = {}

####################################### Placing Purkinje cells ###########################################################################
dist_x = 130.0
dist_z = 3.5
radius_pc = 7.5

# Extract parasagittal sections
# Number of voxels corresponding to the required x distance
vox_dist_x = math.floor(dist_x / VOXEL_SIZE)  # To be sure of not having overlapping
vox_dist_z = math.ceil(dist_z / VOXEL_SIZE)

index_pc_voxels = np.nonzero(mask["Purkinje layer"])
# print(index_pc_voxels)
pc_matrix = np.zeros(
    (
        max(index_pc_voxels[0]) + 1,
        max(index_pc_voxels[1]) + 1,
        max(index_pc_voxels[2]) + 1,
    )
)
pc_matrix[index_pc_voxels] = 1
# print(pc_matrix)

# PC sections are obtained cutting along the z axis
dim = pc_matrix.shape  # 448 x 232 x 372
nb_sections = len(index_pc_voxels[2])
cut_dir = 2
sections = []
started_placing = False
slice_placing_num = 0
for sd in range(dim[cut_dir]):  # FOR EACH SLICE

    # Matrix cutted in the parasagittal direction
    pc_slice = np.take(pc_matrix, sd - 1, axis=cut_dir)
    # Tuple of arrays with indices of PC voxels/pixels in the current slice, ordered
    # in ascending order of first dimension
    pc_slice_indexes = np.nonzero(pc_slice)
    pc_voxels_in_slice = np.nonzero(pc_slice)
    # print("non zero ",(np.nonzero(pc_slice)))
    # If there are more than 1 voxel labelled as PC in the current slice, we start placing PCs
    print("*******SLICE NUMBER : ", sd, "-- num PC indexes: ", len(pc_slice_indexes[0]))

    past_node = []
    past_end_node = []
    pc_voxels_in_slice_mod = []
    count_over_dist_two = 0

    if len(pc_slice_indexes[0]) > 1:
        first_pc_in_slice = 1  # è la prima pc che posiziono nella slice nella slice

        # PLACING THE FIRST PC (at all)
        if not started_placing:
            # print("prima pc da posizionare nella slice num ", sd )
            start_x = np.min(pc_slice_indexes[0])
            # print("pc_slice_indexes[0] = ", pc_slice_indexes[0], "min value : ", start_x )

            # Take a random voxel from the ones at x=start_x
            start_y = random.choice(np.nonzero(pc_slice[start_x, :]))[0]
            # print ("np.nonzero(pc_slice[start_x,:]) = ", np.nonzero(pc_slice[start_x,:]))
            # print ("random.choice(np.nonzero(pc_slice[start_x,:]))", random.choice(np.nonzero(pc_slice[start_x,:])) )
            # print("starting voxel ",start_x, " ", start_y)
            PC_placed = np.array([[start_x, start_y, sd]])
            # print("prima PC placed in  ", start_x, " ", start_y )
            started_placing = True
            slice_placing_num += 1
            # print ("pc_slice_indexes : ", pc_slice_indexes)
            # Remove placed voxel from tuple of indices
            pc_slice_indexes = (
                np.delete(pc_slice_indexes[0], np.argmin(pc_slice_indexes[0])),
                np.delete(pc_slice_indexes[1], np.argmin(pc_slice_indexes[0])),
            )
            # print ("pc_slice_indexes dopo il delete : ", pc_slice_indexes)
            pc_voxels_in_slice = pc_slice_indexes

        else:  # PLACING THE FIRST PC IN THE FOLLOWING SLICES
            # Take the closest voxel to the first placed in the previous slice
            # print("pc_slice_indexes : ",pc_slice_indexes , " start_x = ", start_x, "start_y", start_y)
            x_dist = np.subtract(pc_slice_indexes[0], start_x)
            y_dist = np.subtract(pc_slice_indexes[1], start_y)
            # print("x_dist = ", x_dist, "y_dist", y_dist )
            # print(np.power(x_dist,2), " ", np.power(y_dist,2))
            distances = np.sqrt(np.add(np.power(x_dist, 2), np.power(y_dist, 2)))
            # print("distances", distances )
            # print(distances, " min ", np.min(distances), " argmin ", np.argmin(distances))
            if PC_SHUFFLED:
                # Add random variability to closest point (the one in the first element of sorted_distance_indexes)
                sorted_distance_indexes = np.argsort(distances)
                # print("sorted_distance_indexes = ", sorted_distance_indexes)
                probabilities = np.linspace(1.0, 0.0, num=len(sorted_distance_indexes))
                probabilities = probabilities / probabilities.sum()
                closest_voxel_index = sorted_distance_indexes[
                    np.random.choice(sorted_distance_indexes, 1, p=probabilities)
                ]
                start_y = pc_slice_indexes[1][closest_voxel_index[0]]
                start_x = pc_slice_indexes[0][closest_voxel_index[0]]
                pc_voxels_in_slice = (
                    np.delete(pc_voxels_in_slice[0], closest_voxel_index[0]),
                    np.delete(pc_voxels_in_slice[1], closest_voxel_index[0]),
                )
            else:  # Take closest PC
                start_x = pc_slice_indexes[0][np.argmin(distances)]
                start_y = pc_slice_indexes[1][np.argmin(distances)]
                pc_voxels_in_slice = (
                    np.delete(pc_voxels_in_slice[0], np.argmin(distances)),
                    np.delete(pc_voxels_in_slice[1], np.argmin(distances)),
                )

            PC_placed = np.concatenate((PC_placed, [[start_x, start_y, sd]]), axis=0)
            print("--> PC placed in  ", start_x, " ", start_y)
            first_pc = 1

        if len(pc_slice_indexes[0]) > 5:  # PLACING ADDITIONAL PC IN SLICE
            # If there are enough voxels in the current slice, we place more than 1 PC array
            # print("Searching additional PC array value...")
            first_pc = 1
            start_x_current_slice = start_x
            start_y_current_slice = start_y
            result = []
            cost = []

            current_min_dist = math.sqrt(
                (start_x_current_slice - pc_voxels_in_slice[0][0]) ** 2
                + (start_y_current_slice - pc_voxels_in_slice[1][0]) ** 2
            )
            ind_near_voxel = 0

            current_node_x = start_x_current_slice
            current_node_y = start_y_current_slice

            tot_voxels = len(pc_voxels_in_slice[0])

            while len(pc_voxels_in_slice[0]) > 0:

                # print("rimanenti pc voxels : ", len(pc_voxels_in_slice[0]))
                # print ("voxels rimanenti sono: ", pc_voxels_in_slice)

                if first_pc_in_slice == 1:  # DEFINING THE FIRST END-NODE IN SLICE
                    empty_list = 0
                    fine = 0
                    dist_two = 0

                    while fine == 0:
                        # searcing voxels of the PCLayer around the current-node
                        adjacent_voxels = [
                            (0, -1),
                            (0, 1),
                            (-1, 0),
                            (1, 0),
                            (-1, -1),
                            (-1, 1),
                            (1, -1),
                            (1, 1),
                        ]
                        if (
                            empty_list == 1
                        ):  # se non ho child che vanno bene a distanza 1, li cerco a distanza 2
                            print("cerco child del primo voxel, a distanza 2 ")
                            adjacent_voxels = [
                                (0, -2),
                                (0, 2),
                                (-2, 0),
                                (2, 0),
                                (-2, -2),
                                (-2, 2),
                                (2, -2),
                                (2, 2),
                                (2, 1),
                                (1, 2),
                                (-1, 2),
                                (-2, 1),
                                (-2, -1),
                                (-1, -2),
                                (1, -2),
                                (2, -1),
                            ]  # Adjacent squares
                            dist_two = 1

                        children = []
                        for new_position in adjacent_voxels:
                            node_position = (
                                current_node_x + new_position[0],
                                current_node_y + new_position[1],
                            )

                            if (
                                pc_slice[node_position[0]][node_position[1]] == 0
                            ):  # il voxel deve appartenere al PC layer
                                continue

                            children.append(node_position)

                        if len(children) == 0:
                            empty_list = 1  # non ho trovato child a distanza 1
                            if dist_two == 0:
                                print("non ho trovato child a distanza 1 adiacenti al primo voxel")
                            else:  # searching closest voxels of PCLayer
                                # print("neanche a distanza 2 non ho trovato children per il primo voxel")
                                # print (" il nodo è ", current_node_x, current_node_y, "della slice", sd)
                                # print ("cerco nei voxels a distanza > di 2 ")
                                min_dist = 1000
                                for ind in range(len(pc_voxels_in_slice[0])):
                                    distance_new = math.sqrt(
                                        (current_node_x - pc_voxels_in_slice[0][ind]) ** 2
                                        + (current_node_y - pc_voxels_in_slice[1][ind]) ** 2
                                    )

                                    if distance_new < min_dist:
                                        min_dist = distance_new
                                        ind_near = (
                                            ind  # ind è l'indice con cui seleziono il nodo di end
                                        )
                                        next_child_x = pc_voxels_in_slice[0][ind_near]
                                        next_child_y = pc_voxels_in_slice[1][ind_near]
                                # print("posizione del prossimo nuovo child =", next_child_x, next_child_y)
                                node_next_child = (next_child_x, next_child_y)
                                children.append(node_next_child)

                                fine = 1
                        else:
                            fine = 1

                    # for child in children:
                    # print ("- child ",child)

                    min_index_y = np.min(pc_slice_indexes[1])
                    max_index_y = np.max(pc_slice_indexes[1])
                    # print("max index y =", max_index_y, "| min index y =", min_index_y)
                    diff_all = []
                    if (max_index_y - current_node_y) > (
                        current_node_y - min_index_y
                    ):  # voglio andare verso max_index_y

                        for child in children:
                            # print("calcolo ", max_index_x, "-", child[0] )
                            diff = max_index_y - child[1]
                            # print ("diff = ", diff )
                            diff_all.append(diff)
                            min_diff = np.amin(diff_all)  # valore minima distanza
                            index_diff_min = diff_all.index(min_diff)

                    else:  # voglio andare verso min_index_y
                        for child in children:
                            # print("calcolo ",  child[0], "-", min_index_x, )
                            diff = child[1] - min_index_y
                            # print ("diff = ", diff )
                            diff_all.append(diff)
                            min_diff = np.amin(diff_all)  # valore minima distanza
                            index_diff_min = diff_all.index(min_diff)

                    # print(" minima diff  = ", min_diff, "che è associato alla posizione: ", index_diff_min)
                    # print("nodo = ", children[index_diff_min])
                    #    distance_start_max_index = math.sqrt(((child.position[0] - current_node.position[0]) ** 2) + ((child.position[1] - current_node.position[1]) ** 2))
                    # print("")
                    current_node_x = children[index_diff_min][0]
                    current_node_y = children[index_diff_min][1]
                    first_pc_in_slice = 0
                    flag_for_delete = 1
                else:  # DEFINING THE FOLLOWING END-NODE
                    flag_for_delete = 0
                    past_end_node_temp = (current_node_x, current_node_y)
                    # print("past end node ", past_end_node)
                    # print("p-ast current node: ", current_node_x, current_node_y)
                    # past_current_node = (current_node_x,current_node_y)
                    current_min_dist = 1000
                    for ind in range(len(pc_voxels_in_slice[0])):
                        distance_new = math.sqrt(
                            (current_node_x - pc_voxels_in_slice[0][ind]) ** 2
                            + (current_node_y - pc_voxels_in_slice[1][ind]) ** 2
                        )
                        # print ("distance bw ", current_node_x, current_node_y, "e ", pc_voxels_in_slice[0][ind], pc_voxels_in_slice[1][ind], "=", distance_new)
                        if distance_new <= current_min_dist:
                            current_min_dist = distance_new
                            ind_near_voxel = ind  # ind è l'indice con cui seleziono il nodo di end
                            next_current_node_x = pc_voxels_in_slice[0][ind_near_voxel]
                            next_current_node_y = pc_voxels_in_slice[1][ind_near_voxel]

                    current_node_x = next_current_node_x
                    current_node_y = next_current_node_y
                    # print (" ***voxel più vicino: ", current_node_x, current_node_y, "          dist=", current_min_dist )

                (path, path_length, count_over_dist_two, pc_voxels_in_slice_mod,) = asearch.astar(
                    pc_slice,
                    (start_x_current_slice, start_y_current_slice),
                    (current_node_x, current_node_y),
                    past_end_node,
                    pc_voxels_in_slice,
                    count_over_dist_two,
                    pc_voxels_in_slice_mod,
                    start_x,
                    start_y,
                )
                # print("path: ", path, "path_lenght = ", path_length )

                if flag_for_delete == 0:
                    # print("elimino dalla lista ",pc_voxels_in_slice[0][ind_near_voxel], pc_voxels_in_slice[1][ind_near_voxel]  )
                    pc_voxels_in_slice = (
                        np.delete(pc_voxels_in_slice[0], ind_near_voxel),
                        np.delete(pc_voxels_in_slice[1], ind_near_voxel),
                    )

                # mi assicuro che di essere a distanza sufficiente dalla prima pc posizionata nella slice
                dist_first_pc = math.sqrt(
                    (current_node_x - start_x) ** 2 + (current_node_y - start_y) ** 2
                )
                # print ("distance bw", current_node_x , current_node_y , "e", start_x, start_y, "=", dist_first_pc  )
                # if ((dist_first_pc<vox_dist_x) and (first_pc == 0)):   #first_pc = 0, --> non sto posizionando la seconda pc
                # print("non poiziono la pc in ", current_node_x, current_node_y, "perchè troppo vicina alla pc in ", start_x, start_y)

                place = 0
                # se devo posizionare la seconda pc, non tengo conto della dist in linea d'aria dall prima pc
                if (first_pc == 1) or (first_pc == 0 and dist_first_pc > vox_dist_x):
                    place = 1

                if (path_length > vox_dist_x) and (place == 1):
                    print("--> Placing PC in ", current_node_x, " ", current_node_y)
                    first_pc = 0
                    # print("past current node = ", past_current_node)
                    # previous_node_list.append(past_current_node)
                    # print(" previous node list = ", previous_node_list)
                    past_end_node = past_end_node_temp

                    PC_placed = np.concatenate(
                        (PC_placed, [[current_node_x, current_node_y, sd]]), axis=0
                    )
                    start_x_current_slice = current_node_x
                    start_y_current_slice = current_node_y
                    count_over_dist_two = 0
    # else:
    # print("no enough pc voxels in slice ", sd )


print("Placed ", len(PC_placed), " Purkinje cells")


# Position PC at center of voxel with variability within VOXEL
positions_all["Purkinje layer"] = (
    PC_placed * VOXEL_SIZE
    + VOXEL_SIZE / 2
    + scipy.stats.truncnorm.rvs(-VOXEL_SIZE / 2, VOXEL_SIZE / 2)
)
radii_pc = np.ones([len(PC_placed), 1]) * cell_properties["purkinje_cell"]["radius"]
print(radii_pc.shape, " ", positions_all["Purkinje layer"].shape)
positions_all["Purkinje layer"] = np.hstack((radii_pc, positions_all["Purkinje layer"]))
print(len(positions_all["Purkinje layer"]))

num_Lugaro = np.round(len(PC_placed) / fac_Lugaro)

####################################### Placing the other neurons ###########################################################################
pruned = {}
pruned_per_type_in_layer = {}
t_start = time()

for layer in cell_per_layer.keys():
    print(" *******    PLACING NEURONS IN THE ", layer.upper())

    neu_in_voxels = {}

    # Indexes of voxels belonging to layer (3-D tuple of arrays, each one for each of the 3 dimensions)
    index_cell = np.nonzero(mask[layer])

    if PLOTTING_FOR_CHECK:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(index_cell[0] * 25, index_cell[1] * 25, index_cell[2] * 25)
        plt.show()

    # Extracting map of neuron count in the current layer
    number_neu = dens_neuron[mask[layer]]
    number_inh = dens_inh[mask[layer]]
    number_exc = number_neu - number_inh
    # print(" number_neu.shape =", number_neu.shape )
    # print(" number_inh.shape =", number_inh.shape )
    # print(" number_exc.shape =", number_exc.shape )

    # print of the total number of neurons, inhibitory and excitatory for check
    number_neu_tot = []
    number_inh_tot = []
    number_exc_tot = []
    number_neu_tot = np.round(np.sum(number_neu))
    number_inh_tot = np.round(np.sum(number_inh))
    number_exc_tot = np.round(np.sum(number_exc))
    print(
        "number neu_tot = ",
        number_neu_tot,
        "// inh= ",
        number_inh_tot,
        " //exc = ",
        number_exc_tot,
    )

    # For each cell type in the layer, extract voxels belonging to that layer and cell counts for that cell type.
    # voxels_cell is computed 2 times, redundant..to improve!
    for cell in cell_per_layer[layer]:
        print("-cell type : ", cell)
        if "type" in cell_properties[cell].keys():  # It is not a glomerulus
            if cell_properties[cell]["type"] == "exc":
                print("   extracting voxels ext type...")
                voxels_cell, neu_in_voxels[cell] = extract_voxels(index_cell, number_exc, 2, 12000)
            else:
                print("   extracting voxels inh type...")
                voxels_cell, neu_in_voxels[cell] = extract_voxels(index_cell, number_inh, 2, 12000)
        else:
            # voxels_cell are the same as just computed for the previous cell of the layer
            print("    extracting for glomerulus....")
            neu_in_voxels[cell] = [25 * 25 * 25 * cell_properties[cell]["density"]] * len(
                neu_in_voxels["golgi_cell"]
            )

        if PLOTTING_FOR_CHECK:
            fig = plt.figure()
            ax = Axes3D(fig)
            for vox in range(len(voxels_cell)):
                ax.scatter(
                    voxels_cell[vox][0][0],
                    voxels_cell[vox][0][1],
                    voxels_cell[vox][0][2],
                )
            plt.show()

    voxel_num = 50

    tot_pruned = []
    pruned_per_type = {}
    per_pruned = []

    print("num vox group ", int(len(voxels_cell) / voxel_num))
    for v in range(int(len(voxels_cell) / voxel_num)):
        print("---- subvolume number : ", v)
        voxels_to_place = voxels_cell[v : v + voxel_num]
        particles_cell = []

        for cell in cell_per_layer[layer]:
            num_to_place[cell] = math.ceil(
                sum(neu_in_voxels[cell][v * voxel_num : (v + 1) * voxel_num])
            )
            print("* cell type = ", cell)
            print("    num to place: ", num_to_place[cell])
            particles_cell.append(
                {
                    "name": cell,
                    "voxels": list(range(v * voxel_num, (v + 1) * voxel_num - 1)),
                    "radius": cell_properties[cell]["radius"],
                    "count": num_to_place[cell],
                }
            )

        t = time()
        system = ParticleSystem()
        system.fill(voxels_cell, particles_cell)
        colliding = system.find_colliding_particles()
        # print("{} particles with packing factor {} have {} collisions.".format(
        #   len(system.particles),
        #   system.get_packing_factor(),
        #   len(system.colliding_particles)
        # ))
        # print("Calculated in {} seconds.".format(time() - t))
        # fig = plot_detailed_system(system)
        t = time()
        system.solve_collisions()
        # fig = plot_detailed_system(system)
        # print("Collisions fixed in {} seconds.".format(time() - t))
        tot, per_type = system.prune(colliding)
        tot_pruned.append(tot)
        print("per type: ", per_type)

        for t in per_type.keys():
            if t not in pruned_per_type.keys():
                pruned_per_type[t] = []
            pruned_per_type[t].append(per_type[t])

        if (
            len(colliding) > 0
        ):  # Due to collisions, some particles have been replaced and could have been pruned
            print(
                "{} particles pruned, {}% of the total count, {}% of the colliding particles.".format(
                    tot_pruned[v],
                    int((tot_pruned[v] / num_to_place[cell]) * 100),
                    int((tot_pruned[v] / len(colliding)) * 100),
                )
            )
        per_pruned.append(int((tot_pruned[v] / (num_to_place[cell])) * 100))
        # fig = plot_detailed_system(system)
        particle_positions = np.array([p.position for p in system.particles])
        particle_radii = np.array([[p.radius for p in system.particles]])
        particle_radii = particle_radii.reshape(np.size(particle_radii), 1)
        positions = np.hstack((particle_radii, particle_positions))
        # positions = np.concatenate((cell_properties[cell]['id']*np.ones((len(system.particles),1)), particle_positions), axis=1)
        # Saving the positions of placed cells together with their radius in a numpy array with the first column as the radius and the last 3 columns as the 3D position
        if v == 0:
            positions_all[layer] = positions
        else:
            positions_all[layer] = np.concatenate((positions_all[layer], positions), axis=0)

    pruned_per_type_in_layer[layer] = {}
    for t in pruned_per_type.keys():
        pruned_per_type_in_layer[layer][t] = sum(pruned_per_type[t])

    pruned[cell] = sum(tot_pruned)
    print("Pruned ", pruned[cell], cell, " cells")
    print("Pruned per type in layers ", pruned_per_type_in_layer)

positions_golgi = positions
t_end = time()
tot_time = t_end - t_start
print("Pruned cells: ", pruned_per_type_in_layer)
print("Total positioning time: ", tot_time)


cell_per_layer["Purkinje layer"] = ["purkinje_cell"]

discretized_phi, discretized_theta = morphology_cache._discretize_orientations([30, 30])
discretized_orientations = np.array(
    [np.cos(discretized_phi), np.sin(discretized_phi), np.sin(discretized_theta)]
)

for layer in cell_per_layer.keys():
    print("Associating orientations in ", layer)
    rotations_all[layer] = np.empty([0, 2])
    voxel_ind = positions_all[layer][:, 1:] / VOXEL_SIZE
    voxel_ind = voxel_ind.astype(int)
    print("voxel indexes shape ", voxel_ind.shape)
    print("orient shape ", orientations.shape)
    # Associate closest orientation to each cell p in the current layer
    for p in range(len(positions_all[layer])):
        # print("positions[p] ", positions[p])

        # print("max voxel ind ", np.amax(voxel_ind,axis=0), " ", np.amax(voxel_ind[1]))
        # print("orientation vector values ",orientations[:,voxel_ind[0],voxel_ind[1],voxel_ind[2]])
        # print("or vec ",orientations[:,voxel_ind[p,0],voxel_ind[p,1],voxel_ind[p,2]])
        current_orientation = np.array(
            [
                orientations[:, voxel_ind[p, 0], voxel_ind[p, 1], voxel_ind[p, 2]],
            ]
            * len(discretized_orientations[0])
        )
        # print("power",np.power((current_orientation.T - discretized_orientations),2))
        # print("sum ",np.sum(np.power((current_orientation.T - discretized_orientations),2),axis=0))
        orientation_distance = np.sqrt(
            np.sum(np.power((current_orientation.T - discretized_orientations), 2), axis=0)
        )
        # print("or distance ",orientation_distance)
        closest_index = np.argmin(orientation_distance)
        # closest_discrete_orientation = discretized_orientations[:,closest_index]
        closest_discrete_orientation = np.array(
            [
                [discretized_phi[closest_index], discretized_theta[closest_index]],
            ]
        )
        # print(closest_discrete_orientation.shape)
        rotations_all[layer] = np.concatenate(
            (rotations_all[layer], closest_discrete_orientation), axis=0
        )

    print("shape rotations ", rotations_all[layer].shape)
    # Build scaffold
    for cell in cell_per_layer[layer]:
        cell_type = scaffold.get_cell_type(cell)
        if cell == "granule_cell":
            index_granule = np.where(positions_all[layer][:, 0] == cell_properties[cell]["radius"])
            positions_granule = positions_all[layer][index_granule, 1:]
            positions_granule = np.squeeze(positions_granule)
            rotations_granule = rotations_all[layer][index_granule, :]
            rotations_granule = np.squeeze(rotations_granule)
            print("pos shape ", positions_granule.shape)
            print("rot shape ", rotations_granule.shape)
            scaffold.place_cells(
                cell_type,
                cell_type.placement.layer_instance,
                positions_granule,
                rotations_granule,
            )
        elif cell == "golgi_cell":
            index_golgi = np.where(positions_all[layer][:, 0] == cell_properties[cell]["radius"])
            positions_golgi = positions_all[layer][index_golgi, 1:]
            positions_golgi = np.squeeze(positions_golgi)
            rotations_golgi = rotations_all[layer][index_golgi, :]
            rotations_golgi = np.squeeze(rotations_golgi)
            scaffold.place_cells(
                cell_type,
                cell_type.placement.layer_instance,
                positions_golgi,
                rotations_golgi,
            )
        elif cell == "glomerulus":
            index_glom = np.where(positions_all[layer][:, 0] == cell_properties[cell]["radius"])
            positions_glom = positions_all[layer][index_glom, 1:]
            positions_glom = np.squeeze(positions_glom)
            rotations_glom = rotations_all[layer][index_glom, :]
            rotations_glom = np.squeeze(rotations_glom)
            scaffold.place_cells(
                cell_type,
                cell_type.placement.layer_instance,
                positions_glom,
                rotations_glom,
            )
        else:
            scaffold.place_cells(
                cell_type,
                cell_type.placement.layer_instance,
                positions_all[layer][:, 1:],
                rotations_all[layer],
            )

# Place entities
mossy_fibers_type = scaffold.get_cell_type("mossy_fibers")
scaffold.place_cell_type(mossy_fibers_type)


scaffold.compile_output()
print("Connecting neurons...")
# cProfile.run("scaffold.connect_cell_types()","profile_results.stat")
scaffold.connect_cell_types()
scaffold.compile_output()
