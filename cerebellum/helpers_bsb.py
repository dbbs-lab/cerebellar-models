# flake8: noqa
from compute_results import get_mask_purkinje
from bsb.config import from_json
from bsb.core import Scaffold, from_hdf5
from helpers import relative_folder, extract_voxels, cell_properties, cell_per_layer

import math
import scipy
import random
import asearch
import bsb.options
import numpy as np

# display
bsb.options.verbosity = 3
PLOTTING_FOR_CHECK = False
PC_SHUFFLED = False

# define
VOXEL_SIZE = 25.0  # um
dist_x = 130.0
dist_z = 3.5
radius_pc = 7.5

config = from_json("../atlas_cerebellum.json")
scaffold = Scaffold(config)

num_to_place = {}
positions_all = {}
rotations_all = {}

####################################### Placing Purkinje cells ###########################################################################
mask = get_mask_purkinje()

# Extract parasagittal sections
# Number of voxels corresponding to the required x distance
vox_dist_x = math.floor(dist_x / VOXEL_SIZE)  # To be sure of not having overlapping
vox_dist_z = math.ceil(dist_z / VOXEL_SIZE)

index_pc_voxels = np.nonzero(mask["Purkinje layer"])
print(index_pc_voxels)
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
