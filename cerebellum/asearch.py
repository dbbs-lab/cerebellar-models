# flake8: noqa
import math
import numpy as np


class Node:
    """A node class for A* Pathfinding"""

    def __init__(self, parent=None, position=None):
        self.parent = parent
        self.position = position

        self.distance_from_parent = 0
        self.g = 0
        self.h = 0
        self.f = 0

    def __eq__(self, other):
        return self.position == other.position


def astar(
    maze,
    start,
    end,
    past_node,
    pc_voxels_in_slice,
    count_over_dist_two,
    pc_voxels_in_slice_mod,
    start_x,
    start_y,
):
    """Returns a list of tuples as a path from the given start to the given end in the given maze"""

    # print("len(maze) = ", len(maze) , "--(len(maze)-1)= ", (len(maze)-1) )
    # print("len(maze[len(maze)-1])= ",  len(maze[len(maze)-1]) , "-- (len(maze[len(maze)-1])-1) = ",( len(maze[len(maze)-1]) -1)   )
    distance_start_end = 0
    # Create start and end node
    start_node = Node(None, start)
    start_node.g = start_node.h = start_node.f = 0
    end_node = Node(None, end)
    end_node.g = end_node.h = end_node.f = 0
    # print("_|_|_|_|_|__START NODE: ", start_node.position , "-- END NODE: ", end_node.position )
    # print ("num voxels = ", num_voxels)

    pc_voxels_current = []
    new_end_node = 1

    # Initialize both open and closed list
    open_list = []
    closed_list = []

    open_list.append(start_node)  # Add the start node

    # Loop until you find the end
    while len(open_list) > 0:

        # i = 0
        # print("             open list = ", end =',' )
        # for open_node in open_list:
        #    print(open_list[i].position, end =',')
        #    i = i+1
        # print(" ")

        # j = 0
        # print("           closed list = ", end =',' )
        # for closed_node in closed_list:
        #    print(closed_list[j].position, end =',')
        #    j = j+1
        # print(" ")

        current_node = open_list[0]  # Get the current node
        # print("current_node.f = ", current_node.f)
        current_index = 0
        for index, item in enumerate(open_list):
            if item.f < current_node.f:
                current_node = item
                current_index = index
        # print("NEW CURRENT NODE POSITION", current_node.position[0], current_node.position[1]  )
        # time.sleep(0.5)

        # Pop current off open list, add to closed list
        # open_list.pop(current_index)
        open_list = []  # for every new node, reset of open list
        closed_list.append(current_node)

        # Found the goal
        if current_node == end_node:

            # print ("past_current_node ", past_current_node.position)
            path = []
            length = 0
            current = current_node

            while current is not None:
                path.append(current.position)
                length += current.distance_from_parent
                current = current.parent

                dist_first_pc = math.sqrt(
                    (current_node.position[0] - start_x) ** 2
                    + (current_node.position[1] - start_y) ** 2
                )
                # print("length ",length)
                if (length > 5) and (dist_first_pc > 5):  # PC positioning
                    pc_voxels_in_slice_mod = []

            return (
                path[::-1],
                length,
                count_over_dist_two,
                pc_voxels_in_slice_mod,
            )  # Return reversed path and path lengt

        end_placing = 0
        empty_openlist = 0

        while end_placing == 0:

            next_child_x = []
            next_child_y = []

            adjacent_voxels = [
                (0, -1),
                (0, 1),
                (-1, 0),
                (1, 0),
                (-1, -1),
                (-1, 1),
                (1, -1),
                (1, 1),
            ]  # Adjacent squares
            if (
                empty_openlist == 1
            ):  # if no children can be found at distance 1, searches at distance 2
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

            # Generate children
            children = []
            for new_position in adjacent_voxels:  # Adjacent squares

                # Get node position
                node_position = (
                    current_node.position[0] + new_position[0],
                    current_node.position[1] + new_position[1],
                )

                # Make sure within range
                if (
                    node_position[0] > (len(maze) - 1)
                    or node_position[0] < 0
                    or node_position[1] > (len(maze[len(maze) - 1]) - 1)
                    or node_position[1] < 0
                ):
                    continue

                # Make sure walkable terrain
                if maze[node_position[0]][node_position[1]] == 0:
                    continue

                if node_position == past_node:
                    continue

                # Create new node
                new_node = Node(current_node, node_position)

                # Append
                children.append(new_node)

            new_child = 0
            # Loop through children
            for child in children:

                ok = 0  # flag to check if children is in closed_list or not; if yes, I don't put it in the open_list

                # print ("- child ",child.position)
                # Child is on the closed list
                for closed_child in closed_list:
                    if child == closed_child:
                        # print("            --> child ", child.position, "is in closed list")
                        ok = 1
                        break  # exit from 'for closed_child in closed_list'

                # Create the f, g, and h values
                distance_from_parent = math.sqrt(
                    ((child.position[0] - current_node.position[0]) ** 2)
                    + ((child.position[1] - current_node.position[1]) ** 2)
                )
                # print("distance parent: ",distance_from_parent)
                child.g = current_node.g + 1
                child.distance_from_parent = distance_from_parent
                # Substituted euclidean distance to compute path length from origin node
                # child.g = current_node.g + distance_from_parent
                # print("child distance from parent: ",distance_from_parent)
                # print(child.position, " ", current_node.position)
                child.h = ((child.position[0] - end_node.position[0]) ** 2) + (
                    (child.position[1] - end_node.position[1]) ** 2
                )
                child.f = child.g + child.h

                # Child is already in the open list
                for open_node in open_list:
                    if child == open_node and child.g > open_node.g:
                        continue

                # Add the child to the open list
                if (
                    ok == 0
                ):  # se il child non è nella close list, allora va bene, e lo aggiungo all open_list
                    open_list.append(child)
                    new_child = +1

            if new_child == 0:  # no children found for current node
                # print(" no child at distance 1")
                if empty_openlist == 0:
                    empty_openlist = 1  #  search for children at distance 2
                else:  # if empty_openlist is already 1, it means no children can be found even at distance 2
                    # ...searching in the closest voxel
                    count_over_dist_two = count_over_dist_two + 1

                    if (
                        count_over_dist_two == 1
                    ):  # first time I search over distance 2 for the last positioned PC
                        pc_voxels_in_slice_mod = pc_voxels_in_slice  # remaining voxels in new list

                    # print( " new end node ", new_end_node)
                    if new_end_node == 1:
                        pc_voxels_current = pc_voxels_in_slice_mod  # back to previous voxel list, i.e. without the cancelled ones
                        new_end_node = 0

                    current_min_dist = 1000
                    # in pc_voxels in slice  i nodi end vengono eliminati, quindi non ruscirei più
                    # ad avere un nodo child in quella posizione. quindi formo un altra lista dove il nodo end non viene tolto
                    # print (" voxels dove cerco il nuovo child: ", pc_voxels_current )
                    for ind in range(len(pc_voxels_current[0])):
                        distance_new = math.sqrt(
                            (current_node.position[0] - pc_voxels_current[0][ind]) ** 2
                            + (current_node.position[1] - pc_voxels_current[1][ind]) ** 2
                        )
                        # print (" - distance bw ", current_node.position[0], current_node.position[1], "e ", pc_voxels_current[0][ind], pc_voxels_current[1][ind], "=", distance_new)
                        if distance_new < current_min_dist:
                            current_min_dist = distance_new
                            ind_near_voxel = ind  # ind è l'indice con cui seleziono il nodo di end
                            next_child_x = pc_voxels_current[0][ind_near_voxel]
                            next_child_y = pc_voxels_current[1][ind_near_voxel]
                    # print("posizione del prossimo nuovo child =", next_child_x, next_child_y)

                    next_child = Node(current_node, (next_child_x, next_child_y))
                    next_child.distance_from_parent = current_min_dist
                    # print( " distanza nuovo chld = ", next_child.distance_from_parent)
                    next_child.g = current_node.g + 1
                    next_child.h = ((next_child.position[0] - end_node.position[0]) ** 2) + (
                        (next_child.position[1] - end_node.position[1]) ** 2
                    )
                    next_child.f = next_child.g + next_child.h
                    open_list.append(next_child)

                    pc_voxels_current = (
                        np.delete(pc_voxels_current[0], ind_near_voxel),
                        np.delete(pc_voxels_current[1], ind_near_voxel),
                    )

                    end_placing = 1

            else:
                end_placing = 1
