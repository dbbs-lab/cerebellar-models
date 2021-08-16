# flake8: noqa
import numpy as np
import json


def return_ids_containing_str_list(str_list):
    id_list = []
    for kk in id_to_region_dictionary_ALLNAME.keys():
        region_is_in = True
        for str1 in str_list:
            if (id_to_region_dictionary_ALLNAME[kk].lower()).find(
                str1.lower()
            ) < 0:  # if any of the regions is not there, do not take
                region_is_in = False
        if region_is_in:
            id_list.append(kk)
    return id_list


def hex_to_rgb(value):
    value = value.lstrip("#")
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def search_children(object_, numiter, lastname_ALL="", lastname="", pos0=0.5, sc=1.0, darken=True):
    global iterTMP, region_keys
    # ~ print " "*numiter, object_["name"], ":", object_["id"]
    # ~ id_to_region_dictionary[ object_["id"] ] = lastname+"/"+object_["name"]
    # ~ id_to_region_dictionary_ALLNAME[ object_["id"] ] = object_["name"]
    # ~ region_dictionary_to_id[ object_["name"] ] = object_["id"]
    # ~ region_dictionary_to_id_ALLNAME[ lastname+"/"+object_["name"] ] = object_["id"]

    regions_ALLNAME_list.append(lastname_ALL + "|" + object_["name"])
    name2allname[object_["name"]] = lastname_ALL + "|" + object_["name"]
    allname2name[lastname_ALL + "|" + object_["name"]] = object_["name"]
    id_to_region_dictionary[object_["id"]] = object_["name"]
    id_to_region_dictionary_ALLNAME[object_["id"]] = lastname_ALL + "|" + object_["name"]
    region_dictionary_to_id[object_["name"]] = object_["id"]
    region_dictionary_to_id_ALLNAME[lastname_ALL + "|" + object_["name"]] = object_["id"]
    region_dictionary_to_id_ALLNAME_parent[lastname_ALL + "|" + object_["name"]] = lastname_ALL
    region_dictionary_to_id_parent[object_["name"]] = lastname

    # coloring process
    clrTMP = np.float32(np.array(list(hex_to_rgb(object_["color_hex_triplet"]))))
    if np.sum(clrTMP) > 255.0 * 3.0 * 0.75 and darken:
        clrTMP *= 255.0 * 3.0 * 0.75 / np.sum(clrTMP)
    region_to_color[lastname_ALL + "|" + object_["name"]] = list(clrTMP)
    id_to_color[object_["id"]] = list(clrTMP)

    regions_pos[object_["name"]] = pos0
    allnameOrder[lastname_ALL + "|" + object_["name"]] = iterTMP
    iterTMP += 1
    region_keys.append(object_["name"])
    try:
        is_leaf[lastname_ALL + "|" + object_["name"]] = 1
        # ~ region_dictionary_to_id_ALLNAME_child[  lastname_ALL+"|"+object_["name"] ] = children
        # ~ id_children[object_["id"]] = object_["children"]
        for ichildren, children in enumerate(object_["children"]):
            # ~ search_children(children, numiter+1, lastname_ALL+"/"+object_["name"])
            numchilds = len(object_["children"]) - 1
            if numchilds == 0:
                posTMP = pos0 + 0.0
            else:
                posTMP = pos0 + sc * (float(ichildren) / float(numchilds) - 0.5)
            search_children(
                children,
                numiter + 1,
                lastname_ALL + "|" + object_["name"],
                object_["name"],
                pos0=posTMP,
                sc=sc / float(len(object_["children"])),
            )
            is_leaf[lastname_ALL + "|" + object_["name"]] = 0
    except:
        print("No children of object")


DATA_PATH = "../Flocculus3.0_Lingula/"
region_name = "Flocculus"


id_to_region_dictionary = {}
id_to_region_dictionary_ALLNAME = {}  # id to complete name
region_dictionary_to_id = {}  # name to id
region_dictionary_to_id_ALLNAME = {}  # complete name to id
region_dictionary_to_id_ALLNAME_parent = {}  # complete name to complete name parent
region_dictionary_to_id_parent = {}  # name to complete name parent
allname2name = {}  # complete name to name
name2allname = {}  # name to complete name
allnameOrder = {}
iterTMP = 0
region_keys = []  # list of regions names
regions_ALLNAME_list = []  # list of complete regions names
is_leaf = {}  # full name to int (! if is leaf, else 0)
regions_pos = {}  # name to position in ordered array (according to depth in tree)
id_to_color = {}  # region id to color in RGB
region_to_color = {}  # complete name to color in RGB


jsontextfile = open(DATA_PATH + "data/brain_regions.json", "r")
jsoncontent = json.loads(jsontextfile.read())

search_children(jsoncontent["msg"][0], 0)


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


with open("tests/fixtures/id_to_region_dictionary_ALLNAME.json", "w+") as f:
    json.dump(id_to_region_dictionary_ALLNAME, f)


with open("tests/fixtures/id_region.json", "w+") as f:
    json.dump(id_region, f)


with open("tests/fixtures/id_current_region.json", "w+") as f:
    json.dump(id_current_region, f)
