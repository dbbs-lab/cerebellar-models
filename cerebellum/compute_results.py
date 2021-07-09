from JSONread import * 
import nrrd

def get_mask_purkinje():
    ann, h = nrrd.read(DATA_PATH + "data/annotations.nrrd")
    region_names = []
    mask = {}
    for id_ in id_region:
        region_names.append(id_to_region_dictionary[id_])
        region, layer = id_to_region_dictionary[id_].split(", ")
        mask[layer] = ann == id_
    mask[region_name] = ann == id_current_region
    return mask