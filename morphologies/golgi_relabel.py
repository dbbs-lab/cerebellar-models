import numpy as np
from bsb.morphologies import Morphology

def relabel_func(morpho : Morphology) -> Morphology:
    cell = morpho
    cell.center()
    branches = [b for b in cell.branches if np.all(b.tags == 16)]
    for b in branches:
        b.label(["dendrites", "basal_dendrites"])
    branches = [b for b in cell.branches if np.all(b.tags == 17)]
    for b in branches:
        b.label(["dendrites", "apical_dendrites"])
    branches = [b for b in cell.branches if np.all(b.tags == 18)]
    for b in branches:
        b.label(["axon", "axonal_initial_segment"])
    print(cell.flatten_labels())
    return cell