import numpy as np
from bsb.morphologies import Morphology

def relabel_func(morpho : Morphology) -> Morphology:
    cell = morpho
    cell.center()
    branches = [b for b in cell.branches if np.all(b.tags == 16)]
    for b in branches:
        b.label(["axon", "AIS"])
    branches = [b for b in cell.branches if np.all(b.tags == 17)]
    for b in branches:
        b.label(["axon", "AIS_K"])
    branches = [b for b in cell.branches if np.all(b.tags == 18)]
    for b in branches:
        b.label(["axon", "axonmyelin"])
    branches = [b for b in cell.branches if np.all(b.tags == 19)]
    for b in branches:
        b.label(["axon", "nodes"])
    branches = [b for b in cell.branches if np.all(b.tags == 20)]
    for b in branches:
        b.label(["dendrites", "basal_dendrites"])
    branches = [b for b in cell.branches if np.all(b.tags == 21)]
    for b in branches:
        b.label(["dendrites", "pf_targets", "sc_targets"])
    branches = [b for b in cell.branches if np.all(b.tags == 22)]
    for b in branches:
        b.label(["dendrites", "aa_targets", "sc_targets"])
    return cell