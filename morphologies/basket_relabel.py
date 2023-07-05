import numpy as np
from bsb.morphologies import Morphology

def relabel_func(morpho : Morphology) -> Morphology:
    cell = morpho
    cell.center()
    branches = [b for b in cell.branches if np.all(b.tags == 16)]
    for b in branches:
        b.label(["axon", "axon_initial_segment"])
    return cell