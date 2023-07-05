
import numpy as np
from bsb.morphologies import Morphology

def relabel_func(morpho : Morphology) -> Morphology:
    cell = morpho
    print(cell)
    #cell = Morphology.from_swc("morphologies/GranuleCell.swc")
    cell.center()
    branches = [b for b in cell.branches if np.any(b.tags == 16)]
    for b in branches:
        b.label(["axon", "axon_hillock"])
    branches = [b for b in cell.branches if np.any(b.tags == 17)]
    for b in branches:
        b.label(["axon", "axon_initial_segment"])
    branches = [b for b in cell.branches if np.any(b.tags == 18)]
    for b in branches:
        b.label(["axon", "ascending_axon"])
    branches = [b for b in cell.branches if np.any(b.tags == 19)]
    for b in branches:
        b.label(["axon", "parallel_fiber"])
    return cell
