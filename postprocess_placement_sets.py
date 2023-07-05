from bsb.core import from_storage
import numpy as np
from bsb.storage import Storage

network = from_storage("MouseCerebellum.hdf5")

for ps in network.get_placement_sets():
    print("numb of ", ps.tag, ": ", len(ps))
