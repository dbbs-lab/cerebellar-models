from bsb.storage import Storage
from bsb.morphologies import Morphology


mr = Storage("hdf5", "morphologies.hdf5").morphologies
mr.save("GranuleCell", Morphology.from_swc("morphologies/GranuleCell.swc"))
mr.save("GolgiCell", Morphology.from_swc("morphologies/GolgiCell.swc"))
mr.save("PurkinjeCell", Morphology.from_swc("morphologies/PurkinjeCell.swc"))
mr.save("BasketCell", Morphology.from_swc("morphologies/BasketCell.swc"))
mr.save("StellateCell", Morphology.from_swc("morphologies/StellateCell.swc"))
