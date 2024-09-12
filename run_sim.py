import nest
from bsb import from_storage, options

options.verbosity = 4
network = from_storage("mouse_cereb_dcn_nest.hdf5")

sim_results = network.run_simulation("basal_activity")
