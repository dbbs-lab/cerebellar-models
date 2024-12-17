import mpi4py
import nest
from bsb import ConfigurationError, SimulationData, from_storage
from bsb_nest import NestAdapter
from bsb_nest.adapter import NestResult

adapter = NestAdapter()
forward_model = from_storage("mouse_cereb_microzones.hdf5")
inverse_model = from_storage("mouse_cereb_microzones.hdf5")
