import nest
import mpi4py

from bsb import from_storage, ConfigurationError, SimulationData
from bsb_nest import  NestAdapter
from bsb_nest.adapter import  NestResult

adapter = NestAdapter()
forward_model = from_storage("mouse_cereb_microzones.hdf5")
inverse_model = from_storage("mouse_cereb_microzones.hdf5")