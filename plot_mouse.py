from bsb.core import from_storage
import numpy as np
from bsb.storage import Storage
from bsb.plotting import plot_network
import csv


network = from_storage("MouseCerebellum.hdf5")

plot_network(network)
