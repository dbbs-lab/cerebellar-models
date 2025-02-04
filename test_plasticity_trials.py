import bsb
import nest
from bsb import from_storage, config
import numpy as np
from neo import io
import matplotlib.pyplot as plt
from elephant.statistics import mean_firing_rate, instantaneous_rate
from quantities import ms, s, Hz
from os import listdir
from os.path import isfile, join

# nest.Install("cerebmodule")
# nest.Install("custom_stdp_module")
#
# network = from_storage("mouse_cereb_dcn_io_nest_custom.hdf5")
# n_trials = 15
#
# res = []
# for i in range(n_trials):
#     sim = network.run_simulation("mf_cf_stimulus")
#     paperino = bsb.simulation.results.SimulationResult(sim)
#     np.append(res, paperino)
#     np.append(time, sim.times)
#
# plt.figure()
# plt.plot(f_pc, time, "b")
# plt.show()


# Read simulation data
#my_file_name = "trials_io.nio"
#total_spikes = []
file = "trials_io_3.nio"
# for f in listdir(nio_folder):
#     file = join(nio_folder, f)
#     if isfile(file) and (".nio" in file):
#         sim = io.NixIO(file, mode="ro")
#         block = sim.read_all_blocks()[0]
#         segment = block.segments[0]
#         my_spiketrains = segment.spiketrains
#         total_spikes.append(my_spiketrains.magnitude)
sim = io.NixIO(file, mode="ro")
block = sim.read_all_blocks()[0]
segment = block.segments[0]
my_spiketrains = segment.spiketrains
# total_spikes.append(my_spiketrains.magnitude)



nb_spike_trains = len(my_spiketrains)
fig, ax = plt.subplots(nb_spike_trains, sharex=True, figsize=(10, nb_spike_trains * 6))
for i, spike_t in enumerate(my_spiketrains):  # Iterate over all spike trains
    name = spike_t.annotations["device"]  # Retrieve the device name
    cell_list = spike_t.annotations["senders"]  # Retrieve the ids of the cells spiking
    spike_times = spike_t.magnitude  # Retrieve the spike times
    ax[i].scatter(spike_times, cell_list, c=f"C{i}")
    ax[i].set_xlabel(f"Time ({spike_t.times.units.dimensionality.string})")
    ax[i].set_ylabel(f"Neuron ID")
    ax[i].set_title(f"Spikes from {name}")
plt.tight_layout()
plt.savefig("raster_plot.png", dpi=200)

# for i in range(10):
fr = instantaneous_rate(my_spiketrains[8], sampling_period=0.1*ms)
time = np.arange(0, 15200, 0.1)
plt.figure(figsize = (30,10))
plt.plot(time, fr)
plt.show()