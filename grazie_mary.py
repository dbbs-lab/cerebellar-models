from os import listdir
from os.path import isfile, join
from neo import io as nio
import numpy as np
from bsb import from_storage
from elephant.statistics import instantaneous_rate
from quantities import ms
from elephant.kernels import GaussianKernel
import matplotlib.pyplot as plt
from neo import SpikeTrain

def extract_ct_device_name(cell_type: str):
    return cell_type.split("_rec")[0]
def extract_spikes_dict(folder_nio, cell_names):
    spikes_res = []
    cell_dict = {}
    current_id = 0
    ignored_ct=["glomerulus", "ubc_glomerulus"]
    for f in listdir(folder_nio):
        file_ = join(folder_nio, f)
        if isfile(file_) and (".nio" in file_):
            block = nio.NixIO(file_, mode="ro").read_all_blocks()[0]
            spiketrains = block.segments[0].spiketrains
            for st in spiketrains:
                cell_type = extract_ct_device_name(st.annotations["device"])
                if cell_type not in cell_names:
                    cell_type += "_cell"
                if cell_type in cell_names and cell_type not in ignored_ct:
                    if cell_type not in cell_dict:
                        cell_dict[cell_type] = {"id": current_id, "senders": []}
                        current_id += 1
                        spikes_res.append([])
                    if isinstance(st.annotations["senders"], np.int64):  # pragma: nocover
                        st.annotations["senders"] = [st.annotations["senders"]]
                    if len(st.annotations["senders"]) > 0:
                        spikes_res[cell_dict[cell_type]["id"]].append(st)
                        cell_dict[cell_type]["senders"].extend(st.annotations["senders"])
    for cell_type in cell_dict:
        cell_dict[cell_type]["senders"] = np.unique(cell_dict[cell_type]["senders"])
    return spikes_res, cell_dict
def load_spikes(scaffold,simulation_name, spikes_res, cell_dict, dt):
    u_gids = []
    u_cell_types = []
    for i, cell_type in enumerate(cell_dict):
        senders = cell_dict[cell_type]["senders"].tolist()
        u_gids.extend(senders)
        u_cell_types.extend([i] * len(senders))
    time_to = scaffold.simulations[simulation_name].duration
    if len(u_gids) == 0:
        return (
            np.zeros((int(time_to / dt), 0), dtype=bool),
            np.array([], dtype=int),
            [],
        )
    sorting = np.argsort(u_gids)
    u_gids = np.array(u_gids)[sorting]
    u_cell_types = np.array(u_cell_types)[sorting]
    inv_convert = np.full(np.max(u_gids) + 1, -1)
    for i, u_gid in enumerate(u_gids):
        inv_convert[u_gid] = i
    tot_num_neuron = len(u_gids)
    all_spikes = np.zeros((int(time_to / dt) + 1, tot_num_neuron), dtype=bool)
    for cell_type in cell_dict:
        for st in spikes_res[cell_dict[cell_type]["id"]]:
            spikes = st.magnitude
            senders = np.array(inv_convert[np.array(st.annotations["senders"])])
            spikes = np.asarray(np.floor(spikes / dt), dtype=int)
            all_spikes[(spikes, senders)] = True
    nb_neurons = np.zeros(len(cell_dict), dtype=int)
    for i, uf in enumerate(cell_dict.keys()):
        nb_neurons[i] = len(cell_dict[uf]["senders"])
    return all_spikes, nb_neurons, list(cell_dict.keys())
def get_filt_spikes(all_spikes, dt, time_from, time_to):
    return all_spikes[int(time_from / dt) : int(time_to / dt)]
def extract_fr(loc_spikes, nb_neurons, time_from, time_to, cell_names):
    firing_rates = {}
    num_filter = len(nb_neurons)
    counts = np.zeros(num_filter + 1)
    counts[1:] = np.cumsum(nb_neurons)
    for i in range(num_filter):
        spikes = loc_spikes[:, int(counts[i]) : int(counts[i + 1])]
        all_fr = np.sum(spikes, axis=0) / ((time_to - time_from) / 1000.0)
        firing_rates[cell_names[i]] = np.mean(all_fr)
    return firing_rates
def plot_fr(loc_spikes, nb_neurons, cell_names,dt, time_from, time_to):
    duration = (time_to-time_from)
    num_filter = len(nb_neurons)
    counts = np.zeros((num_filter + 1))
    counts[1:] = np.cumsum(nb_neurons)
    for i, ct in enumerate(cell_names):
        times, newIds = np.where(loc_spikes[:, int(counts[i]) : int(counts[i + 1])])
        spike_train = SpikeTrain((times*dt) * ms, t_stop=time_to*ms)
        ist_rate = instantaneous_rate(spike_train, sampling_period= 0.1* ms,
                                  kernel=GaussianKernel(10 * ms), border_correction=True)
        plt.figure(figsize=(20,15))
        plt.title(f"Firing rate of {cell_names[i]}")
        plt.plot(ist_rate.times, ist_rate.magnitude / nb_neurons[i], color="red")
        plt.xlabel("Time (ms)")
        plt.ylabel("Firing rate (Hz)")
        plt.xlim(0,duration)
        plt.show()
        plt.close()
def compute_firing_rates(scaffold, folder_nio, simulation_name,time_from, time_to, dt, plot=False):
    cell_names = [scaffold.get_cell_types()[i].name for i in range(len(scaffold.cell_types))]
    spike_res, cell_dict = extract_spikes_dict(folder_nio, cell_names)
    all_spikes, nb_neurons, populations = load_spikes(scaffold,simulation_name, spike_res, cell_dict, dt)
    loc_spikes = get_filt_spikes(all_spikes, dt, time_from=time_from, time_to=time_to)
    if plot:
        plot_fr(loc_spikes, nb_neurons, populations,dt, time_from, time_to)
    return extract_fr(loc_spikes, nb_neurons, time_from, time_to, populations)
if __name__ == "__main__":
    folder_nio = "nio_files/trials_no_io"
    scaffold_name = "mouse_cereb_io_trials.hdf5"
    scaffold = from_storage(scaffold_name)
    simulation_name = "mf_cf_stimulus"
    n_trials = 15
    pc_mean_fr = []
    dcn_p_fr = []
    for i in range(n_trials):
        time_from = (i*1000)+500
        time_to = (i*1000)+760
        dt = 0.1 # nest resolution
        fr = compute_firing_rates(scaffold, folder_nio, simulation_name,time_from, time_to,dt, plot=False)
        print(fr)
        pc_mean_fr.append(fr['purkinje_cell'])
        dcn_p_fr.append(fr['dcn_p'])

    # time_from = 0
    # time_to = 15200
    # dt = 0.1
    # fr_long = compute_firing_rates(scaffold, folder_nio, simulation_name, time_from, time_to, dt, plot=True)

    plt.figure()
    plt.plot(np.arange(1,n_trials+1), pc_mean_fr, '--', marker='o', color='blue', markersize = 5, label ='PC')
    plt.plot(np.arange(1, n_trials + 1), dcn_p_fr, '--', marker='o', color='red', markersize=5, label ='DCNp')
    plt.xlabel("Trials")
    plt.ylabel("Mean firing rate")
    plt.title("PC and DCNp mean firing rate over trials")
    plt.legend()
    plt.xlim(1,15)
    plt.show()

