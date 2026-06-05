import numpy as np
import matplotlib.pyplot as plt
import nest
from cerebellar_models.analysis import StructureReport
from bsb import from_storage
import time


GrC_PARAMS = {
    "A1": 0.01, "A2": -0.94, "C_m": 7, "E_L": -62,
    "E_rev1": 0, "E_rev2": -80, "I_e": -0.888,
    "V_m": -62.0, "V_min": -150, "V_reset": -70, "V_th": -41,
    "k_1": 0.311, "k_2": 0.041407868, "k_adap": 0.022,
    "lambda_0": 1.0, "t_ref": 1.5, "tau_V": 0.3, "tau_m": 24.15,
    "tau_syn1": 5.8, "tau_syn2": 13.61,
}

N_POINTS  = 20
N_REPS    = 10

DURATION  = 10000.0   # 10 s
TRANSIENT = 1000.0    # 1 s

BASE_SEED = 1234

SHOW_RASTERS = False
RASTER_MAX_POINTS = 10


scaffold = from_storage('mouse_cerebellum.hdf5')
report = StructureReport(scaffold)
conn_table = report.plots['connectivity_table']
conn_table.update()

convergences = {
    ps.tag: (np.mean(conv_arr), np.std(conv_arr))
    for ps, conv_arr in zip(
        scaffold.get_connectivity_sets(),
        conn_table.get_convergences().values(),
    )
}

SIM_INFO = {
    'glomerulus_to_granule': {
        'label': 'mf',
        'convergence': convergences['glomerulus_to_granule'][0],
        'convergence_std': convergences['glomerulus_to_granule'][1],
        'weight': 0.23,
        'delay': 1.0,
        'receptor_type': 1,
        'rate_range': [0, 80],
    },
    'golgi_to_glomerulus': {
        'label': 'GoC',
        'convergence': convergences['golgi_to_glomerulus'][0],
        'convergence_std': convergences['golgi_to_glomerulus'][1],
        'weight': 0.36,
        'delay': 2.0,
        'receptor_type': 2,
        'rate_range': [0, 60],
    },
}



def nest_sim(cell_params, sim_info, rates_dict, seed,
             duration=DURATION, transient=TRANSIENT):

    nest.ResetKernel()
    nest.set(resolution=0.1, print_time=False, rng_seed=seed)
    nest.Install('cerebmodule')

    neuron = nest.Create('eglif_cond_alpha_multisyn', 1, params=cell_params)

    # random Vm initialization
    vm0 = np.random.uniform(
        cell_params["V_reset"],
        cell_params["V_th"] - 1.0
    )
    nest.SetStatus(neuron, {"V_m": float(vm0)})

    rng = np.random.default_rng(seed)

    for conn_tag, info in sim_info.items():

        rate = float(rates_dict[conn_tag])
        if rate <= 0:
            continue

        mean_conv = info['convergence']
        std_conv  = info['convergence_std']

        if std_conv > 0:
            n_syn = max(1, int(round(rng.normal(mean_conv, std_conv))))
        else:
            n_syn = max(1, int(round(mean_conv)))

        pgs = nest.Create('poisson_generator', n_syn, params={'rate': rate})

        nest.Connect(
            pgs,
            neuron,
            conn_spec={'rule': 'all_to_all'},
            syn_spec={
                'synapse_model': 'static_synapse',
                'weight': float(info['weight']),
                'delay': float(info['delay']),
                'receptor_type': int(info['receptor_type']),
            },
        )

    sr = nest.Create('spike_recorder')
    nest.Connect(neuron, sr)

    nest.Simulate(duration)

    spikes = sr.get('events')['times']
    spikes = spikes[spikes > transient]

    rate = len(spikes) / ((duration - transient) / 1000.0)

    return rate, spikes



def plot_raster(spikes, rate, rates_dict, rep, idx=None):

    fig, ax = plt.subplots(figsize=(10, 2))

    ax.eventplot(
        spikes,
        lineoffsets=1,
        linelengths=0.8
    )

    rate_str = " ".join(
        f"{k}={v:.1f}Hz" for k, v in rates_dict.items()
    )

    ax.set_title(
        f"rep={rep} | idx={idx} | FR={rate:.2f} Hz | {rate_str}"
    )

    ax.set_xlim(0, DURATION)
    ax.set_yticks([])
    ax.set_xlabel("time (ms)")

    plt.tight_layout()
    plt.show()



def compute_tf(cell_params, sim_info,
               n_points=N_POINTS, n_reps=N_REPS):

    tags = list(sim_info.keys())

    freq_axes = [
        np.linspace(info['rate_range'][0],
                    info['rate_range'][1],
                    n_points)
        for info in sim_info.values()
    ]

    shape = tuple(len(ax) for ax in freq_axes)

    tf_mean = np.zeros(shape)
    tf_std  = np.zeros(shape)

    total = tf_mean.size
    done = 0

    for idx in np.ndindex(*shape):

        rates_dict = {
            tags[d]: freq_axes[d][i]
            for d, i in enumerate(idx)
        }

        rep_rates = np.zeros(n_reps)

        for rep in range(n_reps):

            seed = BASE_SEED + rep

            rate, spikes = nest_sim(
                cell_params,
                sim_info,
                rates_dict,
                seed
            )

            rep_rates[rep] = rate

            if SHOW_RASTERS:
                if done < RASTER_MAX_POINTS:
                    plot_raster(spikes, rate, rates_dict, rep, idx)

        tf_mean[idx] = rep_rates.mean()
        tf_std[idx]  = rep_rates.std()

        done += 1

        print(
            f"[{done}/{total}] {rates_dict} "
            f"-> {tf_mean[idx]:.2f} ± {tf_std[idx]:.2f} Hz",
            end="\r"
        )

    return tf_mean, tf_std, freq_axes, tags


def plot_tf(tf_mean, tf_std, freq_axes, tags, sim_info, save = True):

    labels = [sim_info[t]['label'] for t in tags]

    extent = [
        freq_axes[0][0], freq_axes[0][-1],
        freq_axes[1][0], freq_axes[1][-1],
    ]

    fig1, ax1 = plt.subplots(figsize=(5, 4))
    im1 = ax1.imshow(
        tf_mean.T,
        origin='lower',
        aspect='auto',
        extent=extent,
    )
    ax1.set_title("Mean firing rate")
    ax1.set_xlabel(labels[0])
    ax1.set_ylabel(labels[1])
    fig1.colorbar(im1, ax=ax1, label="Hz")
    fig1.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(5, 4))
    im2 = ax2.imshow(
        tf_std.T,
        origin='lower',
        aspect='auto',
        cmap='plasma',
        extent=extent,
    )
    ax2.set_title("STD across reps")
    ax2.set_xlabel(labels[0])
    ax2.set_ylabel(labels[1])
    fig2.colorbar(im2, ax=ax2, label="Hz")
    fig2.tight_layout()

    if save:
        fig1.savefig('TF.png', dpi=150)
        fig2.savefig('TF_std.png', dpi=150)

    return fig1, fig2


if __name__ == "__main__":

    print("Running TF simulation...")
    t_start = time.time()

    tf_mean, tf_std, freq_axes, tags = compute_tf(
        GrC_PARAMS,
        SIM_INFO
    )

    t_end = time.time()
    elapsed = t_end - t_start
    print(f"\nRun duration: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    fig = plot_tf(tf_mean, tf_std, freq_axes, tags, SIM_INFO)

    plt.show()