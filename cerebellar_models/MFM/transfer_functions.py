"""
New Implementation for automatic derivation of MFM transfer functions.
"""

import json
import multiprocessing as mp
import os

import nest
import numpy as np
import yaml
import matplotlib.pyplot as plt
from bsb import from_storage
from cerebellar_models.analysis import StructureReport
from matplotlib.colors import LinearSegmentedColormap


def _connect_neuron(neuron, sim_info, rates_dict, rng):
    """Wire one neuron to its Poisson generators."""
    for conn_tag, info in sim_info['connections'].items():
        rate = float(rates_dict[conn_tag])
        if rate <= 0:
            continue
        mean_conv = info['convergence']
        std_conv = info['convergence_std']
        n_syn = max(1, int(round(rng.normal(mean_conv, std_conv) if std_conv > 0 else mean_conv)))
        n_syn_per_pair = max(1, int(round(info.get('synapses_per_pair', 1.0))))
        # n_syn_per_pair contacts from the same pre-neuron → scale weight, keep PG count = n_syn
        effective_weight = float(info['weight']) * n_syn_per_pair
        pgs = nest.Create('poisson_generator', n_syn, params={'rate': rate})
        nest.Connect(
            pgs, neuron,
            conn_spec={'rule': 'all_to_all'},
            syn_spec={
                'synapse_model': 'static_synapse',
                'weight': effective_weight,
                'delay': float(info['delay']),
                'receptor_type': int(info['receptor_type']),
            },
        )


def _run_single_sim(sim_info, rates_dict, seed, duration, transient):
    """Single-neuron sim — used only for raster plots."""
    cell_params = sim_info['cell_params']
    nest.ResetKernel()
    nest.set(resolution=0.1, print_time=False, rng_seed=seed)
    nest.Install('cerebmodule')
    neuron = nest.Create('eglif_cond_alpha_multisyn', 1, params=cell_params)
    rng = np.random.default_rng(seed)
    nest.SetStatus(neuron, {"V_m": float(rng.uniform(cell_params["V_reset"], cell_params["V_th"] - 1.0))})
    _connect_neuron(neuron, sim_info, rates_dict, rng)
    sr = nest.Create('spike_recorder')
    nest.Connect(neuron, sr)
    nest.Simulate(duration)
    spikes = sr.get('events')['times']
    spikes = spikes[spikes > transient]
    return len(spikes) / ((duration - transient) / 1000.0), spikes


def _compute_point(args):
    """
    One NEST sim with n_reps independent neurons — avoids n_reps ResetKernel calls.
    Each neuron gets its own PG set and initial V_m drawn from a per-rep seed.
    """
    sim_info, rates_dict, n_reps, base_seed, duration, transient = args
    cell_params = sim_info['cell_params']

    nest.ResetKernel()
    nest.set(resolution=0.1, print_time=False, rng_seed=base_seed)
    nest.Install('cerebmodule')

    neurons = nest.Create('eglif_cond_alpha_multisyn', n_reps, params=cell_params)
    srs = nest.Create('spike_recorder', n_reps)

    for i in range(n_reps):
        rng = np.random.default_rng(base_seed + i)
        neuron = neurons[i:i + 1]
        nest.SetStatus(neuron, {"V_m": float(rng.uniform(cell_params["V_reset"], cell_params["V_th"] - 1.0))})
        _connect_neuron(neuron, sim_info, rates_dict, rng)
        nest.Connect(neuron, srs[i:i + 1])

    nest.Simulate(duration)

    rep_rates = np.zeros(n_reps)
    for i in range(n_reps):
        spikes = srs[i].get('events')['times']
        rep_rates[i] = np.sum(spikes > transient) / ((duration - transient) / 1000.0)

    return rep_rates.mean(), rep_rates.std()



class TransferFunctionSimulator:
    """
    Computes the transfer function of a single-compartment NEST neuron
    driven by Poisson inputs, given a BSB scaffold and circuit configuration.
    Supports multiprocessing: set ``n_workers`` to use a parallel pool.
    """

    def __init__(
        self,
        scaffold_path,
        yaml_path,
        tf_dict_path,
        cell_name,
        n_points=20,
        n_reps=10,
        duration=10000.0,
        transient=1000.0,
        base_seed=1234,
        n_workers=1,
        show_rasters=False,
        raster_max_points=10,
    ):
        self.n_points = n_points
        self.n_reps = n_reps
        self.duration = duration
        self.transient = transient
        self.base_seed = base_seed
        self.n_workers = n_workers if n_workers > 0 else os.cpu_count()
        self.show_rasters = show_rasters
        self.raster_max_points = raster_max_points

        self.tf_mean = None
        self.tf_std = None
        self.freq_axes = None
        self.tags = None

        scaffold = from_storage(scaffold_path)
        report = StructureReport(scaffold)
        conn_table = report.plots['connectivity_table']
        conn_table.update()

        self._convergences = {
            ps.tag: (np.mean(conv_arr), np.std(conv_arr))
            for ps, conv_arr in zip(
                scaffold.get_connectivity_sets(),
                conn_table.get_convergences().values(),
            )
        }
        self._synapses_per_pair = {
            ps.tag: np.mean(syn_arr)
            for ps, syn_arr in zip(
                scaffold.get_connectivity_sets(),
                conn_table.get_nb_synapse_per_pair().values(),
            )
        }

        self.cell_name = cell_name
        self.sim_info = self._build_sim_info(yaml_path, tf_dict_path, cell_name)

    _CELL_COLORS = {
        "granule_cell":   [0.7, 0.15, 0.15, 1.0],
        "golgi_cell":     [0, 0.45, 0.7, 1.0],
        "purkinje_cell":  [0.275, 0.800, 0.275, 1.0],
        "basket_cell":    [1, 0.647, 0, 1.0],
        "stellate_cell":  [1, 0.84, 0, 1.0],
        "dcn_p":          [0.3, 0.3, 0.3, 1.0],
        "dcn_i":          [0.635, 0, 0.145, 1.0],
        "io":             [0.46, 0.376, 0.54, 1.0],
    }

    def _cell_cmap(self):
        color = self._CELL_COLORS.get(self.cell_name, [0, 0, 0, 1.0])
        return LinearSegmentedColormap.from_list(self.cell_name, [(1, 1, 1, 1), color])

    def _build_sim_info(self, yaml_path, tf_dict_path, cell_name):
        with open(yaml_path) as f:
            circuit = yaml.safe_load(f)
        with open(tf_dict_path) as f:
            tf_dict = json.load(f)

        cell_data = tf_dict[cell_name]
        conn_models = circuit['simulations']['nest_basal_activity']['connection_models']
        cell_params = circuit['simulations']['nest_basal_activity']['cell_models'][cell_name]['constants']

        connections = {}
        for conn_tag, conn_info in cell_data['connections'].items():
            s = conn_models[conn_tag]['synapses'][0]
            connections[conn_tag] = {
                'label': conn_info['label'],
                'convergence': self._convergences[conn_tag][0],
                'convergence_std': self._convergences[conn_tag][1],
                'synapses_per_pair': self._synapses_per_pair[conn_tag],
                'weight': s['weight'],
                'delay': s['delay'],
                'receptor_type': s.get('receptor_type', 1),
                'rate_range': conn_info['rate_range'],
            }

        label_ranges = {}
        for conn_tag, info in connections.items():
            label, rate_range = info['label'], tuple(info['rate_range'])
            if label in label_ranges:
                if label_ranges[label] != rate_range:
                    raise ValueError(
                        f"Connections with label '{label}' have conflicting rate_ranges: "
                        f"{list(label_ranges[label])} vs {info['rate_range']}"
                    )
            else:
                label_ranges[label] = rate_range

        return {'cell_params': cell_params, 'connections': connections}

    def _setup_grid(self):
        label_to_range = {}
        for info in self.sim_info['connections'].values():
            if info['label'] not in label_to_range:
                label_to_range[info['label']] = info['rate_range']
        unique_labels = list(label_to_range.keys())
        label_to_idx = {lbl: i for i, lbl in enumerate(unique_labels)}
        freq_axes = [
            np.linspace(label_to_range[l][0], label_to_range[l][1], self.n_points)
            for l in unique_labels
        ]
        return unique_labels, label_to_idx, freq_axes

    def _rates_dict(self, idx, label_to_idx, freq_axes):
        return {
            conn_tag: freq_axes[label_to_idx[info['label']]][idx[label_to_idx[info['label']]]]
            for conn_tag, info in self.sim_info['connections'].items()
        }


    def _run_nest_sim(self, rates_dict, seed):
        return _run_single_sim(self.sim_info, rates_dict, seed, self.duration, self.transient)



    def compute_tf(self, save=True):
        unique_labels, label_to_idx, freq_axes = self._setup_grid()
        shape = tuple(len(ax) for ax in freq_axes)
        n_points_total = int(np.prod(shape))

        print(
            f"TF grid: {dict(zip(unique_labels, shape))} "
            f"→ {n_points_total} points  ({self.n_reps} neurons/sim)"
            + (f"  [{self.n_workers} workers]" if self.n_workers > 1 else "")
        )

        all_indices = list(np.ndindex(*shape))
        all_args = [
            (self.sim_info, self._rates_dict(idx, label_to_idx, freq_axes), self.n_reps, self.base_seed, self.duration, self.transient)
            for idx in all_indices
        ]

        tf_mean = np.zeros(shape)
        tf_std = np.zeros(shape)
        total = len(all_indices)

        if self.n_workers == 1:
            for done, (idx, args) in enumerate(zip(all_indices, all_args), start=1):
                mean, std = _compute_point(args)
                tf_mean[idx] = mean
                tf_std[idx] = std
                if self.show_rasters and done <= self.raster_max_points:
                    rates_dict = self._rates_dict(idx, label_to_idx, freq_axes)
                    _, spikes = self._run_nest_sim(rates_dict, self.base_seed)
                    self.plot_raster(spikes, mean, rates_dict, 0, idx)
                print(f"[{done}/{total}] -> {mean:.2f} ± {std:.2f} Hz", end="\r")
        else:
            with mp.Pool(processes=self.n_workers) as pool:
                for done, (idx, result) in enumerate(
                        zip(all_indices, pool.imap(_compute_point, all_args)), start=1
                ):
                    tf_mean[idx] = result[0]
                    tf_std[idx] = result[1]
                    print(f"[{done}/{total}] -> {result[0]:.2f} ± {result[1]:.2f} Hz", end="\r")

        print()
        self.tf_mean = tf_mean
        self.tf_std = tf_std
        self.freq_axes = freq_axes
        self.tags = unique_labels

        if save:
            out_dir = "output_TF"
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"TF_{self.cell_name}.npz")
            np.savez(
                out_path,
                tf_mean=tf_mean,
                tf_std=tf_std,
                tags=np.array(unique_labels),
                **{f"freq_axis_{lbl}": ax for lbl, ax in zip(unique_labels, freq_axes)},
            )
            print(f"TF saved → {out_path}")



    def plot_tf(self, save=True, n_slices=5):
        ndim = len(self.tags)
        if ndim == 2:
            return self._plot_tf_2d(save)
        else:
            return self._plot_tf_3d(save, n_slices)

    def _plot_tf_2d(self, save):
        labels = self.tags
        extent = [
            self.freq_axes[0][0], self.freq_axes[0][-1],
            self.freq_axes[1][0], self.freq_axes[1][-1],
        ]
        cmap = self._cell_cmap()

        fig1, ax1 = plt.subplots(figsize=(5, 4))
        im1 = ax1.imshow(self.tf_mean.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
        ax1.set_title("Mean firing rate")
        ax1.set_xlabel(labels[0])
        ax1.set_ylabel(labels[1])
        fig1.colorbar(im1, ax=ax1, label="Hz")
        fig1.tight_layout()

        fig2, ax2 = plt.subplots(figsize=(5, 4))
        im2 = ax2.imshow(self.tf_std.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
        ax2.set_title("STD across reps")
        ax2.set_xlabel(labels[0])
        ax2.set_ylabel(labels[1])
        fig2.colorbar(im2, ax=ax2, label="Hz")
        fig2.tight_layout()

        if save:
            fig1.savefig(f'TF_{self.cell_name}.png', dpi=150)
            fig2.savefig(f'TF_std_{self.cell_name}.png', dpi=150)

        return fig1, fig2

    def _plot_tf_3d(self, save, n_slices):
        labels = self.tags          # [l0, l1, l2]
        ax0, ax1, ax2 = self.freq_axes
        cmap = self._cell_cmap()

        slice_indices = np.round(np.linspace(0, len(ax2) - 1, n_slices)).astype(int)
        extent = [ax0[0], ax0[-1], ax1[0], ax1[-1]]

        vmin_mean = self.tf_mean.min()
        vmax_mean = self.tf_mean.max()
        vmin_std  = self.tf_std.min()
        vmax_std  = self.tf_std.max()

        fig, axes = plt.subplots(2, n_slices, figsize=(4 * n_slices, 7))
        fig.suptitle(f"{self.cell_name} — Transfer Function", fontsize=13)

        for col, k in enumerate(slice_indices):
            val = ax2[k]

            im1 = axes[0, col].imshow(
                self.tf_mean[:, :, k].T,
                origin='lower', aspect='auto', extent=extent,
                cmap=cmap, vmin=vmin_mean, vmax=vmax_mean,
            )
            axes[0, col].set_title(f"{labels[2]}={val:.1f} Hz", fontsize=9)
            axes[0, col].set_xlabel(labels[0])
            if col == 0:
                axes[0, col].set_ylabel(labels[1])
            fig.colorbar(im1, ax=axes[0, col], label="Hz" if col == n_slices - 1 else "")

            im2 = axes[1, col].imshow(
                self.tf_std[:, :, k].T,
                origin='lower', aspect='auto', extent=extent,
                cmap='plasma', vmin=vmin_std, vmax=vmax_std,
            )
            axes[1, col].set_xlabel(labels[0])
            if col == 0:
                axes[1, col].set_ylabel(labels[1])
            fig.colorbar(im2, ax=axes[1, col], label="Hz" if col == n_slices - 1 else "")

        axes[0, 0].annotate("Mean FR", xy=(0, 0.5), xycoords='axes fraction',
                            fontsize=10, rotation=90, va='center', ha='right')
        axes[1, 0].annotate("STD", xy=(0, 0.5), xycoords='axes fraction',
                            fontsize=10, rotation=90, va='center', ha='right')

        fig.tight_layout()

        if save:
            fig.savefig(f'TF_{self.cell_name}.png', dpi=150)

        return fig

    def plot_raster(self, spikes, rate, rates_dict, rep, idx=None):
        fig, ax = plt.subplots(figsize=(10, 2))
        ax.eventplot(spikes, lineoffsets=1, linelengths=0.8)
        rate_str = " ".join(f"{k}={v:.1f}Hz" for k, v in rates_dict.items())
        ax.set_title(f"rep={rep} | idx={idx} | FR={rate:.2f} Hz | {rate_str}")
        ax.set_xlim(0, self.duration)
        ax.set_yticks([])
        ax.set_xlabel("Time [ms]")
        plt.tight_layout()
        plt.show()