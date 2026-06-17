"""
New Implementation for automatic derivation of MFM transfer functions.
"""

import json
import multiprocessing as mp
import os

os.environ.setdefault("PYNEST_QUIET", "1")
import nest
nest.set_verbosity("M_WARNING")
import numpy as np
import yaml
import matplotlib.pyplot as plt
from tqdm import tqdm
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
        if sim_info.get('stochastic_convergence', False):
            std_conv = info['convergence_std']
            n_syn = max(1, int(round(rng.normal(mean_conv, std_conv) if std_conv > 0 else mean_conv)))
        else:
            n_syn = max(1, int(round(mean_conv)))
        n_syn_per_pair = max(1, int(round(info.get('synapses_per_pair', 1.0))))
        effective_weight = float(info['weight']) * n_syn_per_pair
        pg = nest.Create('poisson_generator', 1, params={'rate': rate * n_syn})
        nest.Connect(
            pg, neuron,
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
    nest.set_verbosity("M_WARNING")
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
    nest.set_verbosity("M_WARNING")
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


_WORKER_SIM_INFO = None


def _worker_init(sim_info):
    """Pool initializer: store sim_info once per worker process."""
    global _WORKER_SIM_INFO
    _WORKER_SIM_INFO = sim_info
    nest.set_verbosity("M_WARNING")


def _compute_point_worker(args):
    """Parallel entry point: reads sim_info from the process global."""
    rates_dict, n_reps, base_seed, duration, transient = args
    return _compute_point((_WORKER_SIM_INFO, rates_dict, n_reps, base_seed, duration, transient))



class TransferFunction:
    """
    Holds a computed or loaded transfer function and provides plot/summary methods.
    """

    _CELL_COLORS = {
        "granule_cell":  [0.7, 0.15, 0.15, 1.0],
        "golgi_cell":    [0, 0.45, 0.7, 1.0],
        "purkinje_cell": [0.275, 0.800, 0.275, 1.0],
        "purkinje_cell_Z+" : [0.275, 0.800, 0.275, 1.0],
        "purkinje_cell_Z-" : [0.275, 0.800, 0.275, 1.0],
        "basket_cell":   [1, 0.647, 0, 1.0],
        "stellate_cell": [1, 0.84, 0, 1.0],
        "MLI":           [1, 0.7435, 0, 1.0],
        "dcn_p":         [0.3, 0.3, 0.3, 1.0],
        "dcn_i":         [0.635, 0, 0.145, 1.0],
        "io":            [0.46, 0.376, 0.54, 1.0],
    }

    def __init__(self, cell_name, tf_mean, tf_std, tags, freq_axes):
        self.cell_name = cell_name
        self.tf_mean   = np.asarray(tf_mean)
        self.tf_std    = np.asarray(tf_std)
        self.tags      = list(tags)
        self.freq_axes = [np.asarray(a) for a in freq_axes]

    @classmethod
    def load(cls, cell_name, output_dir='output_TF'):
        """Load a previously saved TF from ``output_dir/TF_{cell_name}.npz``."""
        path = os.path.join(output_dir, f'TF_{cell_name}.npz')
        data = np.load(path, allow_pickle=True)
        tags = data['tags'].tolist()
        freq_axes = [data[f'freq_axis_{t}'] for t in tags]
        return cls(cell_name, data['tf_mean'], data['tf_std'], tags, freq_axes)

    def save(self, output_dir='output_TF'):
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f'TF_{self.cell_name}.npz')
        np.savez(
            path,
            tf_mean=self.tf_mean,
            tf_std=self.tf_std,
            tags=np.array(self.tags),
            **{f'freq_axis_{t}': ax for t, ax in zip(self.tags, self.freq_axes)},
        )
        print(f"TF saved → {path}")


    def summary(self):
        shape = self.tf_mean.shape
        n_total = self.tf_mean.size

        def _fmt(idx):
            return ', '.join(
                f"{t}={self.freq_axes[i][idx[i]]:.1f}"
                for i, t in enumerate(self.tags)
            )

        idx_max     = np.unravel_index(self.tf_mean.argmax(), shape)
        idx_min     = np.unravel_index(self.tf_mean.argmin(), shape)
        idx_max_std = np.unravel_index(self.tf_std.argmax(),  shape)

        print(f"Cell    : {self.cell_name.replace('_', ' ')}")
        print(f"Inputs  : " + ', '.join(
            f"{t} [{ax[0]:.1f}–{ax[-1]:.1f} Hz, {len(ax)} pts]"
            for t, ax in zip(self.tags, self.freq_axes)
        ))
        print(f"Grid    : {' × '.join(str(s) for s in shape)} = {n_total} points")
        print(f"FR mean : {self.tf_mean.mean():.2f} ± {self.tf_mean.std():.2f} Hz")
        print(f"FR max  : {self.tf_mean.max():.2f} Hz  @ ({_fmt(idx_max)})")
        print(f"FR min  : {self.tf_mean.min():.2f} Hz  @ ({_fmt(idx_min)})")
        print(f"Std mean: {self.tf_std.mean():.2f} Hz  |  max {self.tf_std.max():.2f} Hz @ ({_fmt(idx_max_std)})")


    def plot_tf(self, save=True, n_slices=5):
        if len(self.tags) == 2:
            return self._plot_tf_2d(save)
        else:
            return self._plot_tf_3d_scatter(save)

    def _cell_cmap(self):
        color = self._CELL_COLORS.get(self.cell_name, [0, 0, 0, 1.0])
        return LinearSegmentedColormap.from_list(self.cell_name, [(1, 1, 1, 1), color])

    def _plot_tf_2d(self, save):
        labels = self.tags
        extent = [
            self.freq_axes[0][0], self.freq_axes[0][-1],
            self.freq_axes[1][0], self.freq_axes[1][-1],
        ]
        cmap = self._cell_cmap()
        cell_label = self.cell_name.replace('_', ' ')

        fig1, ax1 = plt.subplots(figsize=(5, 4))
        im1 = ax1.imshow(self.tf_mean.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
        ax1.set_title("Mean firing rate")
        ax1.set_xlabel(f"{labels[0]} [Hz]")
        ax1.set_ylabel(f"{labels[1]} [Hz]")
        fig1.colorbar(im1, ax=ax1, label=f"{cell_label} [Hz]")
        fig1.tight_layout()

        mean_std = self.tf_std.mean()
        fig2, ax2 = plt.subplots(figsize=(5, 4))
        im2 = ax2.imshow(self.tf_std.T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
        ax2.set_title(f"STD across reps (mean: {mean_std:.2f} Hz)")
        ax2.set_xlabel(f"{labels[0]} [Hz]")
        ax2.set_ylabel(f"{labels[1]} [Hz]")
        fig2.colorbar(im2, ax=ax2, label=f"{cell_label} [Hz]")
        fig2.tight_layout()

        if save:
            fig1.savefig(f'TF_{self.cell_name}.png', dpi=150)
            fig2.savefig(f'TF_std_{self.cell_name}.png', dpi=150)

        return fig1, fig2

    def _make_3d_scatter(self, data, cmap, colorbar_label, save_path=None):
        ax0, ax1, ax2 = self.freq_axes
        labels = self.tags
        n0, n1, n2 = len(ax0), len(ax1), len(ax2)

        x = np.arange(n0)[:, None, None] * np.ones((1, n1, n2))
        y = np.arange(n1)[None, :, None] * np.ones((n0, 1, n2))
        z = np.arange(n2)[None, None, :] * np.ones((n0, n1, 1))
        c = data.ravel()

        fig = plt.figure()
        ax3d = fig.add_subplot(111, projection='3d')
        sc = ax3d.scatter(x.ravel(), y.ravel(), z.ravel(), s=50, c=c, cmap=cmap)
        ax3d.grid(False)

        ax3d.set_xlabel(f'{labels[0]} [Hz]', fontsize=10)
        ax3d.set_ylabel(f'{labels[1]} [Hz]', fontsize=10)
        ax3d.set_zlabel(f'{labels[2]} [Hz]', fontsize=12, labelpad=10)

        def _ticks(freq_axis, n_ticks=4):
            idx = np.round(np.linspace(0, len(freq_axis) - 1, n_ticks)).astype(int)
            return idx.tolist(), [f'{freq_axis[i]:.1f}' for i in idx]

        xt, xl = _ticks(ax0); ax3d.set_xticks(xt); ax3d.set_xticklabels(xl, fontsize=9)
        yt, yl = _ticks(ax1); ax3d.set_yticks(yt); ax3d.set_yticklabels(yl, fontsize=9)
        zt, zl = _ticks(ax2); ax3d.set_zticks(zt); ax3d.set_zticklabels(zl, fontsize=9)

        ticks = np.linspace(c.min(), c.max(), 5, endpoint=True)
        cb = fig.colorbar(sc, ax=ax3d, ticks=ticks, shrink=0.5, pad=0.15)
        cb.set_label(colorbar_label, labelpad=10, rotation=270 + 180)
        cb.ax.tick_params(labelsize=9)

        if save_path:
            fig.savefig(save_path, dpi=150)

        return fig

    def _plot_tf_3d_scatter(self, save):
        cell_label = self.cell_name.replace('_', ' ')
        fig_mean = self._make_3d_scatter(
            self.tf_mean, self._cell_cmap(), f'{cell_label} [Hz]',
            f'TF_{self.cell_name}.png' if save else None,
        )
        fig_std = self._make_3d_scatter(
            self.tf_std, self._cell_cmap(), f'{cell_label} std [Hz]',
            f'TF_std_{self.cell_name}.png' if save else None,
        )
        return fig_mean, fig_std



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
        stochastic_convergence=False,
    ):
        self.n_points = n_points
        self.n_reps = n_reps
        self.duration = duration
        self.transient = transient
        self.base_seed = base_seed
        self.n_workers = n_workers if n_workers > 0 else os.cpu_count()
        self.show_rasters = show_rasters
        self.raster_max_points = raster_max_points
        self.stochastic_convergence = stochastic_convergence

        self._tf = None

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

    @staticmethod
    def _load_yaml(yaml_path):
        with open(yaml_path) as f:
            data = yaml.safe_load(f)
        if '$import' not in data:
            return data
        import_info = data.pop('$import')
        ref = import_info['ref'].split('#')[0]
        base_path = ref if os.path.isabs(ref) else os.path.join(os.path.dirname(os.path.abspath(yaml_path)), ref)
        base = TransferFunctionSimulator._load_yaml(base_path)
        return TransferFunctionSimulator._deep_merge(base, data)

    @staticmethod
    def _deep_merge(base, override):
        result = dict(base)
        for k, v in override.items():
            if k in result and isinstance(result[k], dict) and isinstance(v, dict):
                result[k] = TransferFunctionSimulator._deep_merge(result[k], v)
            else:
                result[k] = v
        return result

    def _build_sim_info(self, yaml_path, tf_dict_path, cell_name):
        circuit = self._load_yaml(yaml_path)
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

        return {
            'cell_params': cell_params,
            'connections': connections,
            'stochastic_convergence': self.stochastic_convergence,
        }

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

    def _args_generator(self, all_indices, label_to_idx, freq_axes):
        """Yield lightweight args tuples one at a time (no sim_info copy per item)."""
        for idx in all_indices:
            rates_dict = self._rates_dict(idx, label_to_idx, freq_axes)
            yield (rates_dict, self.n_reps, self.base_seed, self.duration, self.transient)

    def compute_tf(self, save=True):
        """Run simulations and return a ``TransferFunction`` object."""
        unique_labels, label_to_idx, freq_axes = self._setup_grid()
        shape = tuple(len(ax) for ax in freq_axes)
        n_points_total = int(np.prod(shape))

        print(
            f"TF grid: {dict(zip(unique_labels, shape))} "
            f"→ {n_points_total} points  ({self.n_reps} neurons/sim)"
            + (f"  [{self.n_workers} workers]" if self.n_workers > 1 else "")
        )

        all_indices = list(np.ndindex(*shape))
        tf_mean = np.zeros(shape)
        tf_std  = np.zeros(shape)

        bar_fmt = "{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"

        if self.n_workers == 1:
            with tqdm(total=len(all_indices), desc=self.cell_name, unit="sim", bar_format=bar_fmt) as pbar:
                for done, idx in enumerate(all_indices, start=1):
                    rates_dict = self._rates_dict(idx, label_to_idx, freq_axes)
                    mean, std = _compute_point(
                        (self.sim_info, rates_dict, self.n_reps, self.base_seed, self.duration, self.transient)
                    )
                    tf_mean[idx] = mean
                    tf_std[idx]  = std
                    if self.show_rasters and done <= self.raster_max_points:
                        _, spikes = self._run_nest_sim(rates_dict, self.base_seed)
                        self.plot_raster(spikes, mean, rates_dict, 0, idx)
                    pbar.update(1)
        else:
            args_gen = self._args_generator(all_indices, label_to_idx, freq_axes)
            with mp.Pool(
                processes=self.n_workers,
                initializer=_worker_init,
                initargs=(self.sim_info,),
                maxtasksperchild=200,
            ) as pool:
                with tqdm(total=len(all_indices), desc=self.cell_name, unit="sim", bar_format=bar_fmt) as pbar:
                    for idx, result in zip(
                        all_indices, pool.imap(_compute_point_worker, args_gen, chunksize=4)
                    ):
                        tf_mean[idx] = result[0]
                        tf_std[idx]  = result[1]
                        pbar.update(1)

        self._tf = TransferFunction(self.cell_name, tf_mean, tf_std, unique_labels, freq_axes)

        if save:
            self._tf.save()

        return self._tf

    def plot_tf(self, save=True, n_slices=5):
        if self._tf is None:
            raise RuntimeError("Run compute_tf() first.")
        return self._tf.plot_tf(save=save, n_slices=n_slices)

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