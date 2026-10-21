"""
SNN basic simulation to validate MFM simulations and fix the range in TF.json
"""
from os import makedirs, remove
from os.path import basename, dirname, join, normpath, splitext
from glob import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from quantities import ms
from bsb import from_storage, parse_configuration_file, MPI
from elephant.statistics import instantaneous_rate
from elephant.kernels import GaussianKernel
from elephant.signal_processing import hilbert
from elephant.phase_analysis import phase_locking_value

from cerebellar_models.analysis.spike_plots import BasicSimulationReport
from cerebellar_models.analysis.spiking_results import SpikingResults


PROTOCOLS = {
    "poisson_step": {
        "device": "poisson_generator",
        "params": {},
    },
    "sinusoidal_poisson": {
        "device": "sinusoidal_poisson_generator",
        "params": {"amplitude": None, "frequency": 10.0, "phase": 0.0},
    },
}

# presynaptic token (connection name prefix before "_to_") -> recorded population column
PRESYN_TO_POP = {
    "glomerulus": "mossy_fibers",
    "golgi": "golgi_cell",
    "ascending_axon": "granule_cell",
    "parallel_fiber": "granule_cell",
    "basket": "basket_cell",
    "stellate": "stellate_cell",
}

class SNNSimulate_forMFM:

    def __init__(self, mf_stims, config_path, scaffold_path, tf_json = None, output_path=None,
                 stim_duration = 1000., transient=1000.):
        self.mf_stims = mf_stims
        self.config_path = config_path
        self.scaffold_path = scaffold_path
        self.transient = transient
        self.stim_duration = stim_duration
        self.output_path = output_path
        self.tf_json = tf_json
        self.scaffold = None
        self.config = None
        self.sim_name = None
        self.simulation = None
        self.windows = []
        self.rate_df = None
        self.spike_data = None

    def configure(self):
        self.scaffold = from_storage(self.scaffold_path)
        self.config = parse_configuration_file(self.config_path)
        self.scaffold.configuration = self.config
        return self.scaffold

    def _select_basal_sim(self):
        sim_name = next(
            name for name in self.scaffold.configuration.simulations
            if "basal_activity" in name
        )
        sim = self.scaffold.configuration.simulations[sim_name]
        self.sim_name = sim_name
        self.simulation = sim

        bg = sim.devices.pop("background_noise")
        self._bg_rate, self._bg_weight, self._bg_delay = bg.rate, bg.weight, bg.delay
        self._targetting = {"strategy": "cell_model", "cell_models": ["mossy_fibers"]}
        return sim

    def _background(self, start, stop):
        return {
            "device": "poisson_generator",
            "rate": self._bg_rate,
            "start": start,
            "stop": stop,
            "targetting": self._targetting,
            "weight": self._bg_weight,
            "delay": self._bg_delay,
        }

    def prepare_sim(self):
        sim = self._select_basal_sim()

        def stimulus(rate, start, stop):
            return {
                "device": "poisson_generator",
                "rate": float(rate),
                "start": start,
                "stop": stop,
                "targetting": self._targetting,
                "weight": self._bg_weight,
                "delay": self._bg_delay,
            }

        # initial transient: only background noise
        t = 0.0
        sim.devices["mf_background_0"] = self._background(t, t + self.transient)
        t += self.transient

        # one stim window per rate, each followed by a background-only gap
        self.windows = []
        for i, rate in enumerate(self.mf_stims):
            sim.devices[f"mf_stim_{i}"] = stimulus(rate, t, t + self.stim_duration)
            self.windows.append((float(rate), t, t + self.stim_duration))
            t += self.stim_duration
            sim.devices[f"mf_background_{i + 1}"] = self._background(t, t + self.transient)
            t += self.transient

        # total length of the concatenated protocol
        sim.duration = t
        return sim

    def load_results(self):
        self.configure()
        self.prepare_sim()
        return self.scaffold

    def _nio_path(self):
        if self.output_path is None:
            raise ValueError("output_path must be set to read/write the simulation results")
        name = basename(normpath(self.output_path))
        if name.endswith(".nio"):
            name = name[:-4]
        return join(self.output_path, name + ".nio")

    def simulate(self, mode="ow"):
        result = self.scaffold.run_simulation(self.sim_name)
        if self.output_path is not None:
            folder = dirname(self._nio_path()) or "."
            stem = splitext(basename(self._nio_path()))[0]
            if MPI.get_rank() == 0:
                makedirs(folder, exist_ok=True)
                for old in glob(join(folder, "*.nio")):
                    remove(old)  # drop stale results so they aren't merged in
            MPI.barrier()
            result.write(join(folder, f"{stem}_rank{MPI.get_rank()}.nio"), mode)
        MPI.barrier()
        return result

    def print_results(self, pdf_dir=None, ignored_ct=None):
        # only rank 0 reads back the results and renders the figures
        if MPI.get_rank() != 0:
            MPI.barrier()
            return
        nio_path = self._nio_path()
        folder_nio = dirname(nio_path) or "."
        pdf_dir = pdf_dir if pdf_dir is not None else folder_nio
        makedirs(pdf_dir, exist_ok=True)
        base = splitext(basename(nio_path))[0]

        # full simulation: all steps
        report = BasicSimulationReport(
            scaffold=self.scaffold,
            simulation_name=self.sim_name,
            folder_nio=folder_nio,
            time_from=0,
            time_to=self.simulation.duration,
            ignored_ct=ignored_ct,
        )
        report.print_report(join(pdf_dir, f"{base}_full.pdf"))

        for i, (rate, start, stop) in enumerate(self.windows):
            report_step = BasicSimulationReport(scaffold=self.scaffold,
                                                simulation_name=self.sim_name,
                                                folder_nio=folder_nio,
                                                time_from=start,
                                                time_to=stop,
                                                ignored_ct=ignored_ct)
            report_step.print_report(join(pdf_dir, f"{base}_step{i}_{rate:g}Hz.pdf"))
            del report_step
        MPI.barrier()

    def compute_firing_rates(self, save_spiketrains=False, ignored_ct=None):
        # only rank 0 reads back the results
        if MPI.get_rank() != 0:
            MPI.barrier()
            return None
        folder_nio = dirname(self._nio_path()) or "."
        sr = SpikingResults(
            scaffold=self.scaffold,
            simulation_name=self.sim_name,
            time_from=0,
            time_to=self.simulation.duration,
            folder_nio=folder_nio,
            ignored_ct=ignored_ct,
        )

        records = []
        self.spike_data = [] if save_spiketrains else None
        for rate, start, stop in self.windows:
            # each stim window is measured in full; the transient of only
            # background noise that precedes every step (initial settling /
            # inter-step gap) is excluded by construction, since it lies
            # outside [start, stop]
            # set time_to before time_from to satisfy the start <= stop check
            sr.time_to = stop
            sr.time_from = start
            spikes = sr.filt_spikes  # one merged SpikeTrain per population
            duration_s = (stop - start) / 1e3

            row = {"mf_rate": rate}
            step_trains = {}
            for st, pop, n in zip(spikes, sr.populations, sr.nb_neurons):
                # mean firing rate = n_spikes / duration / n_cells
                row[pop] = st.size / duration_s / n if n > 0 else np.nan
                if save_spiketrains:
                    step_trains[pop] = st
            records.append(row)
            if save_spiketrains:
                self.spike_data.append(step_trains)

        self.rate_df = pd.DataFrame.from_records(records)
        if self.tf_json is not None:
            self.update_tf_json()
        MPI.barrier()
        return self.rate_df

    def update_tf_json(self, path=None):
        if MPI.get_rank() != 0:
            return
        if self.rate_df is None:
            raise ValueError("run compute_firing_rates() before updating the TF json")
        path = path if path is not None else self.tf_json
        with open(path) as f:
            tf = json.load(f)

        df = self.rate_df
        mf_min, mf_max = round(min(self.mf_stims)), round(max(self.mf_stims))
        # mean firing response of each population at the lowest and highest mf drive
        row_min = df.loc[df["mf_rate"].idxmin()]
        row_max = df.loc[df["mf_rate"].idxmax()]

        for target in tf.values():
            for conn_name, conn in target.get("connections", {}).items():
                if conn.get("label") == "mf":
                    conn["rate_range"] = [mf_min, mf_max]
                    continue
                presyn = conn_name.split("_to_")[0]
                pop = PRESYN_TO_POP.get(presyn)
                if pop is None or pop not in df.columns:
                    continue
                r0, r1 = float(row_min[pop]), float(row_max[pop])
                conn["rate_range"] = [round(min(r0, r1)), round(max(r0, r1))]

        with open(path, "w") as f:
            json.dump(tf, f, indent=2)
        return tf

    def plot_raster(self, pops=None, step=0, output_path=None, markersize=1.0):
        if MPI.get_rank() != 0:
            return None
        if not self.spike_data:
            raise ValueError(
                "No spike trains saved; run compute_firing_rates(save_spiketrains=True) first"
            )
        step_trains = self.spike_data[step]
        pops = list(step_trains.keys()) if pops is None else pops

        fig, ax = plt.subplots(figsize=(12, 6))
        offset = 0
        yticks, ylabels = [], []
        for pop in pops:
            st = step_trains[pop]
            senders = st.array_annotations["senders"]
            _, inv = np.unique(senders, return_inverse=True)
            n = inv.max() + 1 if inv.size else 0
            ax.scatter(st.times.magnitude, inv + offset, s=markersize, label=pop)
            yticks.append(offset + n / 2)
            ylabels.append(pop)
            offset += n

        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels)
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("population")
        ax.set_title(f"raster step {step} - mf {self.windows[step][0]:g} Hz")
        fig.tight_layout()
        if output_path is not None:
            fig.savefig(output_path)
        return fig, ax

    def _load_pop_spikes(self, ignored_ct=None):
        # cell type -> (merged SpikeTrain, 0-based senders) over the whole simulation
        folder_nio = dirname(self._nio_path()) or "."
        total = max(stop for _, _, stop in self.windows)
        sr = SpikingResults(
            scaffold=self.scaffold,
            simulation_name=self.sim_name,
            time_from=0,
            time_to=total,
            folder_nio=folder_nio,
            ignored_ct=ignored_ct if ignored_ct is not None else [],
        )
        pop_spikes = {}
        for pop, st in zip(sr._populations, sr._all_spikes):
            senders = np.asarray(st.array_annotations["senders"])
            pop_spikes[pop] = (st, senders - senders.min())
        return pop_spikes

    def _presyn_phases(self, pop_spikes, pre_type, pre_ids, start, stop,
                       dt, kernel_sigma, min_spikes):
        # Hilbert phase of each presynaptic neuron's instantaneous rate in the window
        merged, senders = pop_spikes[pre_type]
        phases = []
        for nid in pre_ids:
            st = merged[np.where(senders == nid)[0]].time_slice(start * ms, stop * ms)
            if st.size < min_spikes:
                continue
            rate = instantaneous_rate(
                st, sampling_period=dt * ms,
                kernel=GaussianKernel(sigma=kernel_sigma * ms), border_correction=True,
            )
            phases.append(np.angle(hilbert(rate)).ravel())
        return phases

    @staticmethod
    def _ensemble_plv(phases):
        # mean PLV over all distinct presynaptic pairs (upper triangle, no diagonal)
        n = len(phases)
        if n < 2:
            return np.nan
        vals = [
            phase_locking_value(phases[i], phases[j])
            for i in range(n) for j in range(i + 1, n)
        ]
        return float(np.mean(vals))

    def compute_plv(self, connections, output_path="plv.json", min_conv=2, min_spikes=2,
                    dt=0.1, kernel_sigma=20., seed=None, ignored_ct=None):
        """
        For each connection, pick at random a significant postsynaptic neuron (with at
        least ``min_conv`` converging partners) and compute the phase locking value of
        its presynaptic ensemble for every stimulation step, plus the mean across steps.
        Results are written to ``output_path`` (plv.json).
        """
        if MPI.get_rank() != 0:
            MPI.barrier()
            return None
        rng = np.random.default_rng(seed)
        pop_spikes = self._load_pop_spikes(ignored_ct)

        results = {}
        for conn in connections:
            cs = self.scaffold.get_connectivity_set(conn)
            pre_locs, post_locs = cs.load_connections().all()
            pre, post = pre_locs[:, 0], post_locs[:, 0]
            if cs.pre_type_name not in pop_spikes:
                continue  # presynaptic population was not recorded
            # map each post cell to its converging presynaptic ids, keep eligible ones
            pre_per_post = {p: pre[post == p] for p in np.unique(post)}
            eligible = [p for p, ids in pre_per_post.items() if len(ids) >= min_conv]
            if not eligible:
                continue
            post_id = int(rng.choice(eligible))
            pre_ids = pre_per_post[post_id]

            per_freq = {}
            for freq, start, stop in self.windows:
                phases = self._presyn_phases(
                    pop_spikes, cs.pre_type_name, pre_ids, start, stop,
                    dt, kernel_sigma, min_spikes,
                )
                per_freq[f"{freq:g}"] = self._ensemble_plv(phases)

            valid = [v for v in per_freq.values() if not np.isnan(v)]
            results[conn] = {
                "pre_type": cs.pre_type_name,
                "post_type": cs.post_type_name,
                "post_neuron_id": post_id,
                "n_pre": int(len(pre_ids)),
                "plv_per_frequency": per_freq,
                "plv_mean": float(np.mean(valid)) if valid else np.nan,
            }

        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)
        MPI.barrier()
        return results


class SNNSimulate_forValidation(SNNSimulate_forMFM):
    def __init__(self, protocol, config_path, scaffold_path, rate,
                 tf_json=None, output_path=None,
                 stim_duration=1000., transient=1000., **protocol_params):
        if protocol not in PROTOCOLS:
            raise ValueError(
                f"protocol must be one of {sorted(PROTOCOLS)}, got {protocol!r}"
            )
        super().__init__(
            mf_stims=[rate], config_path=config_path, scaffold_path=scaffold_path,
            tf_json=tf_json, output_path=output_path,
            stim_duration=stim_duration, transient=transient,
        )
        self.protocol = protocol
        self.rate = float(rate)

        allowed = PROTOCOLS[protocol]["params"]
        unknown = set(protocol_params) - set(allowed)
        if unknown:
            raise ValueError(
                f"{protocol!r} accepts only {sorted(allowed)}, "
                f"got unexpected {sorted(unknown)}"
            )
        self.protocol_params = {}
        for name, default in allowed.items():
            val = protocol_params.get(name, default)
            self.protocol_params[name] = float(val if val is not None else self.rate)

    def _stim_device(self, start, stop):
        return {
            "device": PROTOCOLS[self.protocol]["device"],
            "rate": self.rate,
            **self.protocol_params,
            "start": start,
            "stop": stop,
            "targetting": self._targetting,
            "weight": self._bg_weight,
            "delay": self._bg_delay,
        }

    def prepare_sim(self):
        sim = self._select_basal_sim()

        t = 0.0
        sim.devices["mf_background_0"] = self._background(t, t + self.transient)
        t += self.transient

        start, stop = t, t + self.stim_duration
        sim.devices["mf_validation"] = self._stim_device(start, stop)
        # window label: the driving frequency if the protocol has one, else the rate
        label = self.protocol_params.get("frequency", self.rate)
        self.windows = [(float(label), start, stop)]
        t = stop

        sim.devices["mf_background_1"] = self._background(t, t + self.transient)
        t += self.transient

        sim.duration = t
        return sim
