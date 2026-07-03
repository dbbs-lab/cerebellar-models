import multiprocessing as mp
import os
import random
from copy import deepcopy

import numpy as np
from deap import algorithms, base, creator, tools
from deap.tools.emo import assignCrowdingDist

from cerebellar_models.optimization.features import *
from cerebellar_models.optimization.utils import _suppress_output

try:
    from deap.base import clone as _deap_clone
except:
    try:
        from deap.tools import clone as _deap_clone
    except:
        _deap_clone = deepcopy

PARALLEL_SETTINGS = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "OMP_PROC_BIND": "FALSE",
    "OMP_WAIT_POLICY": "PASSIVE",
    "NEST_NUM_THREADS": "1",
    "NEST_NO_BANNER": "1",
    "OMPI_MCA_pml": "ob1",
    "OMPI_MCA_btl": "self,tcp",
    "OMPI_MCA_btl_base_verbose": "0",
}


def _worker_init():
    for k in list(os.environ):
        if k.startswith(("OMPI_", "PMI_", "PMIX_", "MPI_")):
            os.environ.pop(k, None)
    os.environ.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    import random as _rnd

    import numpy as _np

    pid = os.getpid()
    _rnd.seed((pid * 69069) & 0xFFFFFFFF)
    _np.random_seed = (pid * 1013904223) & 0xFFFFFFFF
    _np.random.seed(_np.random_seed)


def _clamp(x, lb, ub):
    if lb is not None and x < lb:
        return lb
    if ub is not None and x > ub:
        return ub
    return x


class Constraint:
    def __init__(self, func, name=None, ctx=None):
        self.func = func
        self.name = name or getattr(func, "__name__", "constraint")
        self.ctx = ctx or {}

    def __call__(self, ind):
        return self.func(ind, self.ctx)


class ArchiveND:
    def __init__(self, eps_vec, tau_vec, cap=1000):
        self.eps_vec = list(eps_vec)
        self.tau_vec = list(tau_vec)
        self.dim = len(self.eps_vec)
        self.cap = cap
        self.items = []

    def _cell(self, fit):
        return tuple(int(fit[i] / self.eps_vec[i]) for i in range(self.dim))

    def _is_feasible(self, fit):
        return all(fit[i] <= self.tau_vec[i] for i in range(self.dim))

    def update(self, candidates):
        feas = [ind for ind in candidates if self._is_feasible(ind.fitness.values)]
        if not feas:
            return
        nd = tools.sortNondominated(feas, k=len(feas), first_front_only=True)[0]
        grid = {}
        for ind in self.items + nd:
            c = self._cell(ind.fitness.values)
            cur = grid.get(c)
            if (cur is None) or ind.fitness.dominates(cur.fitness):
                grid[c] = ind
        self.items = list(grid.values())
        if len(self.items) > self.cap:
            self.items = tools.selNSGA2(self.items, self.cap)

    def sample(self, k=None):
        if not self.items:
            return []
        if k is None or k >= len(self.items):
            return list(self.items)
        return random.sample(self.items, k)


def _apply_constraints(ind_, constraints):
    for c in constraints:
        ind_ = c(ind_)
    return ind_


class Repair:
    def __init__(self, constraints):
        self.constraints = list(constraints) if constraints else []

    def __call__(self, ind):
        return _apply_constraints(ind, self.constraints) if self.constraints else ind


class EvaluateDecorator:
    def __init__(self, repair, evaluate):
        self.repair = repair
        self.evaluate = evaluate

    def __call__(self, ind):
        ind = self.repair(ind)
        return self.evaluate(ind)


class MateDecorator:
    def __init__(self, repair, mate):
        self.repair = repair
        self.mate = mate

    def __call__(self, ind1, ind2):
        off1, off2 = self.mate(ind1, ind2)
        off1 = self.repair(off1)
        off2 = self.repair(off2)
        return off1, off2


class MutateDecorator:
    def __init__(self, repair, mutate):
        self.repair = repair
        self.mutate = mutate

    def __call__(self, ind):
        res = self.mutate(ind)
        if isinstance(res, tuple):
            return tuple(self.repair(x) for x in res)
        else:
            return self.repair(res)


class Optimizer(object):
    def __init__(
        self,
        multicomp_data,
        protocol: dict,
        nest_model: str,
        model_params: dict,
        opt_params: list,
        fitness: dict,
        bounds,
        archive=None,
        init_params_fn=None,
        evaluate_fn=None,
        constraints=None,
        knee_weights=None,
    ):

        self.multicomp_data = multicomp_data
        self.protocol = protocol
        self.nest_model = nest_model
        self.model_params = model_params
        self.opt_params = opt_params
        self.fitness = fitness
        self.init_params_fn = init_params_fn
        self.evaluate_fn = evaluate_fn
        self.archive = archive
        self.constraints = list(constraints) if constraints else []
        self.bounds = bounds
        self.POP_SIZE = 400
        self.NO_IMPROVE_PATIENCE = 15
        self.N_GEN = 120
        self.ETA_CROSSOVER = 10.0
        self.ETA_MUTATION = 20.0
        self.CROSSOVER_PB = 0.6
        self.MUTATION_PB = 0.4
        self.BEST_HISTORY_FILE = "best_history.pkl"
        self.print_evolution = False
        self._fit_names = list(self.fitness.keys())
        self._fit_hist = {name: [] for name in self._fit_names}
        self._hv_hist = []
        self.knee_weights = np.array(
            knee_weights if knee_weights is not None else np.ones(len(fitness)), dtype=float
        )
        self.output_path = None
        self.return_pareto = False
        self.FI_SNAPSHOT_INTERVAL = 10

    def _extract_multicomp_features(self):
        neuron_features = multicomp_features(
            self.multicomp_data,
            threshold=self.protocol["threshold"],
            start_stim=self.protocol["start_stim"],
            end_stim=self.protocol["end_stim"],
        )
        return neuron_features

    def _nest_single_sim(self, current, params=None):

        with _suppress_output():
            import nest

            try:
                nest.set_verbosity("M_FATAL")
            except:
                try:
                    nest.SetKernelStatus({"print_time": False})
                except:
                    pass
            nest.ResetKernel()
            nest.SetKernelStatus({"local_num_threads": 1})
            try:
                nest.Install("cerebmodule")
            except:
                raise Exception("NEST module not installed.")
            nest.resolution = 0.1

        import nest

        cell_params = self.model_params if params is None else params
        try:
            cell = nest.Create(self.nest_model, 1, params=cell_params)
        except:
            raise Exception("NEST model not found.")
        cell.V_m = cell_params["E_L"]
        gen = nest.Create("dc_generator")
        nest.SetStatus(
            gen,
            {
                "amplitude": current,
                "start": self.protocol["start_stim"],
                "stop": self.protocol["end_stim"],
            },
        )
        sr = nest.Create("spike_recorder")
        nest.Connect(gen, cell)
        nest.Connect(cell, sr)
        nest.Simulate(self.protocol["duration"])
        return nest.GetStatus(sr, "events")[0]["times"]

    def _nest_protocol(self, cell_params=None, targ_df=None):
        if targ_df is None:
            targ_df = self._extract_multicomp_features()
        return {I: self._nest_single_sim(I, params=cell_params) for I in targ_df["current"].values}

    def _extract_nest_features(self, targ_df=None):
        if targ_df is None:
            targ_df = self._extract_multicomp_features()
        nest_results = self._nest_protocol(targ_df=targ_df)
        spike_trains = [np.array(nest_results[I]) for I in targ_df["current"].values]
        return point_neuron_features(
            spike_trains=spike_trains,
            currents=targ_df["current"].values,
            start_stim=self.protocol["start_stim"],
            end_stim=self.protocol["end_stim"],
            duration=self.protocol.get("duration"),
        )

    def _set_parallel(self):
        os.environ.update(PARALLEL_SETTINGS)
        mp.set_start_method("spawn", force=True)
        avail = os.cpu_count() or 1
        default_workers = max(1, avail - 2)
        N_WORKERS = int(os.environ.get("NWORKERS", default_workers))
        N_WORKERS = max(1, min(N_WORKERS, avail))
        pool = mp.Pool(processes=N_WORKERS, initializer=_worker_init)
        return pool

    @staticmethod
    def init_params_wrapper(optimizer, *args, **kwargs):
        if optimizer.init_params_fn is not None:
            return optimizer.init_params_fn(optimizer, *args, **kwargs)
        return None

    @staticmethod
    def evaluate_wrapper(optimizer, individual, *args, **kwargs):
        if optimizer.evaluate_fn is not None:
            return optimizer.evaluate_fn(optimizer, individual, *args, **kwargs)
        return None

    def constrained_mate(self, ind1, ind2):
        n = len(ind1)
        assert n == len(self.opt_params)
        lower_bounds = [self.bounds[name][0] for name in self.opt_params]
        upper_bounds = [self.bounds[name][1] for name in self.opt_params]
        p1 = list(ind1)
        p2 = list(ind2)
        for j in range(n):
            p1[j] = _clamp(p1[j], lower_bounds[j], upper_bounds[j])
            p2[j] = _clamp(p2[j], lower_bounds[j], upper_bounds[j])
        tools.cxSimulatedBinaryBounded(
            p1, p2, eta=self.ETA_CROSSOVER, low=lower_bounds, up=upper_bounds
        )
        for j in range(n):
            ind1[j] = p1[j]
            ind2[j] = p2[j]
        return ind1, ind2

    def constrained_mutate(self, ind, indpb=0.12):
        n = len(ind)
        assert n == len(self.opt_params)
        lower_bounds = [self.bounds[name][0] for name in self.opt_params]
        upper_bounds = [self.bounds[name][1] for name in self.opt_params]
        p = list(ind)
        for j in range(n):
            p[j] = _clamp(p[j], lower_bounds[j], upper_bounds[j])
        (p_mut,) = tools.mutPolynomialBounded(
            p, eta=self.ETA_MUTATION, low=lower_bounds, up=upper_bounds, indpb=indpb
        )
        for j in range(n):
            ind[j] = p_mut[j]
        return (ind,)

    def pareto_selectorND(self, pop, k):
        feas = [ind for ind in pop if self.archive._is_feasible(ind.fitness.values)]
        infeas = [ind for ind in pop if ind not in feas]
        arch_feas = [a for a in self.archive.items if self.archive._is_feasible(a.fitness.values)]
        tot_feas = feas + arch_feas
        if len(tot_feas) >= k:
            return tools.selNSGA2(tot_feas, k)

        def _evaluate_violation(ind):
            return sum(
                max(ind.fitness.values[i] - self.archive.tau_vec[i], 0.0)
                for i in range(self.archive.dim)
            )

        infeas.sort(key=_evaluate_violation)
        infeas_num = k - len(tot_feas)
        infeas_pool = infeas[: max(3 * infeas_num, infeas_num)]
        chosen = tot_feas + (tools.selNSGA2(infeas_pool, infeas_num) if infeas_pool else [])
        return tools.selNSGA2(chosen, k)

    def add_constraint(self, constraint: Constraint):
        self.constraints.append(constraint)

    def pick_knee(self, front, rho=1e-4):
        F = np.array([ind.fitness.values for ind in front], dtype=float)
        w = self.knee_weights
        asf = np.max(F / (w + 1e-12), axis=1) + rho * np.sum(F / (w + 1e-12), axis=1)
        i = int(np.argmin(asf))
        return front[i], float(asf[i])

    def jitter(self, pop, toolbox, alpha=0.05):
        lower_bounds = [self.bounds[name][0] for name in self.opt_params]
        upper_bounds = [self.bounds[name][1] for name in self.opt_params]
        R = [upper_bounds[j] - lower_bounds[j] for j in range(len(lower_bounds))]
        elite = tools.selBest(pop, len(pop) // 5)
        jittered = [toolbox.clone(ind) for ind in elite]
        for ind in jittered:
            for j in range(len(ind)):
                sigma = alpha * R[j]
                ind[j] = _clamp(ind[j] + random.gauss(0.0, sigma), lower_bounds[j], upper_bounds[j])
        pop = toolbox.select(pop + jittered, self.POP_SIZE)
        return pop

    def log_fitness_history(self, pop):
        vals = [ind.fitness.values for ind in pop if ind.fitness.valid]
        if not vals:
            return
        A = np.array(vals, dtype=float)
        A = np.nan_to_num(A, nan=np.nan, posinf=np.nan, neginf=np.nan)
        means = np.nanmean(A, axis=0)
        for i, name in enumerate(self._fit_names):
            self._fit_hist[name].append(float(means[i]))

    def log_hypervolume(self, pareto_front):
        """Compute and store hypervolume of the Pareto front.

        Reference point is [1.0] * n_obj (worst possible for all objectives).
        HV increases as the front improves in convergence or diversity.
        """
        try:
            ref = [1.0 + 1e-6] * len(self._fit_names)
            pts = [list(ind.fitness.values) for ind in pareto_front]
            pts = [[max(0.0, v) for v in p] for p in pts]
            hv = tools.hypervolume(pts, ref)
            self._hv_hist.append(float(hv))
        except Exception:
            self._hv_hist.append(float("nan"))

    def plot_fitness_history(self, fname="fitness_history.png"):
        from itertools import cycle

        import matplotlib.pyplot as plt

        L = min(len(v) for v in self._fit_hist.values() if v)
        gens = np.arange(1, L + 1, dtype=int)
        markers = cycle(["o", "s", "^", "d", "x", "*", "v", ">", "<", "p", "h"])

        fig, ax1 = plt.subplots()
        for name in self._fit_names:
            y = np.asarray(self._fit_hist[name][:L], dtype=float)
            ax1.plot(gens, y, "-" + next(markers), label=name)
        ax1.set_xlabel("Generation")
        ax1.set_ylabel("Mean error")
        ax1.set_title("Population mean fitness per generation")

        lines1, labs1 = ax1.get_legend_handles_labels()
        ax1.legend(lines1, labs1, fontsize=7, loc="upper right")
        fig.tight_layout()
        if self.output_path:
            plt.savefig(os.path.join(self.output_path, fname), dpi=100)
        else:
            plt.savefig(fname, dpi=100)
        plt.close()

    def log_best(self, gen, ind, reset=False):
        import pickle

        if getattr(self, "output_path", None):
            path = os.path.join(self.output_path, "best_history.pkl")
        else:
            path = getattr(self, "BEST_HISTORY_FILE", "best_history.pkl")
        mode = "wb" if reset or (not os.path.exists(path)) else "ab"

        rec = {
            "gen": int(gen),
            "params": list(ind),
            "fitness": tuple(ind.fitness.values),
        }
        with open(path, mode) as f:
            pickle.dump(rec, f)

    def plot_fI_snapshot(
        self, snapshots, fname="fI_snapshots.png", use_global=False, print_params=True
    ):
        import matplotlib.pyplot as plt
        import numpy as np

        def _rate_for_plot(spike_times, t0, t1, T, use_global_flag):
            st = np.asarray(spike_times if spike_times is not None else [], dtype=float)
            if use_global_flag:
                n = st.size
                win_s = T / 1000.0
            else:
                n = np.sum((st >= t0) & (st < t1))
                win_s = (t1 - t0) / 1000.0
            return float(n) / win_s if win_s > 0 else 0.0

        neuron_data = self._extract_multicomp_features()
        currents = neuron_data["current"].values.astype(float)
        t0, t1, T = (
            self.protocol["start_stim"],
            self.protocol["end_stim"],
            self.protocol["duration"],
        )

        def _peak_to_rate(peak_val, curr):
            st = []
            if peak_val is None:
                pass
            elif isinstance(peak_val, (list, tuple, np.ndarray)):
                st.extend(np.asarray(peak_val, dtype=float).ravel().tolist())
            else:
                try:
                    st.append(float(peak_val))
                except Exception:
                    pass
            st = np.asarray(st, dtype=float)
            st = st[np.isfinite(st)]
            if curr == 0.0:
                win_s = T / 1000.0
            else:
                st = st[(st >= t0) & (st <= t1)]
                win_s = (t1 - t0) / 1000.0
            return float(st.size) / win_s if win_s > 0 else 0.0

        target_freqs = np.array(
            [
                _peak_to_rate(row["peak_time"], float(row["current"]))
                for _, row in neuron_data.iterrows()
            ],
            dtype=float,
        )

        plt.figure(figsize=(7, 5))
        plt.plot(currents, target_freqs, "o-", label="Target")

        last_gen = None
        last_params = None
        for i, (gen, params) in enumerate(snapshots):
            sim = self._nest_protocol(params, targ_df=neuron_data)
            nest_freqs = np.array(
                [_rate_for_plot(sim[I], t0, t1, T, use_global) for I in currents], dtype=float
            )
            is_last = i == len(snapshots) - 1
            plt.plot(
                currents,
                nest_freqs,
                "s-" if is_last else "-",
                linewidth=2.0 if is_last else 1.0,
                alpha=1.0 if is_last else 0.25,
                label=(f"NEST best (gen {gen})" if is_last else None),
            )
            if is_last:
                last_gen = gen
                last_params = params

        title = f"f–I snapshots until gen {last_gen}" if last_gen is not None else "f–I snapshots"
        plt.xlabel("Injected current [pA]")
        plt.ylabel("Firing rate [Hz]")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        if self.output_path:
            savefig_path = os.path.join(self.output_path, fname)
            plt.savefig(savefig_path, dpi=300)
        else:
            plt.savefig(fname, dpi=300)
        plt.close()

        if print_params and (last_params is not None):
            msg = [f"[fI snapshot] gen {last_gen} – optimal params:"]
            for k in self.opt_params:
                if k in last_params:
                    try:
                        v = float(last_params[k])
                        msg.append(f"  {k} = {v:.6g}")
                    except Exception:
                        msg.append(f"  {k} = {last_params[k]}")
            print("\n".join(msg))

    def plot_pareto_evolution(self, pareto_snapshots, fname="pareto_evolution.png"):
        """Radar chart (knee evolution) + parallel coordinates (full Pareto front).

        pareto_snapshots: list of (gen, front_fit, knee_fit)
          front_fit : np.array (n_front, n_obj) — fitness of all Pareto front members
          knee_fit  : np.array (n_obj,)         — fitness of the knee solution

        How to read
        -----------
        Radar  — each spoke = one objective (centre=0 perfect, edge=1 worst).
                 Each polygon = knee at one generation snapshot.
                 Polygon shrinks toward centre as optimisation converges.
                 Spokes that stop shrinking = stuck objectives.

        Parallel coordinates — each vertical axis = one objective (0 bottom,
                 1 top). Each line = one Pareto front individual; colour maps
                 generation (dark=early, bright=late).
                 Crossing lines between two axes → trade-off between them.
                 Parallel lines → compatible objectives (improve together).
                 Lines bunched at bottom of an axis → that objective solved.
                 Lines spread top-to-bottom → main source of diversity.
                 Thick dashed line = current knee.
        """
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if not pareto_snapshots:
            return

        n_obj = len(self._fit_names)
        short_pc = [n.replace("_error", "").replace("_", " ") for n in self._fit_names]
        short_rd = [n.replace("_error", "").replace("_", "\n") for n in self._fit_names]
        n_snaps = len(pareto_snapshots)
        cmap = cm.get_cmap("plasma_r", max(n_snaps, 2))
        x = np.arange(n_obj)

        fig = plt.figure(figsize=(16, 6))
        gs = GridSpec(
            1,
            3,
            figure=fig,
            left=0.05,
            right=0.96,
            top=0.87,
            bottom=0.13,
            wspace=0.38,
        )
        ax_radar = fig.add_subplot(gs[0, 0], polar=True)
        ax_pc = fig.add_subplot(gs[0, 1:])

        # ── Radar chart ────────────────────────────────────────────────
        angles = np.linspace(0, 2 * np.pi, n_obj, endpoint=False).tolist()
        angles_closed = angles + angles[:1]

        ax_radar.set_theta_offset(np.pi / 2)
        ax_radar.set_theta_direction(-1)
        ax_radar.set_rlim(0, 1)
        ax_radar.set_rgrids([0.25, 0.5, 0.75, 1.0], fontsize=6, alpha=0.4)
        ax_radar.set_thetagrids(np.degrees(angles), short_rd, fontsize=7)
        ax_radar.set_title("Knee — radar", fontsize=9, pad=14)

        step = max(1, n_snaps // 8)
        radar_idx = list(range(0, n_snaps, step))
        if n_snaps - 1 not in radar_idx:
            radar_idx.append(n_snaps - 1)

        for ri, si in enumerate(radar_idx):
            gen, _, knee_fit = pareto_snapshots[si]
            t = si / max(n_snaps - 1, 1)
            color = cmap(t)
            is_last = si == n_snaps - 1
            vals = list(knee_fit) + [knee_fit[0]]
            label = f"gen {gen}" if (ri == 0 or is_last) else None
            ax_radar.plot(
                angles_closed,
                vals,
                color=color,
                linewidth=2.5 if is_last else 1.0,
                alpha=1.0 if is_last else 0.5,
                label=label,
            )
            ax_radar.fill(angles_closed, vals, color=color, alpha=0.15 if is_last else 0.05)

        ax_radar.legend(
            loc="lower left",
            bbox_to_anchor=(-0.3, -0.22),
            fontsize=7,
            ncol=2,
        )

        # ── Parallel coordinates ───────────────────────────────────────
        for h in [0.25, 0.5, 0.75]:
            ax_pc.axhline(h, color="0.90", linewidth=0.7, zorder=0)
        for xi in x:
            ax_pc.axvline(xi, color="0.78", linewidth=1.0, zorder=0)

        for snap_idx, (gen, front_fit, knee_fit) in enumerate(pareto_snapshots):
            t = snap_idx / max(n_snaps - 1, 1)
            color = cmap(t)
            is_last = snap_idx == n_snaps - 1
            alpha_front = 0.07 + 0.25 * t

            for row in front_fit:
                ax_pc.plot(x, row, color=color, alpha=alpha_front, linewidth=0.6, zorder=1)

            ax_pc.plot(
                x,
                knee_fit,
                color=color,
                linewidth=2.8 if is_last else 1.2,
                alpha=1.0 if is_last else 0.6,
                linestyle="--",
                marker="o",
                markersize=5 if is_last else 2,
                zorder=3,
                label=f"knee gen {gen}" if is_last else None,
            )

        ax_pc.set_xticks(x)
        ax_pc.set_xticklabels(short_pc, fontsize=8)
        ax_pc.set_ylim(-0.02, 1.05)
        ax_pc.set_ylabel("Error  [0 = perfect · 1 = worst]", fontsize=9)
        ax_pc.set_title("Pareto front — parallel coordinates", fontsize=9)
        ax_pc.legend(fontsize=8, loc="upper right")
        ax_pc.spines["top"].set_visible(False)
        ax_pc.spines["right"].set_visible(False)

        sm = cm.ScalarMappable(
            cmap="plasma_r",
            norm=plt.Normalize(pareto_snapshots[0][0], pareto_snapshots[-1][0]),
        )
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax_pc, pad=0.01, aspect=30, shrink=0.8)
        cb.set_label("Generation", fontsize=8)

        last_gen = pareto_snapshots[-1][0]
        fig.suptitle(
            f"Pareto front evolution — up to gen {last_gen}",
            fontsize=11,
        )
        if self.output_path:
            plt.savefig(os.path.join(self.output_path, fname), dpi=150, bbox_inches="tight")
        else:
            plt.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close()

    def set_optimizer(self, init_kwargs=None, eval_kwargs=None):
        init_kwargs = init_kwargs or {}
        eval_kwargs = eval_kwargs or {}
        weights = tuple(self.fitness.values())
        if hasattr(creator, "Fitness"):
            old_len = len(getattr(creator.Fitness, "weights", []))
            new_len = len(weights)
            if old_len != new_len:
                del creator.Fitness
                if hasattr(creator, "Individual"):
                    del creator.Individual
        creator.create("Fitness", base.Fitness, weights=weights)
        creator.create("Individual", list, fitness=creator.Fitness)

        toolbox = base.Toolbox()
        pool = self._set_parallel()
        toolbox.register("map", pool.map)
        toolbox.register("clone", _deap_clone)
        toolbox.register("initialize_params", self.init_params_wrapper, self, **init_kwargs)
        toolbox.register(
            "individual", tools.initIterate, creator.Individual, toolbox.initialize_params
        )
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", self.evaluate_wrapper, self, **eval_kwargs)
        toolbox.register("mate", self.constrained_mate)
        toolbox.register("mutate", self.constrained_mutate)
        toolbox.register("select", lambda pop, k: self.pareto_selectorND(pop, k))

        repair = Repair(self.constraints)
        toolbox.evaluate = EvaluateDecorator(repair, toolbox.evaluate)
        toolbox.mate = MateDecorator(repair, toolbox.mate)
        toolbox.mutate = MutateDecorator(repair, toolbox.mutate)
        return toolbox

    def optimize(self):
        toolbox = self.set_optimizer()
        pop = [creator.Individual(toolbox.initialize_params()) for _ in range(self.POP_SIZE)]
        for ind, fit in zip(pop, toolbox.map(toolbox.evaluate, pop)):
            fit = np.nan_to_num(np.array(fit, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
            ind.fitness.values = tuple(fit)

        best_individual, best_fitness, best_knee_score = None, np.inf, np.inf
        best_history = []
        snapshots = []
        pareto_front_snapshots = []
        best_gen = []
        stale = 0

        for gen in range(self.N_GEN):
            if self.print_evolution:
                print(f"Executing generation {gen+1}/{self.N_GEN}")

            self.archive.update(pop)
            parents_source = pop[:] + self.archive.sample(
                k=min(len(pop) // 4, len(self.archive.items))
            )
            assignCrowdingDist(parents_source)
            n_parents = len(pop) if len(pop) % 2 == 0 else len(pop) - 1
            parents = tools.selTournamentDCD(parents_source, n_parents)
            offspring = algorithms.varOr(
                parents,
                toolbox,
                lambda_=self.POP_SIZE,
                cxpb=self.CROSSOVER_PB,
                mutpb=self.MUTATION_PB,
            )

            for ind, fit in zip(offspring, toolbox.map(toolbox.evaluate, offspring)):
                fit = np.nan_to_num(np.array(fit, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
                ind.fitness.values = tuple(fit)

            pop = toolbox.select(pop + offspring, self.POP_SIZE)
            self.log_fitness_history(pop)
            if self.print_evolution:
                self.plot_fitness_history()

            pareto_front = tools.sortNondominated(pop, k=len(pop), first_front_only=True)[0]
            self.log_hypervolume(pareto_front)
            best_gen, knee_score = self.pick_knee(pareto_front)

            if knee_score < best_knee_score:
                best_knee_score = knee_score
                best_individual = toolbox.clone(best_gen)
                best_fitness = tuple(best_individual.fitness.values)
                best_history.append(best_fitness)
                stale = 0
            else:
                stale += 1
            self.log_best(gen + 1, toolbox.clone(best_gen), reset=(gen == 0))

            if stale >= self.NO_IMPROVE_PATIENCE:
                if self.print_evolution:
                    print("Add small jitter to top individuals to escape stagnation.")
                pop = self.jitter(pop, toolbox, alpha=0.1)
                stale = 0
            if (gen + 1) % self.FI_SNAPSHOT_INTERVAL == 0:
                genes = best_individual if best_individual is not None else best_gen
                cell_params_best = {
                    **self.model_params,
                    **{k: float(v) for k, v in zip(self.opt_params, genes)},
                }
                snapshots.append((gen + 1, cell_params_best))
                self.plot_fI_snapshot(snapshots)
                front_fit = np.array([ind.fitness.values for ind in pareto_front], dtype=float)
                knee_fit = np.array(best_gen.fitness.values, dtype=float)
                pareto_front_snapshots.append((gen + 1, front_fit, knee_fit))
                if self.print_evolution:
                    self.plot_pareto_evolution(pareto_front_snapshots)
            if (gen + 1) % 50 == 0 and self.return_pareto:
                if self.output_path is not None:
                    save_dir = os.path.join(self.output_path, "pareto_snapshots")
                else:
                    save_dir = "pareto_snapshots"
                os.makedirs(save_dir, exist_ok=True)
                pareto_front = tools.sortNondominated(pop, k=len(pop), first_front_only=True)[0]
                genes = np.array([list(ind) for ind in pareto_front], dtype=float)
                np.save(os.path.join(save_dir, f"pareto_genes_{gen+1}.npy"), genes)
                fitns = np.array([ind.fitness.values for ind in pareto_front], dtype=float)
                np.save(os.path.join(save_dir, f"pareto_fitness_{gen+1}.npy"), fitns)
        return best_individual, best_fitness
