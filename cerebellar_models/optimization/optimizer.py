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
        self.knee_weights = np.array(
            knee_weights if knee_weights is not None else np.ones(len(fitness)), dtype=float
        )

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

    def _nest_protocol(self, cell_params=None):
        multicomp_features = self._extract_multicomp_features()
        return {
            I: self._nest_single_sim(I, params=cell_params)
            for I in multicomp_features["current"].values
        }

    def _extract_nest_features(self):
        multicomp_features = self._extract_multicomp_features()
        nest_results = self._nest_protocol()
        spike_trains = []
        for I in multicomp_features["current"].values:
            spike_trains.append(np.array(nest_results[I]))
        nest_features = point_neuron_features(
            spike_trains=spike_trains,
            currents=multicomp_features["current"].values,
            start_stim=self.protocol["start_stim"],
            end_stim=self.protocol["end_stim"],
        )
        return nest_features

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

    def log_fitness_hystory(self, pop):
        vals = [ind.fitness.values for ind in pop if ind.fitness.valid]
        if not vals:
            return
        A = np.array(vals, dtype=float)
        A = np.nan_to_num(A, nan=np.nan, posinf=np.nan, neginf=np.nan)
        means = np.nanmean(A, axis=0)
        for i, name in enumerate(self._fit_names):
            self._fit_hist[name].append(float(means[i]))

    def plot_fitness_hystory(self, fname="fitness_hystory.png"):
        from itertools import cycle

        import matplotlib.pyplot as plt

        L = min(len(v) for v in self._fit_hist.values() if v)
        gens = np.arange(1, L + 1, dtype=int)
        markers = cycle(["o", "s", "^", "d", "x", "*", "v", ">", "<", "p", "h"])
        plt.figure()
        for name in self._fit_names:
            y = np.asarray(self._fit_hist[name][:L], dtype=float)
            plt.plot(gens, y, "-" + next(markers), label=name)
        plt.xlabel("Generation")
        plt.ylabel("Fitness")
        plt.title("Population mean fitness per generation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(fname, dpi=100)
        plt.close()

    def log_best(self, gen, ind, reset=False):
        import pickle

        path = getattr(
            self, "BEST_HISTORY_FILE", getattr(self, "BEST_HISTORY_FILE", "best_history.pkl")
        )
        mode = "wb" if reset or (not os.path.exists(path)) else "ab"
        rec = {
            "gen": int(gen),
            "params": list(ind),
            "fitness": tuple(ind.fitness.values),
        }
        with open(path, mode) as f:
            pickle.dump(rec, f)

    def plot_fI_snapshot(self, snapshots, out_png=None, use_global=False, print_params=True):
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

        target_freqs = (
            neuron_data["mean_frequency"]
            .apply(
                lambda v: (
                    float(np.nanmean(v))
                    if isinstance(v, (list, np.ndarray))
                    else float(v) if v is not None else 0.0
                )
            )
            .values
        )

        plt.figure(figsize=(7, 5))
        plt.plot(currents, target_freqs, "o-", label="Target")

        last_gen = None
        last_params = None
        for i, (gen, params) in enumerate(snapshots):
            sim = self._nest_protocol(params)
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
        plt.savefig(out_png or getattr(self, "FI_SNAPSHOTS_FILE", "fI_snapshots.png"), dpi=300)
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
            self.log_fitness_hystory(pop)
            if self.print_evolution:
                self.plot_fitness_hystory()

            pareto_front = tools.sortNondominated(pop, k=len(pop), first_front_only=True)[0]
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
                pop = self.jitter(pop, toolbox)
                stale = 0
            if (gen + 1) % 10 == 0:
                genes = best_individual if best_individual is not None else best_gen
                cell_params_best = {
                    **self.model_params,
                    **{k: float(v) for k, v in zip(self.opt_params, genes)},
                }
                snapshots.append((gen + 1, cell_params_best))
                self.plot_fI_snapshot(snapshots)

        return best_individual, best_fitness
