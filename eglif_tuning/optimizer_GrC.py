import logging
import random
import warnings

from cerebellar_models.optimization.features import *
from cerebellar_models.optimization.fitness import *
from cerebellar_models.optimization.optimizer import (
    ArchiveND,
    Constraint,
    Optimizer,
    _apply_constraints,
)
import os

# === WARNING SETTING ===
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
logging.getLogger().setLevel(logging.ERROR)
for name in ("deap", "matplotlib", "numexpr"):
    logging.getLogger(name).setLevel(logging.ERROR)
np.seterr(all="ignore")
os.environ.setdefault("PYTHONWARNINGS", "ignore")


cell_params = {
    "t_ref": 1.5,
    "V_min": -150,
    "C_m": 7,
    "V_th": -41,
    "V_reset": -70,
    "E_L": -62,
    "I_e": -0.888,
    "V_m": -62.0,
    "tau_m": 24.15,
    "k_adap": 0.022,
    "k_1": 0.311,
    "k_2": 1.0 / 24.15,  # oscillatory regime
    "A1": 0.01,
    "A2": -0.94,
}

protocol = {
    "start_stim": 100.0,
    "end_stim": 600.0,
    "duration": 700.0,
    "threshold": cell_params["V_th"],
}


def constrain_bounds(ind, ctx):
    bounds = ctx["bounds"]
    for pos, (a, b) in bounds.items():
        x = float(ind[pos]) if np.isfinite(ind[pos]) else (a + b) / 2.0
        if a > b:
            a, b = b, a
        ind[pos] = min(max(x, a), b)
    return ind


def kadap_constrain(ind, ctx):
    """
    k_adap >= f(k2)
    f(k2) = (C_m/4) * (k2 + 1/tau_m)^2
    """
    idx_kadap = ctx["idx_kadap"]
    C_m = ctx["C_m"]
    tau_m = ctx["tau_m"]
    k2 = 1.0 / tau_m
    kad_low = (C_m / 4.0) * (k2 + 1.0 / tau_m) ** 2
    if ind[idx_kadap] < kad_low:
        ind[idx_kadap] = kad_low
    return ind


def random_init(optimizer):
    ind = [random.uniform(*optimizer.bounds[p]) for p in optimizer.opt_params]
    ind = _apply_constraints(ind, optimizer.constraints)
    return ind


def evaluate_GrC(opt, ind, cell_params=None):
    params = dict(cell_params or opt.model_params)
    for i, p in enumerate(opt.opt_params):
        params[p] = float(ind[i])
    opt.model_params = params

    # Feature extraction
    targ = opt._extract_multicomp_features()
    nest = opt._extract_nest_features()

    # 1. Rheobase error
    rheo, thr = rheobase_loss(nest, targ)
    sel = currents_above_thr(targ, thr)
    # slope = slope_loss(targ, nest, thr, sel)
    curv = curvature_loss(targ, nest, sel, use_max_slope_penalty=True)
    gap = gap_loss(targ, nest, sel, weighted="gaussian")
    return (rheo, curv, gap)


if __name__ == "__main__":
    # Example for GrC
    data_folder = os.path.join(os.path.dirname(__file__), "tofit_eglif/results_tofitEglif/GrC/")

    bounds = {
        "I_e": (-5.0, 5.0),
        "k_adap": (0, 1),
        "k_1": (0.01, 10.0),
        "A1": (0.01, 50.0),
        "A2": (-10.0, 30.0),
    }

    fitness = {
        "rheobase_error": -1,
        "curvature_error": -1,
        "gap_error": -1,
    }

    archive = ArchiveND(eps_vec=[0.05] * 3, tau_vec=[0.15] * 3, cap=1000)

    optimizer = Optimizer(
        multicomp_data=data_folder,
        protocol=protocol,
        nest_model="eglif_multirec_opt",  # not-stochastic version
        opt_params=list(bounds.keys()),
        model_params=cell_params,
        fitness=fitness,
        bounds=bounds,
        archive=archive,
        knee_weights=[1.0, 1.0, 1.0],
    )

    param_positions = {name: i for i, name in enumerate(optimizer.opt_params)}
    constraint1 = Constraint(
        func=constrain_bounds,
        name="bounds",
        ctx={"bounds": {param_positions[k]: v for k, v in bounds.items()}},
    )
    constraint2 = Constraint(
        func=kadap_constrain,
        name="kadap>=f(k2)",
        ctx={
            "idx_kadap": param_positions["k_adap"],
            "C_m": cell_params["C_m"],
            "tau_m": cell_params["tau_m"],
        },
    )
    optimizer.add_constraint(constraint1)
    optimizer.add_constraint(constraint2)
    optimizer.N_GEN = 300
    optimizer.init_params_fn = random_init
    optimizer.evaluate_fn = evaluate_GrC
    optimizer.print_evolution = True
    output_folder = 'GrC_opt'
    os.makedirs(output_folder, exist_ok=True)
    optimizer.output_path = output_folder
    toolbox = optimizer.set_optimizer()

    best_individual, best_fitness = optimizer.optimize()
    print("Optimization complete!")
    print("Best individual:", best_individual)
    print("Fitness:", best_fitness)
