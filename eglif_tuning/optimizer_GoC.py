import logging
import os
import warnings

from optimizer_GrC import *

from cerebellar_models.optimization.features import *
from cerebellar_models.optimization.fitness import *
from cerebellar_models.optimization.optimizer import (
    ArchiveND,
    Constraint,
    Optimizer,
    _apply_constraints,
)

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)
for name in ("deap", "matplotlib", "numexpr"):
    logging.getLogger(name).setLevel(logging.ERROR)
np.seterr(all="ignore")
os.environ.setdefault("PYTHONWARNINGS", "ignore")

cell_params = {
    "t_ref": 2,
    "C_m": 145,
    "V_th": -55,
    "V_reset": -75,
    "E_L": -62,
    "V_min": -150,
    "V_m": -60.0,
    "lambda_0": 0.15,
    "tau_V": 0.09,
    "tau_m": 44,
    "I_e": 16.0,
    "k_adap": 0.2,
    "k_1": 0.03,
    "k_2": 1.0 / 44.0,
    "A1": 260.0,
    "A2": 180.0,
}

protocol = {
    "start_stim": 200.0,
    "end_stim": 700.0,
    "duration": 2000.0,
    "threshold": cell_params["V_th"],
}

bounds = {
    "I_e": (5.0, 50.0),
    "k_adap": (0.0, 0.50),
    "k_1": (0.0, 2.0),
    #"k_2": (1.0 / cell_params["tau_m"] - 1e-9,2.0),  # ensure oscillatory/damped oscillations regimes (k2>1/tau_m)
    "A1": (0.0, 500.0),
    "A2": (0.0, 500.0),
}

data_folder = os.path.join(os.path.dirname(__file__), "tofit_eglif/results_tofitEglif/GoC/")


# EVALUATION function
def evaluate_GoC(opt, ind, cell_params=None):
    params = dict(cell_params or opt.model_params)
    for i, p in enumerate(opt.opt_params):
        params[p] = float(ind[i])
    opt.model_params = params

    targ = opt._extract_multicomp_features()
    nest = opt._extract_nest_features()

    pacemaking = pacemaking_loss(targ, nest)
    cv_err = cv_loss(targ, nest, regularity=True, alpha=0.1)
    rheo, thr = rheobase_loss(nest, targ)
    sel = currents_above_thr(targ, thr)
    curv = curvature_loss(targ, nest, sel, use_max_slope_penalty=True)
    gap = gap_loss(targ, nest, sel, weighted='gaussian')
    #pos_loss = post_first_spike_loss(targ, nest, protocol, thr)
    neg_loss = post_rebound_loss(targ, nest, protocol, thr=thr, sign='neg', window=100.0, scale=5)

    return (
        float(pacemaking),
        float(cv_err),
        float(rheo),
        float(curv),
        float(gap),
        #float(pos_loss),
        float(neg_loss),
    )


if __name__ == "__main__":

    fitness = {
        k: -1
        for k in [
            "pacemaking_error",
            "cv_error",
            "rheobase_error",
            "curvature_error",
            "gap_error",
            #"pos_error",
            "neg_error",
        ]
    }

    archive = ArchiveND(eps_vec=[0.05] * 6, tau_vec=[0.15] * 6, cap=1000)

    opt = Optimizer(
        multicomp_data=data_folder,
        protocol=protocol,
        nest_model="eglif_multirec_opt",
        opt_params=list(bounds.keys()),
        model_params=cell_params,
        fitness=fitness,
        bounds=bounds,
        archive=archive,
        knee_weights=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] #1.0],
    )

    # --- Constraints ---
    opt.add_constraint(
        Constraint(
            func=constrain_bounds,
            name="bounds",
            ctx={"bounds": {i: v for i, v in enumerate(bounds.values())}},
        )
    )
    opt.add_constraint(
        Constraint(
            func=kadap_constrain,
            name="kadap>=f(k2)",
            ctx={
                "idx_kadap": opt.opt_params.index("k_adap"),
                "C_m": cell_params["C_m"],
                "tau_m": cell_params["tau_m"],
            },
        )
    )

    # --- Initialization & Evaluation ---
    opt.POP_SIZE = 600
    opt.NO_IMPROVE_PATIENCE = 10
    opt.ETA_MUTATION = 10.
    opt.ETA_CROSSOVER = 5.0
    opt.init_params_fn = random_init
    opt.evaluate_fn = evaluate_GoC
    opt.print_evolution = True
    opt.N_GEN = 300
    output_folder = 'GoC_opt'
    os.makedirs(output_folder, exist_ok=True)
    opt.output_path = output_folder
    opt.return_pareto = True
    opt.set_optimizer()

    best, fit, best_gen = opt.optimize()
    print("Optimization complete")
    print("Best individual:", best)
    print("Fitness:", fit)
