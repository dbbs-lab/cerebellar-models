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
    "A1": 157.622,
    "A2": 172.622,
    "C_m": 334,
    "E_L": -59,
    "I_e": 590.0,
    "V_m": -59.0,
    "V_min": -350,
    "V_reset": -69,
    "V_th": -43,
    "k_1": 0.195,
    "k_2": 0.041,
    "k_adap": 1.491,
    "t_ref": 0.5,
    "tau_m": 47,
}

protocol = {
    "start_stim": 200.0,
    "end_stim": 700.0,
    "duration": 2000.0,
    "threshold": cell_params["V_th"],
}

bounds = {
    "I_e": (0.01, 1500.0),
    "k_adap": (0.0, 0.50),
    "k_2": (
        1.0 / cell_params["tau_m"] - 1e-9,
        2.0,
    ),  # ensure oscillatory/damped oscillations regimes (k2>1/tau_m)
    "A1": (0.01, 2000.0),
    "A2": (0.01, 1500.0),
    "k_1": (0.0, 2.0),
}

data_folder = os.path.join(os.path.dirname(__file__), "tofit_eglif/results_tofitEglif/PC/")


# EVALUATION function
def evaluate_PC(opt, ind, cell_params=None):
    params = dict(cell_params or opt.model_params)
    for i, p in enumerate(opt.opt_params):
        params[p] = float(ind[i])
    opt.model_params = params

    targ = opt._extract_multicomp_features()
    nest = opt._extract_nest_features()

    pacemaking = pacemaking_loss(targ, nest)
    cv_err = cv_loss(targ, nest)
    rheo, thr = rheobase_loss(nest, targ)
    sel = currents_above_thr(targ, thr)
    #slope = slope_loss(targ, nest, thr, sel)
    curv = curvature_loss(targ, nest, sel, use_max_slope_penalty=True)
    gap = gap_loss(targ, nest, sel, weighted='gaussian')
    pos_loss = post_first_spike_loss(targ, nest, protocol, thr=thr)
    neg_loss = post_rebound_loss(targ, nest, protocol, thr=thr, sign='neg', window=70.0, scale = 2)

    return (
        float(pacemaking),
        float(cv_err),
        float(rheo),
        float(curv),
        float(gap),
        float(pos_loss),
        float(neg_loss),
    )


if __name__ == "__main__":

    fitness = {
        k: -1
        for k in [
            "pacemaking_error",
            "cv_error",
            "rheobase_error",
            "curv_error",
            "gap_error",
            "pos_lat_error",
            "neg_lat_error",
        ]
    }

    archive = ArchiveND(eps_vec=[0.05] * 7, tau_vec=[0.15] * 7, cap=1000)

    opt = Optimizer(
        multicomp_data=data_folder,
        protocol=protocol,
        nest_model="eglif_multirec_opt",
        opt_params=list(bounds.keys()),
        model_params=cell_params,
        fitness=fitness,
        bounds=bounds,
        archive=archive,
        knee_weights=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
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
    opt.init_params_fn = random_init
    opt.evaluate_fn = evaluate_PC
    opt.print_evolution = True
    opt.N_GEN = 500
    output_folder = 'PC_opt'
    os.makedirs(output_folder, exist_ok=True)
    opt.output_path = output_folder
    opt.set_optimizer()

    best, fit = opt.optimize()
    print("Optimization complete")
    print("Best individual:", best)
    print("Fitness:", fit)
