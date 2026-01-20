import json
import os
import numpy as np
from cerebellar_models.optimization.fitness import *
from utils import *


def load_cell_params_from_json(cell_name: str, in_dir: str = "./results_opt") -> dict:
    path = os.path.join(in_dir, f"{cell_name}_opt.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing params json: {path}")
    with open(path, "r") as f:
        params = json.load(f)
    return params


def compute_errors(target_features, nest_features, protocol):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)

    gap = gap_loss(target_features, nest_features, sel, weighted="gaussian")
    curv = curvature_loss(target_features, nest_features, sel, use_max_slope_penalty=True)
    pacemaking_error = pacemaking_loss(target_features, nest_features)
    cv_error = cv_loss(target_features, nest_features, regularity=True, alpha=0.1)
    pos_loss = post_first_spike_loss(target_features, nest_features, protocol, thr)
    neg_loss = post_rebound_loss(target_features, nest_features, protocol, thr=thr, sign="neg", window=100.0, scale=5)

    return {
         "Rheobase": float(rheobase_error),
         "FI curvature": float(curv),
         "FI gap": float(gap),
         "Coefficient of Variation": float(cv_error),
         "Pacemaking": float(pacemaking_error),
         "Pause after pos stim": float(pos_loss),
         "Rebound after neg stim": float(neg_loss),
    }


def save_outputs(cell_name: str, cell_params_best: dict, error_dict: dict, out_dir: str = "./results_opt"):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, f"{cell_name}_opt.json"), "w") as f:
        json.dump(cell_params_best, f, indent=4)
    with open(os.path.join(out_dir, f"{cell_name}_opt_err.json"), "w") as f:
        json.dump(error_dict, f, indent=4)


if __name__ == "__main__":
    cell_name = "GoC"

    #load best params from json
    cell_params_best = load_cell_params_from_json(cell_name, in_dir="./results_opt")
    print(f"Best cell parameters for {cell_name}: {cell_params_best}")

    protocol = {
         "start_stim": 200.0,
         "end_stim": 700.0,
         "duration": 2000.0,
         "threshold": float(cell_params_best["V_th"]),
    }

    data_folder = f"tofit_eglif/results_tofitEglif/{cell_name}/"

    target_features = extract_multicomp_features(multicomp_data=data_folder, protocol=protocol)
    _, nest_features = extract_nest_features(
         multicomp_features=target_features,
         cell_params=cell_params_best,
         protocol=protocol,
    )

    error_dict = compute_errors(target_features, nest_features, protocol)
    print(f"Final errors: {error_dict}")

    # plots + traces
    output_fig_dir = f"./figures/{cell_name}"
    os.makedirs(output_fig_dir, exist_ok=True)

    results = nest_protocol(target_features, cell_params_best, protocol, multimeter=True)
    plot_results(
         target_features,
         results,
         cell_params_best,
         protocol["start_stim"],
         protocol["end_stim"],
         cell_name,
         output_dir=output_fig_dir,
    )

    # error plot (assumes you have this helper in utils.py)
    plot_error_bar(error_dict, cell_name, output_fig_dir)

    # save JSON outputs (same format/paths as before)
    save_outputs(cell_name, cell_params_best, error_dict, out_dir="./results_opt")
