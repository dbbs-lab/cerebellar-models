import json
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

from cerebellar_models.optimization.fitness import *
from utils import *


def load_best_record(pkl_path: str) -> dict:
    records = []
    with open(pkl_path, "rb") as f:
        while True:
            try:
                records.append(pickle.load(f))
            except EOFError:
                break
    if len(records) == 0:
        raise RuntimeError(f"Empty best_history file: {pkl_path}")
    return records[-1]


def build_cell_params_PC(best_params):
    # keep EXACT PC-specific fixed params and mapping you used
    return {
        "A1": float(best_params[3]),
        "A2": float(best_params[4]),
        "C_m": 334,
        "E_L": -59,
        "I_e": float(best_params[0]),
        "V_m": -59.0,
        "V_min": -350,
        "V_reset": -69,
        "V_th": -43,
        "k_1": float(best_params[5]),
        "k_2": float(best_params[2]),
        "k_adap": float(best_params[1]),
        "t_ref": 0.5,
        "tau_m": 47,
    }


def build_protocol_PC(cell_params_best):
    return {
        "start_stim": 200.0,
        "end_stim": 700.0,
        "duration": 2000.0,
        "threshold": cell_params_best["V_th"],
    }


def compute_errors_PC(target_features, nest_features, protocol):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)

    slope_error = slope_loss(target_features, nest_features, thr, sel)
    gap_error = gap_loss(target_features, nest_features, sel, weighted="inverse")
    pacemaking_error = pacemaking_loss(target_features, nest_features)
    cv_error = cv_loss(target_features, nest_features)
    pos_loss = post_first_spike_loss(target_features, nest_features, protocol, thr=thr)
    neg_loss = post_rebound_loss(target_features, nest_features, protocol, thr=thr, sign="neg", window=100.0)    # dati per ora a 50. window

    return {
        "Rheobase": float(rheobase_error),
        "FI gap error": float(slope_error),
        "gap_error": float(gap_error),
        "cv_error": float(cv_error),
        "pacemaking_error": float(pacemaking_error),
        "post_stim_pos_error": float(pos_loss),
        "post_stim_neg_error": float(neg_loss),
    }


def save_outputs(cell_name, cell_params_best, error_dict, out_dir="./results_opt"):
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/{cell_name}_opt.json", "w") as f:
        json.dump(cell_params_best, f, indent=4)
    with open(f"{out_dir}/{cell_name}_opt_err.json", "w") as f:
        json.dump(error_dict, f, indent=4)




if __name__ == "__main__":
    cell_name = "PC"

    best_hist_path = "neuron_opt/PC_opt/best_history.pkl"
    rec = load_best_record(best_hist_path)

    best_params = rec["params"]
    cell_params_best = build_cell_params_PC(best_params)
    print(f"Best cell parameters for {cell_name}: {cell_params_best}")

    protocol = build_protocol_PC(cell_params_best)
    data_folder = f"tofit_eglif/results_tofitEglif/{cell_name}/"

    target_features = extract_multicomp_features(multicomp_data=data_folder, protocol=protocol)
    _, nest_features = extract_nest_features(
        multicomp_features=target_features,
        cell_params=cell_params_best,
        protocol=protocol,
    )

    error_dict = compute_errors_PC(target_features, nest_features, protocol)
    print(f"Final errors: {error_dict}")

    # traces + plots
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

    # error plot
    plot_error_bar(error_dict, cell_name, output_fig_dir)

    # save JSON outputs in the SAME format/paths as before
    save_outputs(cell_name, cell_params_best, error_dict, out_dir="./results_opt")
