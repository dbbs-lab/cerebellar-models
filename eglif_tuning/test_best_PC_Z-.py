from utils import *

from cerebellar_models.optimization.fitness import *


def build_cell_params_PC(best_params):
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


def compute_errors_PC(target_features, nest_features, protocol):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)

    pacemaking = pacemaking_loss(target_features, nest_features)
    cv_err = cv_loss(target_features, nest_features)
    rheo, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)
    curv = curvature_loss(target_features, nest_features, sel, use_max_slope_penalty=True)
    gap = gap_loss(target_features, nest_features, sel)
    adap = adaptation_index_loss(target_features, nest_features, protocol, sel=sel)
    neg_loss = post_rebound_loss(
        target_features, nest_features, protocol, thr=thr, sign="neg", window=150.0, scale=7
    )

    return {
        "Rheobase": float(rheobase_error),
        "fI curvature error": float(curv),
        "fI gap_error": float(gap),
        "cv_error": float(cv_err),
        "pacemaking_error": float(pacemaking),
        "adaptation_index_error": float(adap),
        "post_stim_neg_error": float(neg_loss),
    }


if __name__ == "__main__":
    cell_name = "PC_Z-"

    best_hist_path = "neuron_opt/PC_Z-_opt/best_history.pkl"
    rec = load_best_record(best_hist_path)

    best_params = rec["params"]
    cell_params_best = build_cell_params_PC(best_params)
    print(f"Best cell parameters for {cell_name}: {cell_params_best}")

    protocol = {
        "start_stim": 1000.5,
        "end_stim": 2000.0,
        "duration": 3000.0,
        "threshold": cell_params_best["V_th"],
    }
    data_folder = f"tofit_eglif/results_tofitEglif/{cell_name}/"

    target_features = extract_multicomp_features(multicomp_data=data_folder, protocol=protocol)
    _, nest_features = extract_nest_features(
        multicomp_features=target_features,
        cell_params=cell_params_best,
        protocol=protocol,
    )

    error_dict = compute_errors_PC(target_features, nest_features, protocol)
    print(f"Final errors: {error_dict}")

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
    plot_error_bar(error_dict, cell_name, output_fig_dir)
    save_outputs(cell_name, cell_params_best, error_dict, out_dir="./results_opt")
