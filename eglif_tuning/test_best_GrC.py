from utils import *

from cerebellar_models.optimization.fitness import *


def build_cell_params_GrC(best_params):
    # keep EXACT GrC-specific fixed params and mapping
    return {
        "t_ref": 1.5,
        "V_min": -150,
        "C_m": 7,
        "V_th": -41,
        "V_reset": -70,
        "E_L": -62,
        "V_m": -62.0,
        "tau_m": 24.15,
        "k_2": 1.0 / 24.15,  # oscillatory regime (fixed)
        "I_e": float(best_params[0]),
        "k_adap": float(best_params[1]),
        "k_1": float(best_params[2]),
        "A1": float(best_params[3]),
        "A2": float(best_params[4]),
    }


def compute_errors_GrC(target_features, nest_features):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)
    curv_error = curvature_loss(target_features, nest_features, sel, use_max_slope_penalty=True)
    gap_error = gap_loss(target_features, nest_features, sel, weighted="inverse")

    return {
        "Rheobase": float(rheobase_error),
        "FI curvature": float(curv_error),
        "FI gap": float(gap_error),
    }


if __name__ == "__main__":
    cell_name = "GrC"

    best_hist_path = "neuron_opt/GrC_opt/best_history.pkl"
    rec = load_best_record(best_hist_path)

    best_params = rec["params"]
    cell_params_best = build_cell_params_GrC(best_params)

    print(f"Best cell parameters for {cell_name}: {cell_params_best}")

    protocol = {
        "start_stim": 100.0,
        "end_stim": 600.0,
        "duration": 700.0,
        "threshold": cell_params_best["V_th"],
    }

    data_folder = f"tofit_eglif/results_tofitEglif/{cell_name}/"

    target_features = extract_multicomp_features(multicomp_data=data_folder, protocol=protocol)
    _, nest_features = extract_nest_features(
        multicomp_features=target_features,
        cell_params=cell_params_best,
        protocol=protocol,
    )

    error_dict = compute_errors_GrC(target_features, nest_features)
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
