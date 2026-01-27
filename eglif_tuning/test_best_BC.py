from utils import *
from cerebellar_models.optimization.fitness import *

def build_cell_params_BC(best_params):
    return {
        "t_ref": 1.59,
        "C_m": 14.6,
        "V_th": -53,
        "V_reset": -78,
        "E_L": -68,
        "V_m": -60.0,
        "tau_m": 9.125,
        "I_e": float(best_params[0]),
        "k_adap": float(best_params[1]),
        "k_1": float(best_params[2]),
        "A1": float(best_params[3]),
        "A2": float(best_params[4]),
        "k_2": float(best_params[5]),
    }

def compute_errors_BC(target_features, nest_features, protocol):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)

    gap = gap_loss(target_features, nest_features, sel, weighted="gaussian")
    curv = curvature_loss(target_features, nest_features, sel, use_max_slope_penalty=True)
    pacemaking_error = pacemaking_loss(target_features, nest_features)
    cv_error = cv_loss(target_features, nest_features, regularity=True, alpha=0.1)
    pos_loss = post_first_spike_loss(target_features, nest_features, protocol, thr)
    neg_loss = poststim_autorhythm_loss(target_features, nest_features, protocol, thr)

    return {
        "Rheobase": float(rheobase_error),
        "FI curvature": float(curv),
        "FI gap": float(gap),
        "Coefficient of Variation": float(cv_error),
        "Pacemaking": float(pacemaking_error),
        "Pacemaking after pos stim": float(pos_loss),
        "Rebound after neg stim": float(neg_loss),
    }


if __name__ == "__main__":
    cell_name = "BC"

    best_hist_path = "neuron_opt/BC_opt/best_history.pkl"
    rec = load_best_record(best_hist_path)

    best_params = rec["params"]
    cell_params_best = build_cell_params_BC(best_params)
    print(f"Best cell parameters for {cell_name}: {cell_params_best}")

    protocol = {
        "start_stim": 500.0,
        "end_stim": 1500.0,
        "duration": 2500.0,
        "threshold": cell_params_best["V_th"],
    }

    data_folder = f"tofit_eglif/results_tofitEglif/{cell_name}/"

    target_features = extract_multicomp_features(multicomp_data=data_folder, protocol=protocol)
    _, nest_features = extract_nest_features(
        multicomp_features=target_features,
        cell_params=cell_params_best,
        protocol=protocol,
    )

    error_dict = compute_errors_BC(target_features, nest_features, protocol)
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
