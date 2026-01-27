from utils import *
from cerebellar_models.optimization.fitness import *

def compute_errors_SC(target_features, nest_features, protocol):
    rheobase_error, thr = rheobase_loss(nest_features, target_features)
    sel = currents_above_thr(target_features, thr)

    gap = gap_loss(target_features, nest_features, sel, weighted="gaussian")
    curv = curvature_loss(target_features, nest_features, sel, use_max_slope_penalty=True)
    pacemaking_error = pacemaking_loss(target_features, nest_features)
    cv_error = cv_loss(target_features, nest_features, regularity=True, alpha=0.1)
    pos_loss = post_first_spike_loss(target_features, nest_features, protocol, thr)
    neg_loss = post_rebound_loss(
        target_features, nest_features, protocol, thr=thr, sign="neg", window=70.0
    )

    return {
        "Rheobase": float(rheobase_error),
        "FI curvature": float(curv),
        "FI gap": float(gap),
        "Coefficient of Variation": float(cv_error),
        "Pacemaking": float(pacemaking_error),
        "Pause after pos stim": float(pos_loss),
        "Rebound after neg stim": float(neg_loss),
    }


if __name__ == "__main__":
    cell_name = "SC"

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

    error_dict = compute_errors_SC(target_features, nest_features, protocol)
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
