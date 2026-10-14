import json
import os
import pickle

import numpy as np
from utils import *

from cerebellar_models.optimization.fitness import *


def gene_to_params(gene):
    gene = np.asarray(gene, dtype=float).reshape(-1)
    if gene.size != 5:
        raise ValueError(f"Gene must be length 6, got {gene.size}")

    return {
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
        "k_2": 1.0 / 44.0,
        "I_e": float(gene[0]),
        "k_adap": float(gene[1]),
        "k_1": float(gene[2]),
        # "k_2": float(gene[3]),
        "A1": float(gene[3]),
        "A2": float(gene[4]),
    }


def load_best_history(path_hist):
    """Load last record from best_history.pkl."""
    records = []
    with open(path_hist, "rb") as f:
        while True:
            try:
                records.append(pickle.load(f))
            except EOFError:
                break
    if not records:
        raise RuntimeError("best_history.pkl is empty or unreadable.")
    last = records[-1]
    best_params = np.asarray(last["params"], dtype=float).reshape(-1)
    best_fitness = np.asarray(last["fitness"], dtype=float).reshape(-1)
    return best_params, best_fitness, last


if __name__ == "__main__":

    cell_name = "GoC"

    path_genes = "neuron_opt/GoC_opt/pareto_snapshots/pareto_genes_300.npy"
    path_fit = "neuron_opt/GoC_opt/pareto_snapshots/pareto_fitness_300.npy"
    path_hist = "neuron_opt/GoC_opt/best_history.pkl"

    TOL_ABS = 0.05  # absolute tolerance for ALL stored objectives

    genes = np.load(path_genes)
    fitness = np.load(path_fit)
    best_params, best_fitness_ref, last_rec = load_best_history(path_hist)

    if genes.ndim != 2 or fitness.ndim != 2:
        raise ValueError(
            f"Expected 2D arrays. genes.ndim={genes.ndim}, fitness.ndim={fitness.ndim}"
        )
    if genes.shape[0] != fitness.shape[0]:
        raise ValueError(f"Mismatch N: genes {genes.shape}, fitness {fitness.shape}")
    if genes.shape[1] != best_params.size:
        raise ValueError(
            f"Mismatch D: genes has D={genes.shape[1]} but best_params has {best_params.size}"
        )
    if fitness.shape[1] != best_fitness_ref.size:
        raise ValueError(
            f"Mismatch M: fitness has M={fitness.shape[1]} but best_fitness has {best_fitness_ref.size}"
        )

    print("Best (from best_history.pkl)")
    print("  params  :", best_params)
    print("  fitness :", best_fitness_ref)

    protocol = {
        "start_stim": 200.0,
        "end_stim": 700.0,
        "duration": 2000.0,
        "threshold": -55,
    }

    data_dir = os.path.abspath(f"./tofit_eglif/results_tofitEglif/{cell_name}/")
    targ = extract_multicomp_features(multicomp_data=data_dir, protocol=protocol)

    delta = np.abs(fitness - best_fitness_ref[None, :])
    admissible_mask = np.all(delta <= TOL_ABS, axis=1)
    admissible_idxs = np.flatnonzero(admissible_mask)

    print("\nAdmissible set")
    print("  Num admissible:", admissible_idxs.size)

    if admissible_idxs.size == 0:
        print("No admissible solutions -> fallback to best_history params.")
        best_gene = best_params
        final_idx = None
    else:
        post_stim_pause = np.empty(admissible_idxs.size, dtype=float)

        for j, idx in enumerate(admissible_idxs):
            print("Processing gene: ", idx)
            gene = genes[idx]
            params = gene_to_params(gene)

            _, nest = extract_nest_features(
                multicomp_features=targ, cell_params=params, protocol=protocol
            )

            _, thr = rheobase_loss(nest, targ)
            post_stim_pause[j] = post_first_spike_loss(targ, nest, protocol, thr, mode="max")

        local_best_j = int(np.nanargmin(post_stim_pause))
        final_idx = int(admissible_idxs[local_best_j])
        best_gene = genes[final_idx]

        print(f"\nSelected admissible gene index: {final_idx}")
        print(f"Min post-first-spike (tie-break): {post_stim_pause[local_best_j]}")

        stored = fitness[final_idx].astype(float).reshape(-1)
        print("Final gene stored fitness:", stored)
        print(
            "Max abs diff vs best_fitness_ref (stored):",
            float(np.max(np.abs(stored - best_fitness_ref))),
        )
        print(
            "All within TOL_ABS (stored)?",
            bool(np.all(np.abs(stored - best_fitness_ref) <= TOL_ABS)),
        )

    cell_params_best = gene_to_params(best_gene)
    print("\nFinal best params:", cell_params_best)

    _, nest = extract_nest_features(
        multicomp_features=targ, cell_params=cell_params_best, protocol=protocol
    )

    pacemaking = pacemaking_loss(targ, nest)
    cv_err = cv_loss(targ, nest, regularity=True, alpha=0.1)
    rheo, thr = rheobase_loss(nest, targ)
    sel = currents_above_thr(targ, thr)
    curv = curvature_loss(targ, nest, sel, use_max_slope_penalty=True)
    gap = gap_loss(targ, nest, sel, weighted="gaussian")
    pos_loss = post_first_spike_loss(targ, nest, protocol, thr, mode="max")
    neg_loss = post_rebound_loss(targ, nest, protocol, thr=thr, sign="neg", window=100.0, scale=5)

    print(
        f"""Final errors:
    Pacemaking:       {pacemaking}
    CV:               {cv_err}
    Rheobase:         {rheo}
    Curvature:        {curv}
    Gap:              {gap}
    Post-first-spike: {pos_loss}
    Rebound:          {neg_loss}
    """
    )

    error_dict = {
        "pacemaking_error": float(pacemaking),
        "cv_error": float(cv_err),
        "rheobase_error": float(rheo),
        "curvature_error": float(curv),
        "gap_error": float(gap),
        "pos_error": float(pos_loss),
        "neg_error": float(neg_loss),
    }

    results = nest_protocol(targ, cell_params_best, protocol, multimeter=True)

    output_dir = f"./figures/{cell_name}"
    os.makedirs(output_dir, exist_ok=True)

    plot_results(
        targ,
        results,
        cell_params_best,
        protocol["start_stim"],
        protocol["end_stim"],
        cell_name,
        output_dir=output_dir,
    )

    os.makedirs("./results_opt", exist_ok=True)

    with open(f"./results_opt/{cell_name}_opt.json", "w") as f:
        json.dump(cell_params_best, f, indent=4)

    with open(f"./results_opt/{cell_name}_opt_err.json", "w") as f:
        json.dump(error_dict, f, indent=4)
