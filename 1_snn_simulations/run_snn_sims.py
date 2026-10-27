import os
from pathlib import Path

import numpy as np

from cerebellar_models.MFM.snn_for_mfm import SNNSimulate_forMFM

SIMULATE = False
COMPUTE_PLV = True

REPO_ROOT = Path(__file__).resolve().parent.parent
os.chdir(REPO_ROOT)

PLV_OUTPUT = "1_snn_simulations/plv.json"
PLV_SEED = 1234


if __name__ == "__main__":
    mf_stims = np.arange(start=4, stop=84, step=4)

    snn_sims_mfm = SNNSimulate_forMFM(
        mf_stims=mf_stims,
        config_path="circuit_Z-.yaml",
        scaffold_path="mouse_cerebellum.hdf5",
        tf_json="TF_dict.json",
        output_path="1_snn_simulations/snn_sims_mfm",
    )

    snn_sims_mfm.configure()
    snn_sims_mfm.prepare_sim()

    if SIMULATE:
        snn_sims_mfm.simulate()

    snn_sims_mfm.print_results(pdf_dir="1_snn_simulations/snn_sims_mfm_pdf")
    results = snn_sims_mfm.compute_firing_rates(save_spiketrains=False)
    #snn_sims_mfm.plot_raster()

    if COMPUTE_PLV:
        plv = snn_sims_mfm.compute_plv(output_path=PLV_OUTPUT, seed=PLV_SEED)
        print(f"PLV saved → {PLV_OUTPUT}")
