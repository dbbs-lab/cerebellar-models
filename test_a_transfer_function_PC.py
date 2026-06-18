import os

import matplotlib.pyplot as plt
import numpy as np

from cerebellar_models.MFM.a_transfer_functions import TransferFunctionFitting


if __name__ == "__main__":
    CELL_NAME = "granule_cell"
    YAML_PATH = "circuit_Z+.yaml"
    TF_DICT_PATH = "TF_dict.json"
    SCAFFOLD_PATH = "mouse_cerebellum.hdf5"
    OUTPUT_DIR = "output_TF"
    NUM_TF_FILE = os.path.join(OUTPUT_DIR, "TF_granule_cell.npz")
    print(NUM_TF_FILE)

    # Instantiate the fitter using saved numerical TF and model parameters.
    fitter = TransferFunctionFitting(
        scaffold_path=SCAFFOLD_PATH,
        yaml_path=YAML_PATH,
        tf_dict_path=TF_DICT_PATH,
        cell_name=CELL_NAME,
        num_tf=None,
    )

    # Basic sanity checks on loaded numerical TF.
    assert fitter.num_tf is not None, "num_tf was not loaded"
    assert fitter.num_tf.size > 0, "num_tf is empty"
    assert isinstance(fitter.num_tf, np.ndarray), "num_tf must be a numpy array"

    # Check that the expected saved file exists for the PC test.
    assert os.path.exists(NUM_TF_FILE), f"Numerical TF file not found: {NUM_TF_FILE}"

    # Validate that E_rev mapping was populated from the YAML constants.
    cell_params = fitter.sim_info["cell_params"]
    assert "Erev" in cell_params, "Erev mapping missing from sim_info.cell_params"
    assert isinstance(cell_params["Erev"], dict), "Erev should be a receptor-type dict"
    assert len(cell_params["Erev"]) in (2, 3), "Purkinje cell should have 2 or 3 reversal potentials"

    # Validate that the loaded E_rev mapping is consistent with connections.
    receptor_types = {
        int(info["receptor_type"])
        for info in fitter.sim_info["connections"].values()
    }
    missing_receptors = receptor_types - set(cell_params["Erev"].keys())
    assert not missing_receptors, (
        f"Connections reference unmapped receptor types: {missing_receptors}"
    )

    # Fit the analytical transfer function and report parameters.
    fitted_params = fitter.fit_tf(alpha=1.0, XX=0.0, maxiter=1000, xtol=1e-5)
    print(f"Fitted parameters: {fitted_params}")

    # Create a quick plot of the fitted analytical TF.
    #fitter.plot_tf(save=False)
    #plt.show()
