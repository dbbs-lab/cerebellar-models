import time
import matplotlib.pyplot as plt
from cerebellar_models.MFM.transfer_functions import TransferFunctionSimulator, TransferFunction
import numpy as np


if __name__ == "__main__":

    SIMULATE = False
    CELL_NAME     = 'purkinje_cell'
    SAVE_NAME     = 'purkinje_cell_Z-'
    YAML_PATH     = 'circuit_Z-.yaml'
    TF_DICT_PATH  = 'TF_dict.json'
    SCAFFOLD_PATH = 'mouse_cerebellum.hdf5'
    OUTPUT_DIR = 'output_TF'

    if SIMULATE:
        sim = TransferFunctionSimulator(
            scaffold_path=SCAFFOLD_PATH,
            yaml_path=YAML_PATH,
            tf_dict_path=TF_DICT_PATH,
            cell_name=CELL_NAME,
            n_workers=14,
        )

        sim.n_reps = 50

        print(f'Summary parameters: {sim.sim_info}')
        print("Running TF simulation...")
        t_start = time.time()

        tf = sim.compute_tf(save=False)
        tf.cell_name = SAVE_NAME
        tf.save(OUTPUT_DIR)

        elapsed = time.time() - t_start
        print(f"\nRun duration: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    else:
        tf = TransferFunction.load(cell_name=SAVE_NAME, output_dir=OUTPUT_DIR)

    # check firing rates
    grc_idx = np.argmin(np.abs(tf.freq_axes[tf.tags.index('GrC')] - 4))
    for mli_rate in [10, 20]:
        mli_idx = np.argmin(np.abs(tf.freq_axes[tf.tags.index('MLI')] - mli_rate))
        fr = tf.tf_mean[grc_idx, mli_idx]
        print(f"GrC=4 Hz, MLI={mli_rate} Hz -> PC firing rate = {fr:.2f} Hz")

    print(f"GrC=0 Hz, MLI=0 Hz -> PC firing rate = {tf.tf_mean[0, 0]:.2f} Hz\n")

    tf.summary()
    tf.plot_tf()
    plt.show()
