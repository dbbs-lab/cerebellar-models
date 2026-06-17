import time
import matplotlib.pyplot as plt
from cerebellar_models.MFM.transfer_functions import TransferFunctionSimulator, TransferFunction

if __name__ == "__main__":

    SIMULATE = False
    CELL_NAME     = 'granule_cell'
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
            n_workers=8,
        )
        sim.n_reps=50

        print(f'Summary parameters: {sim.sim_info}')
        print("Running TF simulation...")
        t_start = time.time()

        tf = sim.compute_tf()

        elapsed = time.time() - t_start
        print(f"\nRun duration: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    else:
        tf = TransferFunction.load(CELL_NAME, OUTPUT_DIR)

    tf.summary()
    tf.plot_tf()
    plt.show()