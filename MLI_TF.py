import matplotlib.pyplot as plt
import numpy as np
from bsb import from_storage

from cerebellar_models.analysis import StructureReport
from cerebellar_models.MFM.transfer_functions import TransferFunction

SCAFFOLD_PATH = "mouse_cerebellum.hdf5"
OUTPUT_DIR = "output_TF"


if __name__ == "__main__":
    scaffold = from_storage(SCAFFOLD_PATH)
    report = StructureReport(scaffold)
    report.plots["density_table"].update()
    counts = report.plots["density_table"].get_counts()

    n_bc = counts["BC"]
    n_sc = counts["SC"]
    total = n_bc + n_sc
    w_bc = n_bc / total
    w_sc = n_sc / total
    print(f"basket_cell: {n_bc}  stellate_cell: {n_sc}  (w_BC={w_bc:.3f}, w_SC={w_sc:.3f})")

    tf_bc = TransferFunction.load("basket_cell", OUTPUT_DIR)
    tf_sc = TransferFunction.load("stellate_cell", OUTPUT_DIR)

    if tf_bc.tags != tf_sc.tags:
        raise ValueError(f"TF axes mismatch: {tf_bc.tags} vs {tf_sc.tags}")
    for ax_b, ax_s, tag in zip(tf_bc.freq_axes, tf_sc.freq_axes, tf_bc.tags):
        if not np.allclose(ax_b, ax_s):
            raise ValueError(f"freq_axis mismatch for label '{tag}'")

    tf_mean = w_bc * tf_bc.tf_mean + w_sc * tf_sc.tf_mean
    tf_std = np.sqrt(w_bc**2 * tf_bc.tf_std**2 + w_sc**2 * tf_sc.tf_std**2)

    tf_mli = TransferFunction("MLI", tf_mean, tf_std, tf_bc.tags, tf_bc.freq_axes)
    tf_mli.save(OUTPUT_DIR)

    tf_mli.summary()
    tf_mli.plot_tf()
    plt.show()
