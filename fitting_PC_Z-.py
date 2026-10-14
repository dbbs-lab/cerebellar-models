import matplotlib.pyplot as plt

from cerebellar_models.MFM.analytical_tf import (
    AnalyticalTransferFunction,
    AnalyticalTransferFunctionFitting,
)
from cerebellar_models.MFM.numerical_tf import NumericalTransferFunction

if __name__ == "__main__":
    FIT = True

    # circuit cell key and numerical template name
    CELL_NAME = "purkinje_cell"
    TEMPLATE_NAME = "purkinje_cell_Z-"
    NUMERICAL_TFS_PATH = "./numerical_TFs/output_TFs"
    YAML_PATH = "./circuit_Z-.yaml"
    TF_DICT_PATH = "./TF_dict.json"
    SCAFFOLD_PATH = "./mouse_cerebellum.hdf5"
    OUTPUT_DIR = "./output_analytical_TFs"

    # Load numerical TF template (Tutorial 2 output)
    num_tf = NumericalTransferFunction.load(
        cell_name=TEMPLATE_NAME, output_dir=NUMERICAL_TFS_PATH
    )

    if FIT:
        fitter = AnalyticalTransferFunctionFitting(
            scaffold_path=SCAFFOLD_PATH,
            yaml_path=YAML_PATH,
            tf_dict_path=TF_DICT_PATH,
            cell_name=CELL_NAME,
            numerical_tf=num_tf,
            cut=0,
        )

        an_tf = fitter.optimize_alpha(output_dir=OUTPUT_DIR)
        fitter.plot_fit(output_dir=OUTPUT_DIR)
        fitter.plot_tf(output_dir=OUTPUT_DIR)
        print("Metrics:", an_tf.metrics(fitter.MEAN))
    else:
        an_tf = AnalyticalTransferFunction.load(
            cell_name=TEMPLATE_NAME, output_dir=OUTPUT_DIR
        )
        an_tf.plot_tf(numerical_tf=num_tf, output_dir=OUTPUT_DIR)

    an_tf.summary()
    plt.show()