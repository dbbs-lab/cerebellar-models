"""
Implementation of the BSB framework for cerebellar cortex reconstructions and simulations.
"""

import os

__version__ = "0.6.0"


def templates():  # pragma: nocover
    """
    :meta private:
    """
    return [os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))]


classmap = {
    "bsb.connectivity.strategy.ConnectionStrategy": {
        "mossy_glom": "cerebellar_models.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
        "glom_gran": "cerebellar_models.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
        "golgi_glom": "cerebellar_models.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
        "glom_golgi": "cerebellar_models.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi",
        "ubc_glom": "cerebellar_models.connectome.to_glomerulus.ConnectomeUBCGlomerulus",
        "glom_ubc": "cerebellar_models.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
        "io_mli": "cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
    },
    "bsb.postprocessing.AfterPlacementHook": {
        "label_microzones": "cerebellar_models.placement.microzones.LabelMicrozones",
    },
    "bsb.postprocessing.AfterConnectivityHook": {
        "struct_report": "cerebellar_models.analysis.structure_analysis.RunStructureReport",
    },
}
