"""
Implementation of the BSB framework for cerebellar cortex reconstructions and simulations.
"""

import os

__version__ = "0.3.0"


def templates():  # pragma: nocover
    """
    :meta private:
    """
    return [os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))]


classmap = {
    "bsb.connectivity.strategy.ConnectionStrategy": {
        "mossy_glom": "cerebellum.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
        "glom_gran": "cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule",
        "golgi_glom": "cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus",
        "glom_golgi": "cerebellum.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi",
        "ubc_glom": "cerebellum.connectome.to_glomerulus.ConnectomeUBCGlomerulus",
        "glom_ubc": "cerebellum.connectome.glomerulus_ubc.ConnectomeGlomerulusUBC",
        "grc_to": "cerebellum.connectome.granule_voxel_intersect.GranuleToMorphologyIntersection",
    },
    "bsb.postprocessing.AfterConnectivityHook": {
        "struct_report": "cerebellum.analysis.structure_analysis.RunStructureReport",
    },
}
