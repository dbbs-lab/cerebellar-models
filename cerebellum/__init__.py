"""
Implementation of the BSB framework for cerebellar cortex reconstructions and simulations.
"""

import os

__version__ = "0.2.0"


def templates():
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
    }
}
