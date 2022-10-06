from imp import get_tag
import unittest
from bsb.unittest import NumpyTestCase
from bsb.config import from_json
from bsb.core import Scaffold
from bsb.storage import Storage
import numpy as np

#RandomStorageFixture, NumpyTestCase, 
class TestGranularLayer(unittest.TestCase, NumpyTestCase):
    def setUp(self):
        self.storage = Storage(engine="hdf5",root="granular")
        #self.storage = self.random_storage(engine="hdf5")
        self.cfg = from_json("cerebellum.json")
        self.cfg.network.x = 100
        self.cfg.network.y = 400
        self.cfg.network.z = 100
        self.network = Scaffold(self.cfg, self.storage)
        #self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","purkinje_layer_placement","molecular_layer_placement","mossy_fibers_to_glomerulus","glomerulus_to_golgi","glomerulus_to_granule","golgi_to_glomerulus","gap_goc","ascending_axon_to_golgi","parallel_fiber_to_golgi"])
        self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","mossy_fibers_to_glomerulus","glomerulus_to_granule","golgi_to_glomerulus"])



    
    def test_convergence(self):
        pass
