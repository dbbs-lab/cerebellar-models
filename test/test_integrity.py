import unittest
import bsb.options
from bsb.core import Scaffold
from bsb.config import from_json
from bsb.morphologies import Morphology
import pathlib
import glob
import os


@unittest.skipIf(not os.getenv("CEREB_INTEGRATION_TESTING", False), "Integration testing not enabled")
class TestCerebellumIntegrity(unittest.TestCase):
    """
    Check if the placement of all cell types is correct
    """
    @classmethod
    def setUpClass(cls):
        bsb.options.verbosity = 4
        cfg = from_json("cerebellum.json")
        cfg.network.x = 100
        cfg.network.z = 100
        cls.small = net = Scaffold(cfg)
        path = pathlib.Path(__file__).parent.parent / "morphologies"
        for f in glob.glob(str(path / "*.swc")):
            m = Morphology.from_swc(f)
            net.morphologies.save(f.split("/")[-1].split(".")[0], m, overwrite=True)
        cls.small.compile(clear=True)

    def test_small_network(self):
        net = self.small
