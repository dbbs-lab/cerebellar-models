import unittest
from cerebellum.brain_atlas_processor import BrainAtlasProcessor
from cerebellum.brain_regions_tree import BrainRegionsTree


class TestBrainAtlasProcessor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brt = BrainRegionsTree(region_name="Flocculus")
        cls.bap = BrainAtlasProcessor(brt=cls.brt)

    def test_pipeline(self):
        self.bap.mask_regions()
        self.bap.fill_regions()
        self.bap.show_regions()
