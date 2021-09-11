import unittest
import json

from anytree.search import findall
from cerebellum.brain_regions_tree import BrainRegionsTree


class TestRegionConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brt = BrainRegionsTree(region_name="Flocculus")
        cls.roi = cls.brt.region_of_interest
        with open("tests/fixtures/id_current_region.json") as f:
            cls.current_region = json.load(f)
        with open("tests/fixtures/id_region.json") as f:
            cls.id_region = json.load(f)
        with open("tests/fixtures/id_to_region_dictionary_ALLNAME.json") as f:
            cls.id_to_region_dictionary_ALLNAME = json.load(f)

    def test_id_current_region(self):
        self.assertEqual(self.brt.region_name, "Flocculus")
        self.assertEqual(self.roi.id, 1049)

    def test_id_region(self):
        # get nodes that are related to Flocculus
        id_region_nodes = findall(self.roi, filter_=lambda node: self.roi in node.path)
        id_region = [i.id for i in id_region_nodes if i != self.roi] + [self.roi.id]
        # #2 Issue
        # id_region.remove(self.roi.id)  # remove Flocculus
        self.assertEqual(id_region, self.id_region)
        self.assertEqual(self.brt.id_region, self.id_region)

    def test_id_to_region_dictionary_ALLNAME(self):
        self.assertEqual(
            self.brt.id_to_region_dictionary_ALLNAME(),
            self.id_to_region_dictionary_ALLNAME,
        )

    def test_get_involved_regions_id(self):
        self.assertEqual(self.brt.get_id_gr_pc_mol(), (10690, 10691, 10692))

    def test_basic_properties(self):
        self.assertEqual(self.brt.region_of_interest in self.brt.involved_regions, True)
