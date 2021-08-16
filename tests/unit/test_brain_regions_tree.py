import unittest
import json

# import pytest

from anytree.search import findall
from cerebellum.brain_regions_tree import BrainRegionsTree


class TestRegionConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brain_regions_tree = BrainRegionsTree()
        cls.roi = cls.brain_regions_tree.region_of_interest
        with open("tests/fixtures/id_current_region.json") as f:
            cls.current_region = json.load(f)
        with open("tests/fixtures/id_region.json") as f:
            cls.id_region = json.load(f)
        with open("tests/fixtures/id_to_region_dictionary_ALLNAME.json") as f:
            cls.id_to_region_dictionary_ALLNAME = json.load(f)

    def test_id_current_region(self):
        self.assertEqual(self.brain_regions_tree.region_name, "Flocculus")
        self.assertEqual(self.roi.id, 1049)

    def test_id_region(self):
        # get nodes that are related to Flocculus
        id_region_nodes = findall(self.roi, filter_=lambda node: self.roi in node.path)
        id_region = [i.id for i in id_region_nodes]
        id_region.remove(self.roi.id)  # remove Flocculus itself
        self.assertEqual(id_region, self.id_region)
        self.assertEqual(self.brain_regions_tree.id_region, self.id_region)

    def test_id_to_region_dictionary_ALLNAME(self):
        self.assertEqual(
            self.brain_regions_tree.all_nodes_dict(), self.id_to_region_dictionary_ALLNAME
        )
