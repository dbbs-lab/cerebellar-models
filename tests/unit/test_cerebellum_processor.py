import unittest
import numpy as np
from anytree import findall, find
from cerebellum.cerebellum_processor import CerebellumProcessor
from cerebellum.brain_regions_tree import BrainRegionsTree


class TestBrainAtlasProcessor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brt = BrainRegionsTree(region_name="Flocculus")
        cls.bap = CerebellumProcessor(brt=cls.brt, nrrd_path="../Flocculus3.0_Lingula/data/")

    def test_pipeline(self):
        self.bap.mask_regions()
        self.bap.fill_regions()
        self.bap.show_regions()

    def test_init(self):
        # print("Annotations shape: " + str(self.bap.ann.shape))
        self.assertEqual(self.bap.ann.shape, self.bap.dens_cell.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.dens_neuron.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.dens_inh.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.orientations.shape[1:])

    def test_render_regions(self):
        # print(self.bap.render_regions())
        pass

    def test_leaf_or_next(self):
        """Is there some part of the annotations that is only
        non-leaf-region?"""

        # is our region (Flocculus) represented in ann?
        roi_id = self.bap.brt.region_of_interest.id
        sum_area = np.sum(self.bap.ann == roi_id)
        self.assertEqual(sum_area, 0)

        # is there any non-leaf (except root) in ann
        tree = self.bap.brt.tree
        base_node = find(tree, lambda node: node.name == "Basic cell groups and regions")
        non_leaf_nodes = findall(base_node, filter_=lambda node: node.children != [])
        passed = False
        for node in non_leaf_nodes:
            # skip root and basic region
            if node.id in [997, 8]:
                continue
            if np.sum(self.bap.ann == node.id) > 0:
                passed = True
                break
        if not passed:
            self.assertTrue(False)
