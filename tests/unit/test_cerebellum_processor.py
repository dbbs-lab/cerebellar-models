import unittest
import numpy as np
from anytree import findall, find
from cerebellum.cerebellum_processor import CerebellumProcessor
from cerebellum.brain_regions_tree import BrainRegionsTree
import pytest


class TestBrainAtlasProcessor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brt = BrainRegionsTree(region_name="Lingula (I)")
        cls.bap = CerebellumProcessor(brt=cls.brt, nrrd_path="../Flocculus3.0_Lingula/data/")

    @pytest.mark.skip
    def test_read_processor(self):
        print()
        print("Testing read/load. Might take some time.")
        # go through all processes
        self.bap.run_pipeline()

        # read ready made
        bap = CerebellumProcessor.read_processor(region_name="Lingula (I)")
        self.assertEqual(type(bap), CerebellumProcessor)

        # check properties needed if they are maintained
        valid = np.array(bap.ann) == np.array(self.bap.ann)
        self.assertTrue(valid.all())

        print("MASK")
        for key in bap.mask.keys():
            valid = np.array(bap.mask[key]) == np.array(self.bap.mask[key])
            self.assertTrue(valid.all())

        print("DENS NEURON")
        valid = np.array(bap.dens_neuron) == np.array(self.bap.dens_neuron)
        self.assertTrue(valid.all())

        print("DENS INH")
        valid = np.array(bap.dens_inh) == np.array(self.bap.dens_inh)
        self.assertTrue(valid.all())

        print("DENS ORIENTATIONS")
        self.assertEqual(len(bap.orientations), len(self.bap.orientations))
        for i in range(len(bap.orientations)):
            for j in range(len(bap.orientations[i])):
                bap_orient_saved = np.array(bap.orientations[i][j])
                bap_orient = np.array(self.bap.orientations[i][j])
                self.assertTrue(np.allclose(bap_orient_saved, bap_orient, equal_nan=True))

    @pytest.mark.skip
    def test_pipeline(self):
        self.bap.mask_regions()
        self.bap.fill_regions()
        self.bap.show_regions()
        self.bap.save_processor()

    @pytest.mark.skip
    def test_init(self):
        # print("Annotations shape: " + str(self.bap.ann.shape))
        self.assertEqual(self.bap.ann.shape, self.bap.dens_cell.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.dens_neuron.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.dens_inh.shape)
        self.assertEqual(self.bap.ann.shape, self.bap.orientations.shape[1:])

    @pytest.mark.skip
    def test_render_regions(self):
        # print(self.bap.render_regions())
        pass

    @pytest.mark.skip
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
