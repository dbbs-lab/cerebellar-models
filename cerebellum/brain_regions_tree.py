from anytree import RenderTree, AnyNode, Node, PreOrderIter
from anytree.importer import DictImporter
from anytree.exporter import DictExporter
from anytree.search import find, findall

import json
import nrrd
from cerebellum.JSONread import (
    DATA_PATH,
    region_name,
    id_to_region_dictionary,
    id_current_region,
    id_region,
)


class BrainNode(Node):
    separator = "|"


class BrainRegionsTree:
    def __init__(
        self,
        file="config/brain_regions.json",
        region_name="Flocculus",
    ):
        self.tree = self._create_tree(file, nodecls=BrainNode)
        self.region_name = region_name
        self._subtree()

    @staticmethod
    def _create_tree(brain_region_file=None, nodecls=AnyNode):
        if not brain_region_file:
            return None
        with open(brain_region_file) as f:
            brain_regions_dict = json.load(f)["msg"][0]

        importer = DictImporter(nodecls)
        brain_regions_tree = importer.import_(brain_regions_dict)
        return brain_regions_tree

    def _subtree(self):
        if not self.region_name:
            raise ("Region name not set.")

        self.region_of_interest = find(
            self.tree, filter_=lambda node: node.name in (self.region_name)
        )

        # previously called id_region
        id_region_nodes = findall(
            self.region_of_interest, filter_=lambda node: self.region_of_interest in node.path
        )
        id_region = [i.id for i in id_region_nodes]
        id_region.remove(self.region_of_interest.id)  # remove Flocculus itself
        self.id_region = id_region

    def as_dict(self):
        exporter = DictExporter()
        return exporter.export(self.tree)

    def all_nodes_dict(self):
        """Get path for each of the nodes ivolved."""

        return {
            str(node.id): "|" + "|".join([i.name for i in list(node.path)])
            for node in PreOrderIter(self.tree)
        }

    def render(self):
        return RenderTree(self.tree)


def get_mask_purkinje():
    ann, h = nrrd.read(DATA_PATH + "data/annotations.nrrd")
    region_names = []
    mask = {}
    for id_ in id_region:
        region_names.append(id_to_region_dictionary[id_])
        region, layer = id_to_region_dictionary[id_].split(", ")
        mask[layer] = ann == id_
    mask[region_name] = ann == id_current_region
    return mask
