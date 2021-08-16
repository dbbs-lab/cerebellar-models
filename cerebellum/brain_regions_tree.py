from anytree import RenderTree, AnyNode
from anytree.importer import DictImporter
from anytree.exporter import DictExporter
import json

from cerebellum.helper_mixins.brt_JSONread_mixins import (
    JSONreadMixin,
    BrainNode,
)


class BrainRegionsTree(JSONreadMixin):
    def __init__(
        self,
        file="config/brain_regions.json",
        region_name="Flocculus",
    ):
        self.tree = self._create_tree(file, nodecls=BrainNode)
        self.region_name = region_name
        self._region_of_interest_updates()

    @staticmethod
    def _create_tree(brain_region_file=None, nodecls=AnyNode):
        if not brain_region_file:
            return None
        with open(brain_region_file) as f:
            brain_regions_dict = json.load(f)["msg"][0]

        importer = DictImporter(nodecls)
        brain_regions_tree = importer.import_(brain_regions_dict)
        return brain_regions_tree

    def as_dict(self):
        exporter = DictExporter()
        return exporter.export(self.tree)

    def render(self):
        return RenderTree(self.tree)
