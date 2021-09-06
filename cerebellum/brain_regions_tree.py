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
        region_name,
        brt_file="config/brain_regions.json",
    ):
        self.tree = self._create_tree(brt_file, nodecls=BrainNode)
        self.region_name = region_name
        self.region_of_interest = None
        self.involved_regions = []
        self.id_region = []
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

    def render_all(self):
        return RenderTree(self.tree)

    def render_roi(self):
        return RenderTree(self.region_of_interest)
