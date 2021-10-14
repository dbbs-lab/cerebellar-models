from anytree import RenderTree
from anytree.importer import DictImporter
from anytree.exporter import DictExporter
from anytree.search import find, findall
from anytree import Node


class RegionNode(Node):
    separator = "|"

    def get_path_str(self):
        return "|" + "|".join([i.name for i in list(self.path)])

    def get_name(self):
        return [i.name for i in list(self.path)][-1]


class OntologyTree:
    """Takes a dict of regions based on TODO: find reference
    and provides the ids of the involved regions, and allows subtree.
    """

    def __init__(self, tree_dict):
        importer = DictImporter(RegionNode)
        self._tree = importer.import_(tree_dict)

    def _get_roi(self, region_name):
        roi = find(self._tree, filter_=lambda node: node.name in region_name)
        if roi is None:
            raise Exception("Region not found.")
        return roi

    @staticmethod
    def node_to_dict(node):
        if isinstance(node, OntologyTree):
            node = node._tree
        exporter = DictExporter()
        return exporter.export(node)

    def render_all(self):
        return RenderTree(self._tree)

    def subtree(self, region_name):
        roi = self._get_roi(region_name)
        json_tree = self.node_to_dict(roi)
        return type(self)(json_tree)

    def get_ids(self, region_name):
        roi = self._get_roi(region_name)

        # get all its children
        involved_regions = findall(roi, filter_=lambda node: roi in node.path)
        return [i.id for i in involved_regions]
