from anytree.importer import DictImporter
from anytree.exporter import DictExporter
from anytree.search import find, findall
from anytree import RenderTree, Node
from abc import ABC, abstractmethod


class OntologyTreeAbstract(ABC):
    """A clean wrapper to the ontology-based region mapping of the brain."""

    @abstractmethod
    def __init__(self, tree_dict):
        """Creates"""
        pass

    @abstractmethod
    def subtree(self, region_name):
        """Returns a subregion with similar properties of the base class."""
        pass

    @abstractmethod
    def get_ids(self, region_name):
        """Returns the ids related to the region name."""
        pass


class RegionNode(Node):
    separator = "|"

    def get_path_str(self):
        return "|" + "|".join([i.name for i in list(self.path)])

    def get_name(self):
        return [i.name for i in list(self.path)][-1]


class OntologyTree(OntologyTreeAbstract):
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
        involved_regions = findall(roi, filter_=lambda node: roi in node.path)
        return [i.id for i in involved_regions]
