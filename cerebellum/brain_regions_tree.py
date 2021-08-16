from anytree.importer import DictImporter
from anytree.exporter import DictExporter
import json


def create_tree(brain_region_file=None):
    if not brain_region_file:
        return None
    with open(brain_region_file) as f:
        brain_regions_dict = json.load(f)["msg"][0]

    importer = DictImporter()
    brain_regions_tree = importer.import_(brain_regions_dict)
    return brain_regions_tree


def export_tree_dict(brain_region_tree):
    exporter = DictExporter()
    return exporter.export(brain_region_tree)
