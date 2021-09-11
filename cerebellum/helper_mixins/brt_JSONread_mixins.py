from anytree.search import find, findall
from anytree import PreOrderIter, Node


class BrainNode(Node):
    separator = "|"

    def get_path_str(self):
        return "|" + "|".join([i.name for i in list(self.path)])


class JSONreadMixin(object):
    """Making possible previous functionalities"""

    def id_to_region_dictionary_ALLNAME(self):
        """Get path for each of the nodes ivolved.
        id_to_dictionary_ALLNAME"""

        return {str(node.id): node.get_path_str() for node in PreOrderIter(self.tree)}

    def id_to_region_dictionary(self):
        """Get path for each of the nodes ivolved.
        id_to_dictionary"""

        return {str(node.id): node.name for node in PreOrderIter(self.tree)}

    def get_id_gr_pc_mol(self):
        id_mol, id_pc, id_gr = -1, -1, -1
        for region_node in self.involved_regions:
            full_name = region_node.get_path_str()
            if self.region_name + "," in full_name:
                if "molecular layer" in full_name:
                    id_mol = region_node.id
                elif "granular layer" in full_name:
                    id_gr = region_node.id
                elif "Purkinje layer" in full_name:
                    id_pc = region_node.id
        return id_gr, id_pc, id_mol

    def get_id_region(self, region_string):
        for region_node in self.involved_regions:
            full_name = region_node.get_path_str()
            if region_string in full_name:
                return region_node.id

    def _region_of_interest_updates(self):
        if not self.region_name:
            raise ("Region name not set.")

        self.region_of_interest = find(
            self.tree, filter_=lambda node: node.name in (self.region_name)
        )

        # previously called id_region
        involved_regions = findall(
            self.region_of_interest, filter_=lambda node: self.region_of_interest in node.path
        )
        # issue #2
        # involved_regions = [
        #    i
        #    for i in involved_regions
        #    if i.id is not self.region_of_interest.id
        # ]  # remove Flocculus itself
        self.involved_regions = involved_regions
        self.id_region = [i.id for i in involved_regions]
