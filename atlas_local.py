from cerebellum.brain_regions_tree import BrainRegionsTree
from cerebellum.brain_atlas_processor import BrainAtlasProcessor


brain_regions_tree = BrainRegionsTree(file="config/brain_regions.json", region_name="Flocculus")
bap = BrainAtlasProcessor(brt=brain_regions_tree)
