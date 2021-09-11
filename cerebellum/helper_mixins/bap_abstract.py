from abc import ABC, abstractmethod


class BrainAtlasProcessor(ABC):
    @abstractmethod
    def __init__(self, brt, nrrd_path):
        """Initializes region of interest, subregions, and nrrd data"""
        pass

    def mask_regions(self):
        """Create a mask for each of the involved regions."""
        pass

    def fill_regions(self):
        """Create a single 3D array containing all the masks."""
        pass

    def show_regions(self):
        """Show a slice of the regions involved."""
        pass
