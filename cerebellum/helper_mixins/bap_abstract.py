from abc import ABC, abstractmethod


class BrainAtlasProcessor(ABC):
    @abstractmethod
    def __init__(self, brt, nrrd_path):
        """Initializes region of interest, subregions, and nrrd data"""
        pass

    @abstractmethod
    def create_mask(self):
        """Create a mask for each of the involved regions."""
        pass

    @abstractmethod
    def apply_mask(self):
        """Create a single 3D array containing all the masks."""
        pass

    @abstractmethod
    def show_regions(self):
        """Show a slice of the regions involved."""
        pass
