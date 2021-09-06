import nrrd

PLOTTING = True
VOXEL_SIZE = 25.0  # um
FAC_LUGARO = 15.0  # Lugaro Cell fraction/purkinje cells


reference_UBC_density = {
    "Lingula (I)": 6000,
    "Uvula": 25600,
    "Nodulus": 98400,
    "Flocculus": 78000,
    "Paraflocculus": 18272.36,
}


class BrainAtlasProcessor:
    def __init__(self, brt, nrrd_path="../Flocculus3.0_Lingula/data/"):
        self.nrrd_path = nrrd_path
        self.brt = brt
        self.mask = {}
        self.__brain_regions_init(self)
        self.__nrrd_read()

    def __brain_regions_init(self, brt):
        """TODO: remove this, re-use brt"""
        self.id_to_region_dictionary = self.brt.id_to_region_dictionary()
        self.id_region = self.brt.id_region
        self.id_gr, self.id_pc, self.id_mol = self.brt.get_id_gr_pc_mol()
        self.id_current_region = self.brt.region_of_interest.id
        self.region_name = self.brt.region_of_interest.name

    def __nrrd_read(self):
        # get the 3-dimensional data of cell distributions of the brain
        data_path = self.nrrd_path
        # voxel - region annotation
        self.ann, h = nrrd.read(data_path + "annotations.nrrd")
        self.dens_cell, h = nrrd.read(data_path + "cell_density.nrrd")
        self.dens_neuron, h = nrrd.read(data_path + "neu_density.nrrd")
        self.dens_inh, h = nrrd.read(data_path + "inh_density.nrrd")
        # orientation voxel
        self.orientations, h = nrrd.read(data_path + "orientations_cereb.nrrd")

    def _init_mask_regions(self):
        """Store information per region"""
        for id_ in self.brt.involved_regions:
            self.mask[str(id)] = self.ann == id

    def render_regions(self):
        return self.brt.render_roi()
