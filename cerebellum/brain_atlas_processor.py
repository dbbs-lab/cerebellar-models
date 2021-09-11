import nrrd
import numpy as np

from cerebellum.helper_mixins.bap_cerebellum_mixins import BAPCerebellum

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
layers_per_cell = {
    "granule": "granular layer",
    "golgi": "granular layer",
    "purkinje": "Purkinje layer",
    "stellate": "Stellate layer",
    "basket": "Basket layer",
}


class BrainAtlasProcessor(BAPCerebellum):
    def __init__(self, brt, nrrd_path="../Flocculus3.0_Lingula/data/"):
        self.nrrd_path = nrrd_path
        self.brt = brt
        self.mask = {}
        self.vox_in_layer = {}
        self.stats = {}
        # self.__brain_regions_init(self)
        self.__nrrd_read()

        self.fill_basket_stellate()

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
        self.bounds, h = nrrd.read(data_path + "boundaries_mo.nrrd")

    def mask_regions(self):
        """Store information per region"""
        for id_ in self.brt.involved_regions:
            self.mask[id_] = self.ann == id_
            self.vox_in_layer[id_] = sum(self.mask[id_])

    def _stats_subregion_masking(self, id_):
        """During masking of region, keeping track of stats.
        TODO: does not deal with parent cumulative values
        skipped for now, non-essential (not used in placement)"""

        self.stats[id_] = {}
        number_cells = np.round(np.sum(self.dens_cell[self.mask[id_]]))
        number_neurons = np.round(np.sum(self.dens_neuron[self.mask[id_]]))
        volumes = sum(self.mask[id_]) / (4.0 ** 3) / 1000.0
        self.stats[id_]["number_cells"] = number_cells
        self.stats[id_]["number_neurons"] = number_neurons
        self.stats[id_]["volumes"] = volumes
        self.stats[id_]["cell_densities"] = number_cells / volumes
        self.stats[id_]["neuron_densities"] = number_neurons / volumes
        self.stats[id_]["number_inhibitory"] = np.round(np.sum(self.dens_inh[self.mask[id_]]))

    def render_regions(self):
        return self.brt.render_roi()
