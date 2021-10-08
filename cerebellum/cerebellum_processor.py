import nrrd
import pickle
import numpy as np
import plotly.graph_objects as go
from cerebellum.helper_mixins.bap_abstract import BrainAtlasProcessor


class CerebellumProcessor(BrainAtlasProcessor):
    def __init__(self, brt, nrrd_path):
        self.nrrd_path = nrrd_path
        self.brt = brt
        self.mask = {}
        self.vox_in_layer = {}
        self.stats = {}
        self._nrrd_read()

    def _nrrd_read(self):
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

    def screenshot_region(self, id_):
        """During masking of region, keeping track of stats.
        TODO: does not deal with parent cumulative values
        skipped for now, non-essential (not used in placement)"""

        stats = {}
        number_cells = np.round(np.sum(self.dens_cell[self.mask[id_]]))
        number_neurons = np.round(np.sum(self.dens_neuron[self.mask[id_]]))
        inhibitory_count = self.dens_inh[self.mask[id_]]
        volumes = sum(self.mask[id_]) / (4.0 ** 3) / 1000.0
        stats[id_]["number_cells"] = number_cells
        stats[id_]["number_neurons"] = number_neurons
        stats[id_]["volumes"] = volumes
        stats[id_]["cell_densities"] = number_cells / volumes
        stats[id_]["neuron_densities"] = number_neurons / volumes
        stats[id_]["number_inhibitory"] = np.round(np.sum(inhibitory_count))
        return stats

    def mask_regions(self):
        """Store information per region node in the mask for filtering."""

        for region_node in self.brt.involved_regions:
            name_region = region_node.get_name()
            id_ = region_node.id
            self.mask[name_region] = self.ann == id_
            self.vox_in_layer[name_region] = sum(self.mask[name_region])
        self.mask_basket_stellate()

    def mask_basket_stellate(self, thickness_ratio=2.0 / 3.0):
        # ratio of molecular layer space for stellate cells
        up_layer_distance = np.abs(self.bounds[0])
        down_layer_distance = np.abs(self.bounds[1])
        relative_layer_distance = np.zeros(up_layer_distance.shape)

        id_mol = self.brt.get_id_region("molecular layer")
        mask_mol = self.ann == id_mol

        up_layer_distance = up_layer_distance[mask_mol]
        down_layer_distance = down_layer_distance[mask_mol]
        relative_layer_distance[mask_mol] = up_layer_distance / (
            up_layer_distance + down_layer_distance
        )

        mask_of_stellate = relative_layer_distance * mask_mol
        mask_of_stellate = (mask_of_stellate > 0) * (mask_of_stellate < thickness_ratio)
        self.mask["Stellate layer"] = mask_of_stellate
        mask_of_basket = ~mask_of_stellate * mask_mol
        self.mask["Basket layer"] = mask_of_basket
        self.vox_in_layer["Stellate layer"] = sum(self.mask["Stellate layer"])
        self.vox_in_layer["Basket layer"] = sum(self.mask["Basket layer"])

    def fill_regions(self):
        """Cutting around the region of interest."""

        id_region = self.brt.id_region
        # Extract all region
        mask_all = np.isin(self.ann, id_region)
        # Scale to have granular layer = 1, PC layer = 2, molecular layer = 3
        region = self.ann - id_region[2] + 1
        # Outside of the region = 0
        region[~mask_all] = 0
        # To differentiate SC and BC in the ML
        region[self.mask["Stellate layer"]] = 4
        # Cut around the region
        self.region_index = np.nonzero(region)
        self.region = region[
            np.amin(self.region_index[0]) - 10 : np.amax(self.region_index[0]) + 10,  # noqa
            np.amin(self.region_index[1]) - 10 : np.amax(self.region_index[1]) + 10,  # noqa
            np.amin(self.region_index[2]) - 10 : np.amax(self.region_index[2]) + 10,  # noqa
        ]

    def show_regions(self, sliding_dir=2):
        # Define colorscale
        color_region = [
            # External to region: gray
            [0, "rgb(220, 220, 220)"],
            [0.2, "rgb(220, 220, 220)"],
            # Granular layer: red
            [0.2, "rgb(256, 0, 0)"],
            [0.4, "rgb(256, 0, 0)"],
            # PC layer: green
            [0.4, "rgb(58, 146, 94)"],
            [0.6, "rgb(58, 146, 94)"],
            # Molecular layer - BC: arancione
            [0.6, "rgb(249, 90, 8)"],
            [0.8, "rgb(249, 90, 8)"],
            # Molecular layer - SC: giallo
            [0.8, "rgb(245, 177, 0)"],
            [1.0, "rgb(245, 177, 0)"],
        ]

        fig_scatter = go.Figure(
            data=go.Heatmap(
                z=np.take(self.region, 39, axis=sliding_dir),
                colorscale=color_region
                # cmin=0, cmax=200,
                # colorbar=dict(thickness=20, ticklen=4)
            )
        )
        fig_scatter.show()

    def render_regions(self):
        return self.brt.render_roi()

    def save_processor(self, keep_arrays=True):
        """"""

        filename = "config/big_data/" + self.brt.region_name + ".pkl"
        with open(filename, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def read_processor(self, region_name="Lingula (I)"):
        filename = "config/big_data/" + region_name + ".pkl"
        with open(filename, "rb") as inp:
            return pickle.load(inp)

    def run_pipeline(self, show_regions=False, save_processor=False):
        self.mask_regions()
        self.fill_regions()
        if show_regions:
            self.show_regions()
        if save_processor:
            self.save_processor()
