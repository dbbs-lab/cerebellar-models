import nrrd
import sparse
import numpy as np
import plotly.graph_objects as go
from cerebellum.helper_mixins.bap_abstract import BrainAtlasProcessor


class Selector(BrainAtlasProcessor):
    def __init__(self, brt, nrrd_path):
        self.nrrd_path = nrrd_path
        self.brt = brt
        self.mask = {}
        self.vox_in_layer = {}
        self.stats = {}

    @staticmethod
    def read_nrrd(filename, nrrd_path, sparce_matrix=False):
        print(filename)
        if sparce_matrix:
            return sparse.COO.from_numpy(nrrd.read(nrrd_path + filename)[0])
        else:
            return nrrd.read(nrrd_path + filename)[0]

    @property
    def ann(self):
        if not hasattr(self, "_ann"):
            self._ann = self.read_nrrd("annotations.nrrd", self.nrrd_path, True)
        return self._ann.todense()

    @property
    def bounds(self):
        if not hasattr(self, "_bounds"):
            self._bounds = self.read_nrrd("boundaries_mo.nrrd", self.nrrd_path, True)
        return self._bounds.todense()

    def create_mask(self, brt=None):
        """Store information per region node in the mask for filtering."""

        for region_node in self.brt.involved_regions:
            name_region = region_node.get_name()
            id_ = region_node.id
            self.mask[name_region] = self.ann == id_
            self.vox_in_layer[name_region] = sum(self.mask[name_region])
        self.create_mask_basket_stellate()

    def create_mask_basket_stellate(self, thickness_ratio=2.0 / 3.0):
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

    def apply_mask(self):
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
