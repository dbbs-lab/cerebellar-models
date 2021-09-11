import numpy as np


class BAPCerebellum(object):
    def fill_basket_stellate(self):
        # ratio of molecular layer space for stellate cells
        thickness_ratio = 2.0 / 3.0
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
