import nrrd
import sparse
import numpy as np


def load_nrrd(file_path, sparse_matrix=False):
    if sparse_matrix:
        return sparse.COO.from_numpy(nrrd.read_data(file_path))
    else:
        return nrrd.read_data(file_path)


class RegionSelector:
    """A clean processor of the nrrd file processing.
    Used in the following manner.

    ```
    with Selector(ann_nrrd, bounds_nrrd) as s:
        mask_1 = s.select(ids_1)
        mask_2 = s.select(ids_2)
    ```
    """

    def __init__(self, annotations_nrrd_file, boundaries_nrrd_file):
        """Reads the annotation nrrd data, prepares for filtering."""

        self._annontations = load_nrrd(annotations_nrrd_file, True)
        self._boundaries = load_nrrd(boundaries_nrrd_file, True)

    def select(self, ids=[], surrounded=False):
        """Select a specific region (with its subregions) in the brain."""

        region = np.isin(self._annontations, ids)

        if surrounded:
            region_index = np.nonzero(region)
            region = region[
                np.amin(region_index[0]) - 10 : np.amax(region_index[0]) + 10,  # noqa
                np.amin(region_index[1]) - 10 : np.amax(region_index[1]) + 10,  # noqa
                np.amin(region_index[2]) - 10 : np.amax(region_index[2]) + 10,  # noqa
            ]
        return region
