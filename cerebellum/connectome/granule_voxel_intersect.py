from itertools import chain

import numpy as np
from bsb import (
    AllenStructure,
    Chunk,
    ConnectivityWarning,
    MorphologySet,
    StoredMorphology,
    VoxelIntersection,
    config,
    warn,
)

from cerebellum.placement.bender_generator import GranuleBender


def load_deformed_morpho(
    placement_set,
    morpho_bender=GranuleBender(
        dict(
            rescale=["ascending_axon"],
            deform=["parallel_fiber"],
            fixed_dimension=0,
        )
    ),
):
    labels_filter = placement_set._morphology_labels
    placement_set.set_morphology_label_filter(None)
    partitions = placement_set.cell_type.scaffold.partitions
    for p in partitions:
        if isinstance(partitions[p], AllenStructure):
            morpho_bender.partition = partitions[p]
            break
    deformed_morpho = morpho_bender.process(
        placement_set.load_positions(),
        np.array([m for m in placement_set.load_morphologies().iter_morphologies()]),
    )
    deformed_morpho = [
        morpho.set_label_filter(labels_filter).as_filtered() for morpho in deformed_morpho
    ]
    u_m, u_ids = np.unique(deformed_morpho, return_inverse=True)

    return MorphologySet([StoredMorphology(i, lambda: m, m.meta) for i, m in enumerate(u_m)], u_ids)


@config.node
class GranuleToMorphologyIntersection(VoxelIntersection, classmap_entry="grc_to"):
    def get_region_of_interest(self, chunk):
        # RoI is here based solely on postsynaptic cell
        post_ps = [ct.get_placement_set() for ct in self.postsynaptic.cell_types]
        lpost, upost = self.postsynaptic._get_rect_ext(tuple(chunk.dimensions))

        bounds = list(np.arange(-u2 + c, -l2 + c + 1) for l2, u2, c in zip(lpost, upost, chunk))
        # Flatten and stack the meshgrid coordinates into a list.
        clist = np.column_stack([a.reshape(-1) for a in np.meshgrid(*bounds, indexing="ij")])
        if not hasattr(self, "_occ_chunks"):
            # Filter by chunks where cells were actually placed
            self._occ_chunks = set(chain.from_iterable(ps.get_all_chunks() for ps in post_ps))
        if not self._occ_chunks:
            warn(
                f"No {', '.join(ps.tag for ps in post_ps)} were placed, skipping {self.name}",
                ConnectivityWarning,
            )
            return []
        else:
            size = next(iter(self._occ_chunks)).dimensions
            return [t for c in clist if (t := Chunk(c, size)) in self._occ_chunks]

    def connect(self, pre, post):
        super(GranuleToMorphologyIntersection, self).connect(pre, post)
