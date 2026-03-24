import numpy as np
from bsb import AfterPlacementHook, config, refs, types


@config.node
class LabelCells(AfterPlacementHook):
    """
    Subdivide cell populations into labelled subpopulations randomly.
    The number of labels defines the number of subpopulations
    """

    cell_types: str = config.reflist(refs.cell_type_ref, required=True)
    """Reference to the cell type."""

    labels: list[str] = config.list(type=str, default=["type1", "type2"])
    """List of labels to assign to each subpopulation."""

    same_size: bool = config.attr(type=bool, default=False)
    """Flag to split in ensembles of same size"""

    def split_indexes(self, ps):
        if self.same_size:
            index_pos = np.random.permutation(len(ps))
            split_indexes = np.asarray(
                np.round(np.linspace(0, len(index_pos), len(self.labels) + 1))[1:],
                dtype=int,
            )
            return np.split(index_pos, split_indexes)[:-1]
        choice = np.random.choice(len(self.labels), len(ps))
        return [np.where(choice == i)[0] for i in range(len(self.labels))]

    def postprocess(self):
        for cell_type in self.cell_types:
            # Load the cell type positions
            ps = self.scaffold.get_placement_set(cell_type)
            split_indexes = self.split_indexes(ps)
            for indexes, label in zip(split_indexes, self.labels):
                ps.label(labels=[label], cells=indexes)


@config.node
class LabelMicrozones(LabelCells):
    """
    Subdivide cell populations into labelled subpopulations of
    same cell counts based on their position along a provided axis.
    The number of labels defines the number of subpopulations
    """

    axis: int = config.attr(type=types.int(min=0, max=2), default=0)
    """Axis along which to subdivide the population."""

    same_size = config.unset()

    def split_indexes(self, ps):
        cell_positions = ps.load_positions()
        # create a filter that split the cells according to
        # the mean of their positions along the chosen axis
        index_pos = np.argsort(cell_positions[:, self.axis])
        split_indexes = np.asarray(
            np.round(np.linspace(0, len(index_pos), len(self.labels) + 1))[1:],
            dtype=int,
        )
        return np.split(index_pos, split_indexes)[:-1]
