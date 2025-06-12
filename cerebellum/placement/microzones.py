import numpy as np
from bsb import AfterPlacementHook, config, refs, types


@config.node
class LabelMicrozones(AfterPlacementHook):
    """
    Subdivide cell populations into labelled subpopulations of
    same cell counts based on their position along a provided axis.
    The number of labels defines the number of subpopulations
    """

    cell_types: str = config.reflist(refs.cell_type_ref, required=True)
    """Reference to the cell type."""

    axis: int = config.attr(type=types.int(min=0, max=2), default=0)
    """Axis along which to subdivide the population."""

    labels: list[str] = config.list(type=str, default=["type1", "type2"])
    """List of labels to assign to each subpopulation."""

    def postprocess(self):
        for cell_type in self.cell_types:
            # Load the cell type positions
            ps = self.scaffold.get_placement_set(cell_type)
            cell_positions = ps.load_positions()

            # create a filter that split the cells according to
            # the mean of their positions along the chosen axis
            index_pos = np.argsort(cell_positions[:, self.axis])
            split_indexes = np.asarray(
                np.round(np.linspace(0, len(index_pos), len(self.labels) + 1))[1:], dtype=int
            )
            last_i = 0
            for i, label in zip(split_indexes, self.labels):
                subpopulation_1 = np.asarray(index_pos[last_i:i], dtype=int)
                # set the cell label according to the filter
                ps.label(labels=[label], cells=subpopulation_1)
                last_i = i
