from bsb import  FractionFilter, Targetting,config
from bsb.simulation.targetting import CellTypeFilter
from bsb.config import refs, types

@config.node
class ByIdTargettingCellTypes(CellTypeFilter, FractionFilter, Targetting, classmap_entry="by_id_cell_types"):
    """
    Targets cell types by id.
    """

    ids: dict[str, list[int]] = config.attr(
        type=types.dict(type=types.list(type=int)), required=True
    )

    @FractionFilter.filter
    def get_targets(self, adapter, simulation, simdata):
        """
        Target cell types by od
        """
        return {
            model: ps.load_ids()[self.ids[model]]
            for model, ps in super().get_targets(adapter, simulation, simdata).items()
        }
