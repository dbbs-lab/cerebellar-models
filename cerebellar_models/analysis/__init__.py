"""
Analysis and plotting tools for cerebellar cortex reconstructions.
"""

import matplotlib

from cerebellar_models.analysis.plots import Legend, Plot, ScaffoldPlot
from cerebellar_models.analysis.report import (
    LIST_CT_INFO,
    BSBReport,
    PlotTypeInfo,
    Report,
)
from cerebellar_models.analysis.spiking_results import (
    BasicSimulationReport,
    FiringRatesPlot,
    FrequencyPlot,
    ISIPlot,
    RasterPSTHPlot,
    SimResultsTable,
    SpikePlot,
    SpikeSimulationReport,
)
from cerebellar_models.analysis.structure_analysis import (
    CellPlacement3D,
    ConnectivityTable,
    PlacementTable,
    RunStructureReport,
    StructureReport,
    TablePlot,
)

matplotlib.use("Agg")

__all__ = [
    "Legend",
    "Plot",
    "ScaffoldPlot",
    "Report",
    "BSBReport",
    "PlotTypeInfo",
    "LIST_CT_INFO",
    "BasicSimulationReport",
    "FiringRatesPlot",
    "FrequencyPlot",
    "SpikePlot",
    "SpikeSimulationReport",
    "SimResultsTable",
    "ISIPlot",
    "RasterPSTHPlot",
    "TablePlot",
    "PlacementTable",
    "ConnectivityTable",
    "CellPlacement3D",
    "StructureReport",
    "RunStructureReport",
]
