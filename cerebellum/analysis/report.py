from typing import List

import matplotlib.backends.backend_pdf
import numpy as np
from bsb import from_storage

from cerebellum.analysis.plots import Plot, ScaffoldPlot


class CellTypeInfo:
    def __init__(self, name: str, color, abbreviation: str = None):
        self.name = name
        self.abbreviation = abbreviation or name
        self.color = color


LIST_CT_INFO = [
    CellTypeInfo("mossy", [0.847, 0, 0.451, 1.0], "mf"),
    CellTypeInfo("glomerulus", [0.847, 0, 0.451, 1.0], "glom"),
    CellTypeInfo("granule", [0.7, 0.15, 0.15, 0.5], "GrC"),
    CellTypeInfo("ascending_axon", [0.7, 0.15, 0.15, 0.5], "aa"),
    CellTypeInfo("parallel_fiber", [0.7, 0.15, 0.15, 0.5], "pf"),
    CellTypeInfo("unipolar_brush", [0.196, 0.808, 0.988, 1.0], "UBC"),
    CellTypeInfo("ubc_glomerulus", [0.196, 0.808, 0.988, 1.0], "ubc_glom"),
    CellTypeInfo("golgi", [0, 0.45, 0.7, 1.0], "GoC"),
    CellTypeInfo("purkinje", [0.275, 0.800, 0.275, 1.0], "PC"),
    CellTypeInfo("basket", [1, 0.647, 0, 1.0], "BC"),
    CellTypeInfo("stellate", [1, 0.84, 0, 1.0], "SC"),
    CellTypeInfo("io", [0.46, 0.376, 0.54, 1.0], "IO"),
    CellTypeInfo("dcn_p", [0.3, 0.3, 0.3, 1.0], "DCN_P"),
    CellTypeInfo("dcn_i", [0.635, 0, 0.145, 1.0], "DCN_I"),
]


class Report:
    def __init__(self, pathname: str, cell_types_info: List[CellTypeInfo] = None):
        self.reports = []
        self.plots = {}
        self.cell_types_info = cell_types_info or LIST_CT_INFO
        self.set_network(pathname)

    @property
    def colors(self):
        return {ct.name: ct.color for ct in self.cell_types_info}

    @property
    def abbreviations(self):
        return {ct.name: ct.abbreviation for ct in self.cell_types_info}

    @property
    def cell_names(self):
        return list(self.network.cell_types.keys())

    def set_subplot_colors(self, plot: Plot):
        for k, v in self.colors.items():
            plot.set_color(k, v)

    def set_color(self, key: str, color: np.ndarray[float]):
        assert len(color) == 3
        self.colors[key] = color
        for plot in self.plots.values():
            self.set_subplot_colors(plot)

    def set_network(self, pathname):
        self.network = from_storage(pathname)
        for plot in self.plots.values():
            if isinstance(plot, ScaffoldPlot):
                plot.set_scaffold(self.network)

    def add_plot(self, name: str, plot: Plot):
        self.set_subplot_colors(plot)
        self.plots[name] = plot
        if isinstance(plot, ScaffoldPlot):
            plot.set_scaffold(self.network)

    def print_report(self, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        pdf = matplotlib.backends.backend_pdf.PdfPages(output_name)
        for name, plot in self.plots.items():
            if not plot.is_plotted:
                plot.plot()
            plot.figure.tight_layout(pad=pad)
            plot.figure.savefig(pdf, format="pdf", dpi=dpi, **kwargs)
        pdf.close()

    def save_plot(self, plot_name: str, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        if plot_name in self.plots:
            self.plots[plot_name].save_figure(output_name, dpi, pad, **kwargs)

    def show_plot(self, plot_name):
        if plot_name in self.plots:
            self.plots[plot_name].show()

    def show(self):
        for plot in self.plots.values():
            plot.show()
