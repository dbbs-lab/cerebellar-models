"""
    Module for the Report class.
"""

from typing import List

import matplotlib.backends.backend_pdf
import numpy as np
from bsb import from_storage

from cerebellum.analysis.plots import Plot, ScaffoldPlot


class PlotTypeInfo:
    """
    Class storing information about element plotted, usually cell types or fibers.
    """

    def __init__(self, name: str, color, abbreviation: str = None):
        self.name = name
        """Name of the element."""
        self.abbreviation = abbreviation or name
        """Abbreviation of the element to print in the plot."""
        self.color = color
        """Color of the element."""


"""List of element info."""
LIST_CT_INFO = [
    PlotTypeInfo("mossy", [0.847, 0, 0.451, 1.0], "mf"),
    PlotTypeInfo("glomerulus", [0.847, 0, 0.451, 1.0], "glom"),
    PlotTypeInfo("granule", [0.7, 0.15, 0.15, 0.5], "GrC"),
    PlotTypeInfo("ascending_axon", [0.7, 0.15, 0.15, 0.5], "aa"),
    PlotTypeInfo("parallel_fiber", [0.7, 0.15, 0.15, 0.5], "pf"),
    PlotTypeInfo("unipolar_brush", [0.196, 0.808, 0.988, 1.0], "UBC"),
    PlotTypeInfo("ubc_glomerulus", [0.196, 0.808, 0.988, 1.0], "ubc_glom"),
    PlotTypeInfo("golgi", [0, 0.45, 0.7, 1.0], "GoC"),
    PlotTypeInfo("purkinje", [0.275, 0.800, 0.275, 1.0], "PC"),
    PlotTypeInfo("basket", [1, 0.647, 0, 1.0], "BC"),
    PlotTypeInfo("stellate", [1, 0.84, 0, 1.0], "SC"),
    PlotTypeInfo("io", [0.46, 0.376, 0.54, 1.0], "IO"),
    PlotTypeInfo("dcn_p", [0.3, 0.3, 0.3, 1.0], "DCN_P"),
    PlotTypeInfo("dcn_i", [0.635, 0, 0.145, 1.0], "DCN_I"),
]


class Report:
    """
    Class interfacing a list of matplotlib plots to save in an external file.
    """

    def __init__(self, pathname: str, cell_types_info: List[PlotTypeInfo] = None):
        self.reports = []
        """List of sub-reports"""
        self.plots = {}
        """Dictionary mapping the report's plots' name to their instance"""
        self.cell_types_info = cell_types_info or LIST_CT_INFO
        """List of PlotTypeInfo for each element to plot."""
        self.set_scaffold(pathname)

    @property
    def colors(self):
        """
        Dictionary from the name of the elements to plot to its color.
        """
        return {ct.name: ct.color for ct in self.cell_types_info}

    @property
    def abbreviations(self):
        """
        Dictionary from the name of the elements to plot to its abbreviation.
        """
        return {ct.name: ct.abbreviation for ct in self.cell_types_info}

    @property
    def cell_names(self):
        """
        List of the name of the elements to plot.
        """
        return list(self.scaffold.cell_types.keys())

    def set_subplot_colors(self, plot: Plot):
        """
        Set the plot color dictionary according to the report's
        """
        for k, v in self.colors.items():
            plot.set_color(k, v)

    def set_color(self, key: str, color: np.ndarray[float]):
        """
        Set a color for an element to plot.
        """
        assert len(color) == 3
        self.colors[key] = color
        for plot in self.plots.values():
            self.set_subplot_colors(plot)

    def set_scaffold(self, pathname):
        """
        Set the scaffold of the report and update each of its subplots
        """
        self.scaffold = from_storage(pathname)
        for plot in self.plots.values():
            if isinstance(plot, ScaffoldPlot):
                plot.set_scaffold(self.scaffold)

    def add_plot(self, name: str, plot: Plot):
        """
        Add a plot to the list of the report's plots.
        """
        self.set_subplot_colors(plot)
        self.plots[name] = plot
        if isinstance(plot, ScaffoldPlot):
            plot.set_scaffold(self.scaffold)

    def print_report(self, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        """
        Print the report and export it in a pdf file.
        """
        pdf = matplotlib.backends.backend_pdf.PdfPages(output_name)
        for name, plot in self.plots.items():
            if not plot.is_plotted:
                plot.plot()
            plot.figure.tight_layout(pad=pad)
            plot.figure.savefig(pdf, format="pdf", dpi=dpi, **kwargs)
        pdf.close()

    def save_plot(self, plot_name: str, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        """
        Save one of the report's plots as a separate figure.
        """
        if plot_name in self.plots:
            self.plots[plot_name].save_figure(output_name, dpi, pad, **kwargs)

    def show_plot(self, plot_name):
        """
        Show one of the report's plots.
        """
        if plot_name in self.plots:
            self.plots[plot_name].show()

    def show(self):
        """
        Show all report's plots one by one.
        """
        for plot in self.plots.values():
            plot.show()
