from bsb.strategy import PlacementStrategy
import math
import numpy as np
from bsb. import config
from bsb.config import types
from bsb.mixins import NotParallel
from bsb.storage import Chunk
from bsb.reporting import report, warn


@config.node
class IOPlacement(NotParallel, PlacementStrategy):
    xy_scale = config.attr(type=float, required=True)
    numbers_of_cells = config.attr(type=int, required=True)



    def place(self, chunk, indicators):
        """
        Cell placement: Create a lattice of parallel arrays/lines in the layer's surface.
        """
        for indicator in indicators.values():
            cell_type = indicator.cell_type
            radius = indicator.get_radius()
            for prt in self.partitions:
                width, depth, height = prt.data.mdc - prt.data.ldc


                #Find the coordinates of the midpoints of the volume on the x and y sides
                x_barycenter = width/2.
                y_barycenter = depth/2.

                #Find the extrema of the rescaled volume on the xy plane
                x_min = x_barycenter - self.xy_scale*width
                y_min = y_barycenter - self.xy_scale*depth
                x_max = x_barycenter + self.xy_scale*width
                y_max = y_barycenter + self.xy_scale*depth

                #Generate an array of 3d positions in [0,1)                             
                positions = np.random.rand((self.numbers_of_cells,3))
                
                #Center at the origin
                positions = positions - 0.5
                
                #Shift and rescale the array to fit in the selected volume
                positions[:,0] = x_barycenter + (x_max-x_min)*positions[:,0]
                positions[:,1] = y_barycenter + (y_max-y_min)*positions[:,1]
                positions[:,1] = y_barycenter + (y_max-y_min)*positions[:,1]
                positions[:,2] = positions[:,2] * height

                # Determine in which chunks the cells must be placed
                cs = self.scaffold.configuration.network.chunk_size
                chunks_list = np.array(
                    [chunk.data + np.floor_divide(p, cs[0]) for p in positions]
                )
                unique_chunks_list = np.unique(chunks_list, axis=0)

                # For each chunk, place the cells
                for c in unique_chunks_list:
                    idx = np.where((chunks_list == c).all(axis=1))
                    pos_current_chunk = positions[idx]
                    self.place_cells(indicator, pos_current_chunk, chunk=c)
                report(f"Placed {len(positions)} {cell_type.name} in {prt.name}", level=3)
