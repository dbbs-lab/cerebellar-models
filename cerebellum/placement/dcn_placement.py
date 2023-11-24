from bsb.strategy import PlacementStrategy
import math
import numpy as np
from bsb. import config
from bsb.config import types
from bsb.mixins import NotParallel
from bsb.storage import Chunk
from bsb.reporting import report, warn
from scipy.spatial.transform import Rotation as R

@config.node
class DCNPlacement(NotParallel, PlacementStrategy):
    xy_scale = config.attr(type=float, required=True)
    number_of_dcnp_cells = config.attr(type=int, required=True)
    satellite_cell_type = config.attr(type=str, required=True)
    dcnp_radius = config.attr(type=float, required=True)
    dcngaba_radius = config.attr(type=float, required=True)



    def place(self, chunk, indicators):
        """
        Cell placement: Create a lattice of parallel arrays/lines in the layer's surface.
        """

        #Both type of DCN cells need to be positioned on a plane which is rotated with respect to the xy plane. 
        #We generate the common rotation angle before processing the two cell_type separately
        rotation_vector = np.random.rand(3,)*np.pi/2.
        rotation = R.from_rotvec(rotation_vector)

        dcnp_pos = np.random.rand((self.number_of_dcnp_cells,3))

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

                #TODO: Divide in two different strategies

                if (cell_type != self.satellite_cell_type):
                    
                    #First we place DCNp cells freely on the xy plane..
                    positions = np.random.rand((self.number_of_dcnp_cells,3))
                    positions[:,0] = x_barycenter + (x_max-x_min)*positions[:,0]
                    positions[:,1] = y_barycenter + (y_max-y_min)*positions[:,1]
                    positions[:,2] = 0.
                    dcnp_pos = positions

                    #..then we rotate them...
                    positions = rotation.apply(positions)

                    #...and we shift them
                    rescaled_vol_lcd = np.array([x_min,y_min,prt.data.ldc[2]])
                    rotated_rescaled_vol_lcd = rotation.apply(rescaled_vol_lcd)
                    shift = rotated_rescaled_vol_lcd - rescaled_vol_lcd
                    positions = positions + shift
                
                else:
                    #We place DCNgaba cells freely on the xy plane
                    dmin = self.dcnp_radius + self.dcngaba_radius
                    mean_dist = np.mean(dcnp_pos)
                    dmax = mean_dist/4. - (self.dcnp_radius + self.dcngaba_radius)

                    radii = dmin + np.random.rand((self.number_of_dcnp_cells,))*(dmax-dmin)

                    phi = np.random.rand((self.number_of_dcnp_cells,))*np.pi
                    theta = np.random.rand((self.number_of_dcnp_cells,))*np.pi*2
                    angles = np.random.rand((self.number_of_dcnp_cells,))*np.pi*2

                    #Generate random points of the surface of a sphere
                    x = np.outer(np.sin(theta), np.cos(phi))
                    y = np.outer(np.sin(theta), np.sin(phi))
                    z = np.outer(np.cos(theta), np.ones_like(phi))

                    #Generate random points inside a sphere
                    x = np.outer(x,radii)
                    y = np.outer(y,radii)
                    z = np.outer(z,radii)

                    #positions = np.full((self.number_of_dcnp_cells,3),dtype=float)
                    positions = np.copy(dcnp_pos)
                    positions[:,0] = positions[:,0] + x
                    positions[:,1] = positions[:,1] + y
                    positions[:,2] = positions[:,1] + z


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
