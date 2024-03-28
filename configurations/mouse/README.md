
# Introduction
This folder contains the configurations for the reconstruction and simulation of the mouse 
cerebellar cortex with BSB v4.\
These reconstructions are based on the iterative work of many researchers distributed in many 
papers.
The role of this file is to make explicit the origin of each value and strategy extracted from 
the literature and integrated into these configurations.\
All the configurations present in these folder are based on the `mouse_cerebellar_cortex.yaml` 
file. It corresponds to the configuration file written for the reconstruction of the cerebellar 
circuit presented in the `De Schepper et al. 2022` paper. This circuit configuration will be later 
referred to as the `canonical circuit`.

We will follow the structure of the BSB configuration files to present each of their sections and 
the data they leverage.

# Circuit configuration
## Introduction
### Biological context
The cerebellar cortex is a subregion of the cerebellum structured as a folded sheet.
This cortical sheet is itself decomposed into 2 major layers: the granular layer and the molecular 
layer, separated by a thin sheet of Purkinje cells, representing the Purkinje layer.
Each layer contains its own cell types composition, placement rules and connectivity strategies.\
In the `canonical circuit`, only a portion of the mouse cerebellar cortex is built in a cubic 
volume with fixed layer thicknesses, to simplify the orientation and scaling of neuronal 
geometrical representation (morphologies) but also to limit the number of neurons and therefore the 
final size of the network.

### Coordinate framework 
By convention, the circuit is oriented so that its layers are stacked vertically with the granular 
layer at the bottom and the molecular layer at the top. We derived a coordinate framework `(x,y,z)` 
from this layout, based on the right-hand orientation convention. Its origin is set at the bottom of
the circuit, the z axis pointing to the top of the molecular layer. The `(y-z)` plane corresponds 
to the para-sagittal sections that are co-planar with the Purkinje dendritic trees and normal to the
granule cells parallel fibers. Finally, unit of distance in the configurations are expressed in 
micrometers `um`. Note that the morphologies provided are oriented by default to match this 
convention.

## Network dimensions
The `canonical circuit` is built in a cubic volume of `300 * 200 * 295 um` in the `(x,y,z)` 
convention (see `network`, `regions` and `partitions` in the configuration). The thickness of each 
of its layer has been determined according to literature findings and to match the size and shape of 
the available morphologies:
- The Purkinje layer corresponds to a one cell thick sheet of Purkinje cells. The Purkinje cell 
soma diameters determine therefore the thickness of this layer. According to 
`Hendelman & Aggerwal, 1980`, Purkinje cell's soma diameters have been estimated to less than 20 um 
in mice. We chose here `15um`.
- The molecular layer total thickness has been calculated to fit the size of the dendritic 
arborization of the Purkinje cell's morphology as `150um`. The molecular layer is itself divided in 
the `canonical circuit` into two sub-layers based on their neuronal composition. 
Here, the bottom part of the molecular layer (hence the part stacked right on top of the Purkinje 
layer) contains the Basket cells (`b_molecular_layer` in the configuration); while the top 
part hold the Stellate cells (`s_molecular_layer`). In fact, according to literature data 
(`J. Kim & Augustine, 2021; Sultan & Bower, 1998`), SCs are more likely located in the outer 
two-third of the Molecular layer. While this distribution of cells is closer to gradient in real 
mice, we assumed a clear separation between the populations. The `basket layer` is therefore `50um` 
thick while the `stellate layer` is `100um` thick.
- The granular layer's thickness has been similarly fitted to match the size of the Golgi cell 
basal dendritic tree, here `130um`. Note that the size of the granule cell ascending axons have been
set to this constraint.

## Cellular composition
The cellular composition of the circuit is determined in the `cell_types` section of the 
configuration file. Each cell type is linked here to a partition in the circuit, and a morphology 
is assigned. Additionally, we introduce here two components used to innervate the circuit: 
the mossy fibers originating from various other brain regions, and the glomerulus that form at their 
terminals. These components are used as building entities to relay stimuli from other regions 
into the circuit and have therefore no morphology attached.\
We will describe here the spatial parameters used in `canonical circuit`:

| Layer             | Cell name            | Type   | Radius (um)   | Density (um^{-3})                   | References              |
|-------------------|----------------------|--------|---------------|-------------------------------------|-------------------------|
| Granular layer    | Glomerulus (glom)    | Exc.   | 1.5           | 0.0003                              | `Solinas et al., 2010`  |
|                   | Mossy fibers (mf)    | Exc.   | /             | count relative to glom. ratio=0.05  | `Billings et al., 2014` |
|                   | Granule Cell (GrC)   | Exc.   | 2.5           | 0.0039                              | `[CITATION]`            |
|                   | Golgi cells (GoC)    | Inh.   | 4.            | 0.000009                            | `[CITATION]`            |
| Purkinje layer    | Purkinje cell (PC)   | Inh.   | 7.5           | planar density: 0.0017              | `[CITATION]`            |
| Molecular layer   | Basket cell (BC)     | Inh.   | 6.            | 0.00005                             | `[CITATION]`            |
|                   | Stellate cell (SC)   | Inh.   | 4.            | 0.00005                             | `[CITATION]`            |


Note that every literature data in this table comes from the rat (except for the Purkinje layer 
planar density).\
The density of `glom` have been calculated in `Solinas et al., 2010` based on the glomerulus to 
granule convergence and divergence ratios (derived from values in `Korbo et al., 1993` and 
`Jakab and Hámari, 1988`).\

## Morphologies
[TODO: Create readme in morphologies and link it here.]

## Placement
Except for Purkinje cells (PC), every entity is supposed to be uniformly distributed in their own 
layer.The bsb `RandomPlacement` strategy is chosen here to place them. In short, this strategy chose 
a random position for each entity within their sub-partition. Note that this does not take into 
account any potential overlapping of cells' soma unlike the `ParticlePlacement`. However, 
comparative analysis conducted in [CITATION] have shown that the latter strategy have a limited 
impact on connectivity and simulation results, while the computational cost of checking soma 
overlapping is not negligible.

PC are placed in arrays, `130 μm` apart from each other along the para-sagittal plane `(xz)` to 
guarantee that their dendritic arborizations do not overlap. Furthermore, each row of PC somas is 
shifted with respect to its predecessor to form a `70` degree angle on the `(xy)` plane.

## Connectivity

| Source | Target | Strategy                           | References                  |
|--------|--------|------------------------------------|-----------------------------|
| mf     | glom   | `ConnectomeMossyGlomerulus`        | `Sultan 2001`               |  
| glom   | GoC    | `ConnectomeGlomerulusGolgi`        | `Kanichay and Silver, 2008` |
| glom   | GrC    | `ConnectomeGlomerulusGranule`      | `Houston et al., 2017`      |
| GoC    | GrC    | `ConnectomeGolgiGlomerulusGranule` | `[CITATION]`                |

### `ConnectomeMossyGlomerulus` 
According to literature data (Billings et al., 2014; Ito, 1984; Sultan, 2001), 
The algorithm selects a mf within the `60` μm along the x-axis and `20` μm along the 
z-axis box surrounding each glom. This selection is random and performed with a truncated 
exponential distribution. Since the placement of mf and glom is uniformed within their partition, 
the convergence and divergence ratios of this connection is guaranteed.

### `ConnectomeGlomerulusGolgi` 
The algorithm selects all glom within the sphere (radius `50` μm) surrounding each GoC soma. 
For each unique glom to connect, the tip of a basal dendrite branch from the golgi morphology is 
selected. This selection is random and performed with a truncated exponential distribution based on 
the distance between the tip of each branch and the glom to connect.

### `ConnectomeGlomerulusGranule` 
The algorithm selects `4` unique glom within the sphere (radius `40` μm) surrounding each GrC soma.
Moreover, each presynaptic glom should belong to a unique mf cluster, i.e. each should be connected 
through the `ConnectomeMossyGlomerulus` strategy to a different mf. The mf cluster, the presynaptic 
glom and the postsynaptic GrC dendrite are all randomly chosen. If not enough glomerulus could be 
found in the `40` μm radius sphere surrounding the GrC soma, the closest glom from the remaining 
cluster are selected to connect. 

### `ConnectomeGolgiGlomerulusGranule` 
The algorithm selects the closest glom (maximum `40`) are within the sphere (radius `150` μm) 
surrounding each GoC soma. For each unique glom selected, the tip of an axon branch from the golgi 
morphology is randomly selected. All GrC connected to the selected glom through the 
`ConnectomeGlomerulusGranule` strategy are also connected to the selected presynaptic GoC axon tip.

## References

[TODO: Add all citations here]