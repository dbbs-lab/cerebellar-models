#!/usr/bin/env python
# coding: utf-8

# In[20]:


import json

contenido={
  "name": "DBBS Mouse Lingula configuration v4.0",
  "storage": {"engine": "hdf5", "root": "lingula.hdf5" },
  
  "network": {"x": 400.0,"y": 400.0, "z": 400.0,"chunk_size": [50,50,50] },

  "regions": {  
    "cerebellar_cortex":{"type": "stack","children": ["granular", "purkinje", "b_molecular", "t_molecular" ]},},
 
  "partitions": {
    "granular": {"cls": "layer","region": "cerebellar_cortex","thickness": 150.0,"stack_index": 1},
    "purkinje": {"cls": "layer","region": "cerebellar_cortex","thickness": 30.0, "stack_index": 2},
    "b_molecular": {"cls": "layer","region": "cerebellar_cortex","thickness": 50.0,"stack_index": 3},
    "t_molecular": {"cls": "layer","region": "cerebellar_cortex","thickness": 100.0,"stack_index": 4},},
 
  "cell_types": {
    "granule_cell": {"spatial": {"radius": 2.5,"density": 3.9e-3,"geometrical": {"pf_height": 180, "pf_height_sd": 20}, "morphological": [{"names": ["GranuleCell"]}]},"plotting": {"display_name": "Granule cell","color": "#e81005", "opacity": 0.3}},
    "mossy_fibers": {"entity": True,"spatial": {"relative_to": "glomerulus","count_ratio": 0.05}},
    "glomerulus": {"relay": True,"spatial": {"radius": 1.5,"density": 3e-4,"plotting": {"display_name": "Glomerulus","color": "#6F6F70"}},},
    "purkinje_cell": {"spatial": { "radius": 7.5,"planar_density": 0.0017, "morphologies": [ {"names": ["PurkinjeCell"]}]},"plotting": {"display_name": "Purkinje cell", "color": "#068f0d"}},
    "golgi_cell": {"spatial": {"radius": 8.0,"density": 9e-6,"morphologies":[{ "names": ["GolgiCell"]}]},"plotting": {"display_name": "Golgi cell","color": "#1009e3"}},
    "stellate_cell": {"spatial": {"radius": 4.0,"density": 0.5e-4,"morphologies":[{"names": ["StellateCell"]}]},"plotting": {"display_name": "Stellate cell","color": "#f5bb1d"}},
    "basket_cell": {"spatial": {"radius": 6.0,"density": 0.5e-4,"morphologies":[{"names": ["BasketCell"]}]},"plotting": {"display_name": "Basket cell","color": "#f5830a"}},}, 

  "placement": {
   "granular_innervation": {"cls": "bsb.placement.Entities","partitions": ["granular"],"cell_types": ["mossy_fibers"]},
   "granular_placement": {"cls": "bsb.placement.ParticlePlacement", "partitions": ["granular"],"cell_types": ["granule_cell", "golgi_cell", "glomerulus"]},
   "purkinje_placement": {"cls": "bsb.placement.ParallelArrayPlacement","partitions": ["purkinje"],"cell_types": ["purkinje_cell"],"spacing_x": 130.0, "angle": 70.0}, 
   "molecular_layer_placement": {"cls": "bsb.placement.ParticlePlacement","partitions": ["b_molecular", "t_molecular"],"cell_types": ["stellate_cell", "basket_cell"],"restrict": {"basket_cell": ["b_molecular"],"stellate_cell": ["t_molecular"]}}},
  
  "after_placement": {
   "microzones": { "cls": "bsb.postprocessing.LabelMicrozones",   "targets": ["purkinje_cell"]},
   "aa_lengths": {"cls": "bsb.postprocessing.AscendingAxonLengths"},},   
    
  "connectivity": {
    "mossy_to_glomerulus": {"cls": "mossy_glomerulus.ConnectomeMossyGlomerulus", "presynaptic": {"cell_types": ["mossy_fibers"],"compartments": ["soma"]},"postsynaptic": {"cell_types": ["glomerulus"],"compartments": ["soma"]},},
    "glomerulus_to_granule": {"cls": "glomerulus_granule.ConnectomeGlomerulusGranule","presynaptic": {"cell_types": ["glomerulus"],"compartments": ["soma"]},"postsynaptic": {"cell_types": ["granule_cell"],"compartments": ["dendrites"]},"divergence": 50, "convergence": 4},
    "glomerulus_to_golgi": {"cls": "glomerulus_golgi.ConnectomeGlomerulusGolgi","presynaptic": { "cell_types": ["glomerulus"],"compartments": ["soma"]},"postsynaptic": {"cell_types": ["golgi_cell"],"compartments": ["dendrites"]},"divergence": 2, "convergence": 6},
    "golgi_to_glomerulus": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["golgi_cell"],"compartments": ["axon" ]},"postsynaptic": {"cell_types": ["glomerulus"],"compartments": ["soma"],}, "divergence": 40,  "convergence": 4},      
    "granule_to_golgi": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["granule_cell"], "compartments": ["parallel_fiber", "ascending_axon" ]},"postsynaptic":{"cell_types": ["golgi_cell"], "compartments": ["dendrites"]}, "tag_aa":"ascending_axon_to_golgi","tag_pf":"parallel_fiber_to_golgi", "aa_convergence": 400, "pf_convergence": 1200},
    "golgi_to_granule": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": { "cell_types": ["golgi_cell"], "compartments": ["axon"]},"postsynaptic": {"cell_types": ["granule_cell"], "compartments": ["dendrites"],"after": ["golgi_to_glomerulus", "glomerulus_to_granule"] },},
    "ascending_axon_to_purkinje": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": { "cell_types": ["granule_cell"], "compartments": ["ascending_axon"] }, "postsynaptic": { "cell_types": ["purkinje_cell"], "compartments": ["aa_targets"]  },},
    "parallel_fiber_to_purkinje": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": { "cell_types": ["granule_cell"], "compartments": ["parallel_fiber"]},"postsynaptic": { "cell_types": ["purkinje_cell"], "compartments": ["pf_targets"] },},
    "parallel_fiber_to_basket": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["granule_cell"],"compartments": ["parallel_fiber"]}, "postsynaptic": { "cell_types": ["basket_cell"], "compartments": ["dendrites"] ,"after": ["granule_to_golgi"] },},
    "parallel_fiber_to_stellate": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": { "cell_types": ["granule_cell"], "compartments": ["parallel_fiber"]},"postsynaptic": {"cell_types": ["stellate_cell"], "compartments": ["dendrites"], "after": ["granule_to_golgi"]  },},
    "basket_to_purkinje": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["basket_cell"],"compartments": ["axon"],},"postsynaptic": { "cell_types": ["purkinje_cell"], "compartments": ["soma"]  },  },
    "stellate_to_purkinje": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["stellate_cell"],"compartments": ["axon"]},"postsynaptic": {"cell_types": ["purkinje_cell"],"compartments": ["soma"]  }, },
    "stellate_to_stellate": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": { "cell_types": ["stellate_cell"], "compartments": ["axon"]},"postsynaptic": {"cell_types": ["stellate_cell"],"compartments": ["dendrites"]  },},
    "basket_to_basket": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["basket_cell"],"compartments": ["axon"]},"postsynaptic": { "cell_types": ["basket_cell"],"compartments": ["dendrites"]},  },
    "golgi_to_golgi": {"cls": "bsb.connectivity.VoxelIntersection","presynaptic": {"cell_types": ["golgi_cell"],"compartments": ["axon"]},"postsynaptic": { "cell_types": ["golgi_cell"], "compartments": ["dendrites"] }}, } ,

}


out_file=open("mouse_lingulaV4.json","w")
json.dump(contenido,out_file,indent=6)
out_file.close()


# In[ ]:




