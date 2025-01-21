For bulding the "custom_stdp_folder", make a "build" folder.

TERMINAL COMMAND FROM THE BUILD FOLDER DIRECTORY:
cmake -Dwith-nest=<path/to/NEST-CONFIG/install/folder> <path/to/custom_stdp> && make
make install

WARNING: remember to empty the build folder everytime you build the neuron and synapsse codes
