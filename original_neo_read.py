from neo import io
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

#Look for all the .nio files and find the more recent one
list_of_files = glob.glob('*.nio')
filename = max(list_of_files, key=os.path.getctime)

#Read sim data
sim = io.NixIO(filename, mode="ro")
blocks = sim.read_all_blocks()
block = blocks[0]
nsegments = len(block.segments[0].analogsignals)

#Plot clamps and recorders data
for i in range(nsegments):
    for seg in block.segments:
        avg = np.mean(seg.analogsignals[i], axis=1)
        #Plot and save figure to file in imgs folder
        plt.figure()
        fn = "imgs/"+str(i)+".png"
        plt.plot(avg)
        plt.savefig(fn)
        plt.close()



