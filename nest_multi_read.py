from neo import io
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.signal import find_peaks

#Look for all the .nio files and find the more recent one
list_of_files = glob.glob('*.nio')
filename = max(list_of_files, key=os.path.getctime)

#Read sim data
sim = io.NixIO(filename, mode="ro")
blocks = sim.read_all_blocks()
print(len(blocks))
block = blocks[0]


nsegments = len(block.segments[0].spiketrains)

simulation_time = 200
resolution = 0.25
time = np.arange(0,simulation_time+resolution,resolution)

#--------------------------------------


mossy_pop_size = 30
gloms_pop_size = 600
gran_pop_size  = 15003

gloms_spikes = []
for i in range(gloms_pop_size):
    gloms_spikes.append([])

gran_spikes = []
for i in range(gran_pop_size):
        gran_spikes.append([])

mossy_spikes = []
for i in range(mossy_pop_size):
    mossy_spikes.append([])

filenames = ["1.nio","2.nio","3.nio","4.nio"]

for filename in filenames:

    sim = io.NixIO(filename, mode="ro")
    blocks = sim.read_all_blocks()
    block = blocks[0]
    nsegments = len(block.segments[0].spiketrains)
    spiketrains = block.segments[0].spiketrains
    dev_gen = None
    rec_mossy = None
    rec_gloms = None
    rec_gran = None


    for st in spiketrains:
        if st.annotations["device"] == "mossy_rec":
            rec_mossy = st
        if st.annotations["device"] == "pg":
            dev_gen = st
        if st.annotations["device"] == "gloms_rec":
            rec_gloms = st
        if st.annotations["device"] == "granule_rec":
            rec_gran = st

    #Extract poisson spike
    gen_spikes = []
    mossy_pop_size = dev_gen.annotations
    print("mossy pop size:",mossy_pop_size)
    print("Unique senders:", len(np.unique(rec_mossy.annotations["senders"])))

    #Extract mossy spikes
    mossy_pop_size = rec_mossy.annotations["pop_size"]
    print("mossy pop size:",mossy_pop_size)
    print("Unique senders:", len(np.unique(rec_mossy.annotations["senders"])))

    ptr = 0
    for i in range(len(rec_mossy)):
        mossy_spikes[rec_mossy.annotations["senders"][i]-1-ptr].append(rec_mossy[i])
        #print(rec_mossy)
        #print(rec_mossy.annotations)

    #Extract mossy spikes
    gloms_pop_size = rec_gloms.annotations["pop_size"]
    print("mossy pop size:",gloms_pop_size)
    print("Unique senders:", len(np.unique(rec_gloms.annotations["senders"])))
    ptr = mossy_pop_size
    for i in range(len(rec_gloms)):
        gloms_spikes[rec_gloms.annotations["senders"][i]-1-ptr].append(rec_gloms[i])

    #Extract granule spikes
    gran_pop_size = rec_gran.annotations["pop_size"]
    print("mossy pop size:",gran_pop_size)
    print("Unique senders:", len(np.unique(rec_gran.annotations["senders"])))
    
    ptr = mossy_pop_size+gloms_pop_size
    for i in range(len(rec_gran)):
        gran_spikes[rec_gran.annotations["senders"][i]-1-ptr].append(rec_gran[i])

pre_neurons = 1
pop_current = 37
global_spikes = []

for i in mossy_spikes:
    global_spikes.append(i)
for i in gloms_spikes:
    global_spikes.append(i)


#spikes = mossy_spikes[pre_neurons-1:pre_neurons+pop_current-1]
#spikes = global_spikes[:1200]
spikes = gloms_spikes[:30]

fig, ax = plt.subplots()
# Loop to plot raster for each trial
for i,sp in enumerate(spikes):
    #print(sp)
    cl = 'blue'
    #if i > mossy_pop_size:
    #    cl = 'red'
    ax.vlines(sp, i - 0.5, i + 0.5, colors=cl)

ax.set_xlim([0, 1000])
#ax.set_xlabel('Time (ms)')

# specify tick marks and label label y axis
ax.set_yticks(range(len(spikes)))
ax.set_ylabel('Neuron ID')

ax.set_title('Spike Times') 

# add shading for stimulus duration)
#ax.axvspan(light_onset_time, light_offset_time, alpha=0.5, color='greenyellow')

plt.show()

#for seg in block.segments:
#   print(seg.spiketrains)
"""

#Plot histogram

tot_peaks = []
for i in range(84,99):
    
    for seg in block.segments:
        avg = np.mean(seg.analogsignals[i], axis=1)
        peaks, _ = find_peaks(avg, height=-50)
        peaks = peaks * 0.25
        #print(peaks)
        tot_peaks = np.concatenate([tot_peaks,peaks])
        #spikes.append(peaks)

count, bins = np.histogram(tot_peaks,bins=50,range=(0,250))
plt.bar(bins[:-1], count, width=bins[1] - bins[0])
plt.show()


#Plot PSTH

spikes = []
for i in range(84,99):
    tot_peaks = []
    for seg in block.segments:
        avg = np.mean(seg.analogsignals[i], axis=1)
        peaks, _ = find_peaks(avg, height=-50)
        peaks = peaks * 0.25
        #print(peaks)
        #tot_peaks = np.concatenate([tot_peaks,peaks])
        spikes.append(peaks)

fig, ax = plt.subplots()

# Loop to plot raster for each trial
for i,sp in enumerate(spikes):
    ax.vlines(sp, i - 0.5, i + 0.5)

#ax.set_xlim([0, len()])
#ax.set_xlabel('Time (ms)')

# specify tick marks and label label y axis
ax.set_yticks(range(len(spikes)))
ax.set_ylabel('Trial Number')

ax.set_title('Neuronal Spike Times') 

# add shading for stimulus duration)
#ax.axvspan(light_onset_time, light_offset_time, alpha=0.5, color='greenyellow')

plt.show()
"""
