import pandas as pd
import matplotlib.pyplot as plt
import os

def show_traces(file_path, synapse, conn, save_path):
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    if synapse=='AMPA':
        df["g"] = df["Current"].values * 1000 / (-40.)
    elif synapse=='GABA':
        df["g"] = -df["Current"].values * 1000 / (-40.)
    else:
        df["g"] = df["Current"].values * 1000 / ((-40.+3.7)  * 0.177)  # 0.177 is the mg_block value associated to -40 mV,
                                                                            # -3.7 is the reversal potential for NMDA

    plt.figure()
    plt.title(f'{synapse} current - {conn} (NEURON sim)')
    plt.plot(df["Time"], df['Current'].values*1000, color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel(f'{synapse} current [pA]')
    plt.xlim(0,1000)
    plt.savefig(save_path + f'/{conn}_{synapse}_current.png', dpi=150)
    plt.show()

    plt.figure()
    plt.title(f'{synapse} conductance - {conn} (NEURON sim)')
    plt.plot(df["Time"], df['g'], color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel(f'{synapse} conductance [nS]')
    plt.xlim(0,1000)
    plt.savefig(save_path + f'/{conn}_{synapse}_conductance.png', dpi=150)
    plt.show()

if __name__ == '__main__':
    os.makedirs('figs/Neuron_Danilo/', exist_ok=True)

    file_path = "../NEURON_traces/Danilo_synapses/SCpfAMPA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_SC'
    synapse = 'AMPA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/SCGABA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'SC'
    synapse = 'GABA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GoCpfAMPA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_GoC'
    synapse = 'AMPA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GoCmfAMPA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'mf_GoC'
    synapse = 'AMPA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GoCGABA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'GoC'
    synapse = 'GABA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GoCaaAMPA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'aa_GoC'
    synapse = 'AMPA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/BCpfNMDA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_BC'
    synapse = 'NMDA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/BCpfAMPA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_BC'
    synapse = 'AMPA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/BCGABA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'BC'
    synapse = 'GABA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GrCGABA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'Grc'
    synapse = 'GABA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/SCpfNMDA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_SC'
    synapse = 'NMDA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/SCpfNMDA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'pf_SC'
    synapse = 'NMDA'
    show_traces(file_path, synapse, conn, save_path)

    file_path = "../NEURON_traces/Danilo_synapses/GoCmfNMDA.txt"
    save_path = 'figs/Neuron_Danilo/'
    conn = 'mf_GoC'
    synapse = 'NMDA'
    show_traces(file_path, synapse, conn, save_path)









