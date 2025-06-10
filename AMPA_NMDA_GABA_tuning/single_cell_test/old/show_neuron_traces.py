import pandas as pd
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    os.makedirs('../figs/Neuron_Masoli/', exist_ok=True)

    file_path1 = "../../NEURON_traces/Masoli_synapses/AMPA_GrC.csv"
    df = pd.read_csv(file_path1, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_AMPA"] = df["Current"].values * 1000 / (-40.)

    plt.figure()
    plt.title('AMPA current - mf_GrC (NEURON sim)')
    plt.plot(df["Time"], df['Current'].values*1000, color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('AMPA current mf_Grc [pA]')
    plt.xlim(0,1000)
    plt.savefig('../figs/Neuron_Masoli/mf_GrC_AMPA_current.png', dpi=150)
    plt.show()

    plt.figure()
    plt.title('AMPA conductance - mf_GrC (NEURON sim)')
    plt.plot(df["Time"], df['g_AMPA'], color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('AMPA conductance mf_GrC [nS]')
    plt.xlim(0,1000)
    plt.savefig('../figs/Neuron_Masoli/mf_GrC_AMPA_conductance.png', dpi=150)
    plt.show()

    file_path2 = "../../NEURON_traces/Masoli_synapses/NMDA_GrC.csv"
    df = pd.read_csv(file_path2, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_NMDA"] = df["Current"].values * 1000 / ((-40.+3.7)  * 0.177) # 0.177 is the mg_block value associated to -40 mV,
                                                                       # -3.7 is the reversal potential for NMDA

    plt.figure()
    plt.title('NMDA current - mf_GrC (NEURON sim)')
    plt.plot(df["Time"], df['Current'].values*1000, color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('NMDA current [pA]')
    plt.xlim(0,1000)
    plt.savefig('../figs/Neuron_Masoli/mf_GrC_NMDA_current.png', dpi=150)
    plt.show()


    plt.figure()
    plt.title('NMDA conductance - mf_GrC (NEURON sim)')
    plt.plot(df["Time"], df['g_NMDA'], color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('NMDA conductance [nS]')
    plt.xlim(0,1000)
    plt.savefig('../figs/Neuron_Masoli/mf_Grc_NMDA_conductace.png', dpi=150)
    plt.show()


