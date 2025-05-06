import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file_path1 = "../NEURON_traces/AMPA_GrC.csv"
    df = pd.read_csv(file_path1, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_AMPA"] = df["Current"].values * 1000 / (-40.)

    plt.figure()
    plt.title('AMPA current (NEURON sim)')
    plt.plot(df["Time"], df['Current'].values*1000, color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('AMPA current [pA]')
    #plt.xlim(200,1000)
    plt.savefig('./figs/AMPA_current_neuron.png', dpi=150)
    plt.show()

    plt.figure()
    plt.title('AMPA conductance (NEURON sim)')
    plt.plot(df["Time"], df['g_AMPA'], color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('AMPA conductance [nS]')
    # plt.xlim(200,1000)
    plt.savefig('./figs/AMPA_conductance_neuron.png', dpi=150)
    plt.show()


    file_path2 = "../NEURON_traces/NMDA_GrC.csv"
    df = pd.read_csv(file_path2, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_NMDA"] = df["Current"].values * 1000 / ((-40.+3.7)  * 0.177) # 0.177 is the mg_block value associated to -40 mV,
                                                                       # -3.7 is the reversal potential for NMDA  

    plt.figure()
    plt.title('NMDA current (NEURON sim)')
    plt.plot(df["Time"], df['Current'].values*1000, color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('NMDA current [pA]')
    #plt.xlim(200,1000)
    plt.savefig('./figs/NMDA_current_neuron.png', dpi=150)
    plt.show()


    plt.figure()
    plt.title('NMDA conductance (NEURON sim)')
    plt.plot(df["Time"], df['g_NMDA'], color='black')
    plt.xlabel('Time [ms]')
    plt.ylabel('NMDA conductance [nS]')
    #plt.xlim(200,1000)
    plt.savefig('./figs/NMDA_conductace_neuron.png', dpi=150)
    plt.show()
