

#%%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)


data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\pressure_winds.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data1 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\temperature_winds.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\density_winds.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data3 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\speed_winds.dat', sep='\s+',
                    engine='python',header=None, comment='#')


#

#print(data.head)
colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']
labels =['t=$2.10^{4}$ yrs','t=$4.10^{4}$ yrs','t=$6.10^{4}$ yrs','t=$8.10^{4}$ yrs','t=$1.10^{5}$ yrs']

fig, ax = plt.subplots(nrows=2,ncols=2, figsize =(14,8))

for i in range(5):
    ax[0,0].plot(data[0], data[i+1], color=colors[i], label=labels[i])
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    ax[0,0].set_xlabel('distance [pc]')
    ax[0,0].set_ylabel('pressure [cgs]')

for i in range(5):
    ax[0,1].plot(data1[0], data1[i+1], color=colors[i])
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlabel('distance [pc]')
    ax[0,1].set_ylabel('temperature [K]')

for i in range(5):
    ax[1,0].plot(data2[0], data2[i+1], color=colors[i])
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    ax[1,0].set_xlabel('distance [pc]')
    ax[1,0].set_ylabel('density [$cm^{-3}$]')

for i in range(5):
    ax[1,1].plot(data3[0], data3[i+1], color=colors[i],label=labels[i])
    ax[1,1].set_xscale('log')
    #ax1.set_yscale('log')
    ax[1,1].set_xlabel('distance [pc]')
    ax[1,1].set_ylabel('speed [cgs]')



ax[1,1].legend()

plt.show()
