

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)


data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\initial_pressure.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data1 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\initial_density.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\initial_temp.dat', sep='\s+',
                    engine='python',header=None, comment='#')

data3 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\initial_velocity.dat', sep='\s+',
                    engine='python',header=None, comment='#')

colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']
labels =['t=$2.10^{4}$ yrs','t=$4.10^{4}$ yrs','t=$6.10^{4}$ yrs','t=$8.10^{4}$ yrs','t=$1.10^{5}$ yrs']


for i in range(5):
    plt.plot(data[0], data[i+1], color=colors[i], label=labels[i])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('distance [pc]')
    plt.ylabel('pressure [cgs]')
plt.title('Pressure')
plt.legend()
plt.show()

for i in range(5):
    plt.plot(data3[0], data3[i+1], color=colors[i], label=labels[i])
    plt.xscale('log')
    #plt.xscale('log')
    plt.xlabel('distance [pc]')
    plt.ylabel('velocity [$km.s^{-1}$]')
plt.title('Velocity')
plt.legend()
plt.show()

for i in range(5):
    plt.plot(data2[0], data2[i+1], color=colors[i], label=labels[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('distance [pc]')
    plt.ylabel('density [$cm^-{3}$]')
    
plt.title('Density')
plt.legend()
plt.show()

for i in range(5):
    plt.plot(data1[0], data1[i+1], color=colors[i], label=labels[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('distance [pc]')
    plt.ylabel('temperature [K]')
    
plt.title('Temperature')
plt.legend()
plt.show()