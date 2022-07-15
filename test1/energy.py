

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\shock_radius_winds2.dat',header=None,sep='\s+')
#data =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov000.dat',header=None,sep='\s+')
data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\shock_radius_winds.dat',header=None,sep='\s+')


colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.plot(data[0], data[1], color=colors[3], label='our project')
plt.plot(data[0], data[2], color=colors[0], label='Sedov')
plt.legend()
plt.xlabel('time [yr]')
plt.ylabel('supernova radius [pc]')
plt.yscale('log')
plt.xscale('log')

plt.show()

plt.plot(data[0], data[3])

plt.xscale('log')
plt.yscale('log')
plt.show()

plt.plot(data[0], data[4], color=colors[0], label='thermal')
plt.plot(data[0], data[5], color=colors[1], label='kinetic')
plt.plot(data[0], data[6], color=colors[3],label='total')
#plt.plot(data[0], data[7], color=colors[4],label='bolometric')
plt.legend()
plt.xlabel('time [pc]')
plt.ylabel('energy [units of $10^{51}$ erg]')
plt.hlines(0,0,1e5, color='black', linestyle='--')
#plt.plot(data[0], data[7], color=colors[5], label='x_energy')
plt.show()





bolometric = np.array(data2[7])
kinetic = np.array(data[5])
thermal = np.array(data[4])
dio = np.array(data[6])

total  = dio - bolometric

plt.plot(data[0], dio)

plt.show()