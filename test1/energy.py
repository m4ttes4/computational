

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\shock_radius_winds2.dat',header=None,sep='\s+')
data2 =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\energ010.dat',header=None,sep='\s+', comment='#')
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\energ000.dat',header=None,sep='\s+')


colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)



plt.plot(data[0], data[1], color=colors[0], label='thermal')
plt.plot(data[0], data[2], color=colors[3], label='kinetic')
plt.plot(data[0], data[3], color=colors[1],label='total')
#plt.plot(data[0], data[7], color=colors[4],label='bolometric')
plt.xscale('log')
plt.legend()
plt.xlabel('time [pc]')
plt.ylabel('log energy [units of $10^{51}$ erg]')
plt.hlines(0,0,1e5, color='black', linestyle='--')
#plt.plot(data[0], data[7], color=colors[5], label='x_energy')
plt.show()





plt.plot(data2[0], data2[1], color=colors[0], label='thermal')
plt.plot(data2[0], data2[2], color=colors[1], label='kinetic')
plt.plot(data2[0], data2[3], color=colors[3],label='total')
#plt.plot(data[0], data[7], color=colors[4],label='bolometric')
plt.xscale('log')
plt.legend()
plt.xlabel('time [pc]')
plt.ylabel('log energy [units of $10^{51}$ erg]')
plt.hlines(0,0,1e5, color='black', linestyle='--')
#plt.plot(data[0], data[7], color=colors[5], label='x_energy')
#plt.show()
plt.show()


bolometric = np.array(data[7])
kinetic = np.array(data[5])
thermal = np.array(data[4])
dio = np.array(data[6])

total  = dio + bolometric


#plt.plot(data[0], bolometric, label='with radiation losses', color=colors[4])