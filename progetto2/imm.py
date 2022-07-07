
#%%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\Sedov010.dat',header=None,sep='\s+') 
data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\Sedov000.dat',header=None,sep='\s+') 

#print(data.info)
# %%

colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.plot(data2[0],data2[6], label='total energy', color='#cb4679')
plt.plot(data2[0],data2[5], label='kinetic energy', color ='#a82296')
plt.plot(data2[0],data2[4], label='thermal energy', color ='#f89441')
#plt.plot(data[0],data[7], label='radiative energy', color ='black')
plt.hlines(0.,xmin=0.,xmax=100000, color='black', linestyle='dotted')
plt.xlabel('time [yrs]')
plt.ylabel('energy [units of $10^{51}$ ergs]')
plt.xscale('log')

plt.legend()
plt.show()

plt.plot(data2[0],data2[2],color ='#a82296',label='Sedov solution')
plt.plot(data[0],data[1],color ='#f89441',label='our project')
plt.xlabel('time [yrs]')
plt.ylabel('supernova radius [pc]')

plt.legend()
plt.show()

plt.plot(data2[0],data2[3], color=colors[0], label ='with cooling')
plt.plot(data2[0], data2[3], color =colors[3], label= 'without cooling' )
plt.title('X ray Luminosity')
plt.xlabel('time [yrs]')
plt.ylabel('X-ray Luminosity [$erg.s^{-1}$]')
plt.yscale('log')
plt.legend()
plt.show()