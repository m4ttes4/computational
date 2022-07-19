import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\sedov000.dat',header=None,sep='\s+')
data2 =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\sedov010.dat',header=None,sep='\s+')



colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.plot(data2[0], data2[1], color=colors[3], label='our project')
plt.plot(data2[0], data2[2], color=colors[0], label='Sedov')
plt.legend()
plt.xlabel('time [yr]')
plt.ylabel('supernova radius [pc]')
plt.yscale('log')
plt.xscale('log')

plt.show()

#plt.plot(data[0], data[1], color=colors[3], label='our project')
plt.plot(data[0], data[3], color=colors[0], label='Winds radius')
plt.legend()
plt.xlabel('time [yr]')
plt.ylabel('supernova radius [pc]')
plt.yscale('log')
plt.xscale('log')

plt.show()
