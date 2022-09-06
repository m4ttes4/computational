import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np




data2 =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\relazioni\sedov000.dat',header=None,sep='\s+')


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

plt.plot(data2[0], data2[5], color=colors[3], label='project, with cooling')



compa=True

if compa == True:
    
    data =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\relazioni\sedov001.dat',header=None,sep='\s+')
    plt.plot(data[0], data[5], color=colors[0], label='project, without cooling')
    plt.plot(data[0], data[4], color='black',linestyle='--', label='Theoretical T [K]')

plt.legend()
plt.xlabel('time [yr]')
plt.ylabel('supernova radius [pc]')
plt.yscale('log')
plt.xscale('log')
plt.title('Post-shock temperature')
plt.show()

