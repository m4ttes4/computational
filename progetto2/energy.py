

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\energ000.dat',header=None,sep='\s+')
data2 =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov000.dat',header=None,sep='\s+')



colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.plot(data[0], data[2])
plt.plot(data[0], data[3])
plt.plot(data[0], data[1])
plt.plot(data[0], data[4])
plt.xscale('log')

plt.show()

plt.plot(data2[0], data2[3])
plt.plot(data2[0], data2[4])
plt.xscale('log')
plt.yscale('log')
plt.show()