

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto1\temp_paper.dat',header=None,sep='\s+',
 comment='#', engine='python')
data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto1\massaDM.dat',header=None,sep='\s+',
 comment='#', engine='python')
data3 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto1\densit.dat',header=None,sep='\s+',
 comment='#', engine='python')
#data =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov000.dat',header=None,sep='\s+')
#data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\shock_radius_winds2.dat',header=None,sep='\s+')
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
labels = ['t = 2 Gyr', 't = 3 Gyr', 't = 5 Gyr', 't = 8 Gyr', 't = 10 Gyr']
colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']


plt.plot(data3[0], data3[1],color = 'black',linestyle='--', label='dark matter density')
plt.plot(data[0], data[3],color = colors[0], label='Rebusco et al.')
plt.plot(data[0], data[4],color = colors[3], label='our project')
plt.xlabel('distance [kpc]')
plt.ylabel('density [cgs]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()


plt.plot(data2[0], data2[1],color = colors[3], label='dark matter numerical')
plt.plot(data2[0], data2[2],color = colors[2], label='dark matter analitical')
plt.plot(data2[0], data2[3],color = colors[0], label='stars mass profile')
plt.xlabel('distance [kpc]')
plt.ylabel('density [cgs]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()