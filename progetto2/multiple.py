

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)


data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sus.dat', sep='\s+',
                    engine='python',header=None, comment='#')


colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']
labels =['t=$2.10^{4}$ yrs','t=$4.10^{4}$ yrs','t=$6.10^{4}$ yrs','t=$8.10^{4}$ yrs','t=$1.10^{5}$ yrs']

fig,ax = plt.subplots(2,2 ,figsize = (14,10))

ax[0,0].plot(data[0], data[2], color=colors[4])
#ax[0,0].set_xscale('log')
#ax[0,0].set_yscale('log')
ax[0,0].set_xlabel('distance')
ax[0,0].set_ylabel('density')
ax[0,0].set_title('density')

ax[0,1].plot(data[0], data[3], color=colors[4])
#ax[0,1].set_xscale('log')
#ax[0,1].set_yscale('log')
ax[0,1].set_xlabel('distance')
ax[0,1].set_ylabel('velocity')
ax[0,1].set_title('velocity')

ax[1,0].plot(data[0], data[4], color=colors[4])
#ax[1,0].set_xscale('log')
#ax[1,0].set_yscale('log')
ax[1,0].set_xlabel('distance')
ax[1,0].set_ylabel('energy')
ax[1,0].set_title('specific energy')

ax[1,1].plot(data[0], data[5], color=colors[4])
#ax[1,1].set_xscale('log')
#ax[1,1].set_yscale('log')
ax[1,1].set_xlabel('distance')
ax[1,1].set_ylabel('pressure')
ax[1,1].set_title('pressure')

fig.suptitle('shock tube test')
plt.show()