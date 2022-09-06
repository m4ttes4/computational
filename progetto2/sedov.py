

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']
nomefile = 'speed_winds.dat'

comparison = False

ism = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\Sedov000.dat', sep='\s+',
                    engine='python',header=None, comment='#')
#ism2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\Sedov010.dat', sep='\s+',
                    #engine='python',header=None, comment='#')

if comparison == True:
    him = pd.read_csv(r'C:\Users\matteo\Desktop\scriptVSC\progetto2\shock_radius_HIM.dat', sep='\s+',
                        engine='python',header=None, comment='#')

    cim = pd.read_csv(r'C:\Users\matteo\Desktop\scriptVSC\progetto2\shock_radius_CIM.dat', sep='\s+',
                        engine='python',header=None, comment='#')

fig, ax = plt.subplots(nrows=1,ncols=1, figsize =(10,6))
ax.plot(ism[0], ism[2], label='standard ISM', color=colors[2], linestyle='-')
plt.plot(ism[0], ism[1], label='our project ISM', color=colors[2],linestyle='--')

#ax.plot(ism2[0], ism2[2], label='with cooling', color=colors[4], linestyle='-')
#plt.plot(ism2[0], ism2[1], label='our project with cooling', color=colors[4],linestyle='--')

if comparison == True:
    ax.plot(him[0],him[2], label = 'HIM',color=colors[3])
    plt.plot(him[0], him[1], label='our project HIM', color=colors[3],linestyle='--')

    ax.plot(cim[0],cim[2], label = 'CIM',color=colors[0])
    plt.plot(cim[0], cim[1], label='our project CIM', color=colors[0],linestyle='--')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('shock radius [Pc]')
ax.set_xlabel('time [yr]')
#ax.get_xaxis().set_visible(False)
plt.title('Sedov solutions comparison')

plt.legend()
plt.show()