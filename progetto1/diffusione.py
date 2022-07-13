
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov010.dat',header=None,sep='\s+') 
data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto1\diff_sorg.dat',header=None,sep='\s+',
 comment='#', engine='python')
#data =  pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto2\sedov000.dat',header=None,sep='\s+')
#data2 = pd.read_csv(r'C:\Users\matteo\Desktop\computational\test1\shock_radius_winds2.dat',header=None,sep='\s+')

labels = ['t = 2 Gyr', 't = 3 Gyr', 't = 5 Gyr', 't = 8 Gyr', 't = 10 Gyr']
colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

fig, ax = plt.subplots(1,1, figsize = (10,5))

for i in range(3,8):
    ax.plot(data[0], data[i], color = colors[i-3], label = labels[i-3])
    ax.set_xscale('log')
ax.plot(data[0], data[1], label='$Z_{fe}$ subcracted', color = 'black')
ax.plot(data[0], data[2], label='$Z_{fe}$ observed', linestyle = '--',color = 'black')
ax.set_xlabel('distance [Mpc]')
ax.set_ylabel('$Z_{fe} / Z_{sol}$')
plt.legend()
fig.suptitle('source and diffusion')
plt.show()