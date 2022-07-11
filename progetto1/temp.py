

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv(r'C:\Users\matteo\Desktop\computational\progetto1\temp_paper.dat',header=None,sep='\s+',
 comment='#', engine='python')

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

labels = ['t = 2 Gyr', 't = 3 Gyr', 't = 5 Gyr', 't = 8 Gyr', 't = 10 Gyr']
colors =['#f89441','#e56b5d','#cb4679','#a82296','#0c0887']

plt.plot(data[0], data[1], color=colors[0], label='our project')
plt.plot(data[0], data[2], color=colors[3], label='Rebusco et al.')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('distance [kpc]')
plt.ylabel('T [K]')
plt.show()