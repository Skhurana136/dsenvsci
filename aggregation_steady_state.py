# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd

path_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/massflux_withbreakthrough_forMartin_v4_complete.csv"
master = pd.read_csv(path_data, sep = "\t")

master.columns
master['%del2massflux'] = master['del2massflux'] * 100
master.Chem.unique()
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate"]
subset = master[master['Chem'].isin (gvarnames)]

Pelist = [0.1, 1, 10]
Reglist = ["Slow", "Medium", "Fast"]
Dalist = [100, 100, 0.1, 10]

for i in range(len(Reglist)):
    subset.loc[subset.Regime == Reglist[i], 'Pe'] = Pelist[i]
    for k in range(len(gvarnames)):
        subset.loc[(subset.Regime == Reglist[i]) & (subset.Chem == gvarnames[k]), 'Da'] = Dalist[k]/Pelist[i]

x = subset.groupby(['Regime', 'Variance', 'Anisotropy', 'Chem'])['%del2massflux'].describe()

subset['datime'] = subset['Da']*subset['fraction']
bins = np.linspace (min(subset['datime']), max(subset['datime']), 1000)

subset['bin'] = pd.cut(subset['datime'], bins)
subsetgroup = subset[['bin','%del2massflux']].groupby('bin').median()
subsetgroup.plot(kind='bar')
dagroup = subset.groupby(['datime'])['%del2massflux'].describe()

import matplotlib.pyplot as plt
plt.scatter(subset['datime'], subset['Da'], s=subset['del2massflux']*100, alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.ylim([0.002,1200])
plt.xlim([0.002,1200])
plt.ylabel("Damkohler number")
plt.xlabel("Damkohler and time produce")

plt.scatter(subset['Da'], subset['Pe'], s=subset['del2massflux']*1000, alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.ylim([0.002,12])
plt.xlim([0.002,1200])
plt.ylabel("Peclet number")
plt.xlabel("Damkohler")

plt.scatter(subset['datime'], subset['Chem'], s=subset['del2massflux']*1000, alpha=0.5)
#plt.yscale('log')
plt.xscale('log')
#plt.ylim([0.002,12])
plt.xlim([0.002,1200])
plt.ylabel("Peclet number")
plt.xlabel("Damkohler and time product")

from pandas.tools.plotting import parallel_coordinates
plt.figure()
parallel_coordinates(subset, 'Chem', colormap=plt.get_cmap("Set2"))


plt.show()

#Descriptive statistice
group.Chem.describe

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
sns.kdeplot(subset['%ofhomogeneous'], subset['%del2massflux'])
#sns.plt.show()
 
# Custom it with the same argument as 1D density plot
sns.kdeplot(subset['Da'], subset['%del2massflux'], cmap="Reds", shade=True, bw=.15)
plt.ylim([30,110])
plt.xlim([0,200])

 
# Some features are characteristic of 2D: color palette and wether or not color the lowest range
sns.kdeplot(df.sepal_width, df.sepal_length, cmap="Blues", shade=True, shade_lowest=True, )
sns.plt.show()


from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

subset = subset.sort_values(by = subset['fraction'])
cont_data = pd.pivot_table(subset, values = 'del2massflux',
                              index = ['fraction'],
                              columns = 'Pe')

x = subset.Pe.unique()
y = subset.fraction.unique()
fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.contour3D (x, y, cont_data, cmap = 'binary')
