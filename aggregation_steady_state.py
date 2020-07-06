# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd
import numpy as np

path_impact_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/massflux_withbreakthrough_forMartin_v4_complete.csv"
path_da_data = "Z:/Saturated_flow/diffusion_transient/rates_ss.csv"
impact = pd.read_csv(path_impact_data, sep = "\t")
da = pd.read_csv(path_da_data, sep = "\t")

impact.columns
impact['%del2massflux'] = impact['del2massflux'] * 100
impact.Chem.unique()
gvarnames = ["DO", "Ammonium", "Nitrate"]
subset = impact[impact['Chem'].isin (gvarnames)]
gratenames = list(da.Rate_type.unique())
Pelist = [0.1, 1, 10]
Reglist = ["Slow", "Medium", "Fast"]

for i in range(len(Reglist)):
    subset.loc[subset.Regime == Reglist[i], 'Pe'] = Pelist[i]

df = pd.merge(
    subset,
    da[["Trial", "Regime", "Rate_type", "Totalrate", "Microbe","Rateperbio", "Da", "Dabio", "Chem"]],
    on=["Trial", "Regime", "Chem"]
)

df["darec"] = 1/df["Da"]

import matplotlib.pyplot as plt
mapping = ['o','^', 's']
plt.figure()
for sp in list(df.Chem.unique()):
    print (sp)
    dfc = df[df['Chem']==sp]
    m = mapping[list(df.Chem.unique()).index(sp)]
    print (m)
    plt.scatter("darec", "%del2massflux", s = 100, c = np.log(dfc["Pe"]), linewidths = 1, alpha = .7, edgecolor = 'k', cmap = "YlGnBu", marker = m, data = dfc, label = sp)
plt.xscale("log")
plt.yscale("log")
plt.xlim(left = 29000)
#plt.ylim(top = 1000)
plt.ylabel ("Log of percentage impact")
plt.xlabel ("Log of Da")
plt.legend()

subset['bin'] = pd.cut(subset['datime'], bins)
subsetgroup = subset[['bin','%del2massflux']].groupby('bin').median()
subsetgroup.plot(kind='bar')
dagroup = subset.groupby(['datime'])['%del2massflux'].describe()

#Descriptive statistice
group.Chem.describe

import seaborn as sns
sns.set_style("white")
sns.kdeplot(df['%ofhomogeneous'], df['%del2massflux'])
#sns.plt.show()
 
# Custom it with the same argument as 1D density plot
sns.kdeplot(df['darec'], df['%del2massflux'], cmap="Reds", shade=True, bw=.15)
plt.ylim([30,110])
plt.xlim([0,200])

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
