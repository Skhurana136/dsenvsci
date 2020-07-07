# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd
import numpy as np

path_impact_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/massflux_withbreakthrough_forMartin_v4_complete.csv"
path_da_data = "Z:/Saturated_flow/diffusion_transient/rates_median_ss.csv"
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

for Reg in list(df.Regime.unique()):
    for c in list(df.Chem.unique()):
        baseDa = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].Da.iloc[0]
        baseDabio = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].Dabio.iloc[0]
        basedarec = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].darec.iloc[0]
#        df = df[(df.Regime == Reg) & (df.Chem == c)]
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactDa'] = df[(df.Regime == Reg) & (df.Chem == c)].Da/baseDa
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactDabio'] = df[(df.Regime == Reg) & (df.Chem == c)].Dabio/baseDabio
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactdarec'] = df[(df.Regime == Reg) & (df.Chem == c)].darec/basedarec
        print(df.shape)

df.to_csv("Z:/Saturated_flow/diffusion_transient/median_da_impact.csv")

import matplotlib.pyplot as plt
mapping = ['o','^', 's']
plt.figure()
for sp in list(df.Chem.unique()):
    print (sp)
    dfc = df[df['Chem']==sp]
    m = mapping[list(df.Chem.unique()).index(sp)]
    print (m)
    plt.scatter("impactDa", "del2massflux", s = 100, c =np.log(dfc["Pe"]), linewidths = 1, alpha = .7, edgecolor = 'k', cmap = "YlGnBu", data = dfc, marker = m, label = sp)
#plt.xscale("log")
#plt.yscale("log")
#plt.xlim(left = 0.01)
#plt.ylim(bottom = 0.001)
#plt.axhline(0.1, linestyle = '--', color = 'gray')
#plt.axhline(1, linestyle = '--', color = 'gray')
#plt.axhline(10, linestyle = '--', color = 'gray')
#plt.axhline(100, linestyle = '--', color = 'gray')
#plt.axhline(1000, linestyle = '--', color = 'gray')
plt.ylabel ("Impact on Damkohler number")
plt.xlabel ("Spatial heterogeneity")
plt.legend()

#histogram

import seaborn as sns
subset = df[["fraction", "%del2massflux", "Da", "impactDa", "Regime", "Chem"]]
fig = sns.pairplot(subset, kind="scatter", hue="Chem", markers=["o", "s", "D"], palette="Set2")
fig.savefig("Z:/Saturated_flow/diffusion_transient/cycling-da.png", pad_inches = 0.1)

subset["DaLog"] = np.log(subset["Da"])
subset["%del2massfluxLog"] = np.log(subset["%del2massflux"])
subset2 = subset[["fraction", "%del2massfluxLog", "DaLog", "impactDa","Regime", "Chem"]]
# Basic correlogram
fig = sns.pairplot(subset2, kind="scatter", hue="Chem", markers=["o", "s", "D"], palette="Set2")
fig.savefig("Z:/Saturated_flow/diffusion_transient/cycling-da-log.png", pad_inches = 0.1)

df['bin'] = pd.cut(df['darec'], bins)
subsetgroup = df[['bin','%del2massflux']].groupby('bin').median()
subsetgroup.plot(kind='bar')
dagroup = df.groupby(['darec'])['%del2massflux'].describe()

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

subset = df.sort_values(by = 'fraction')
data = subset[['del2massflux', 'Pe', 'fraction']]
cont_data = pd.pivot_table(data, values = 'del2massflux',
                              index = 'Pe',
                              columns = 'fraction')

x = list(df.fraction.unique())
y = list(df.Pe.unique())
fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.contour3D (x, y, cont_data, cmap = 'binary')
