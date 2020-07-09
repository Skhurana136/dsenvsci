# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd
import numpy as np

path_impact_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/massflux_withbreakthrough_forMartin_v4_complete.csv"
path_da_data = "Z:/Saturated_flow/diffusion_transient/Da_mean_ss.csv"
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

da.columns

df = pd.merge(
    subset,
    da[["Trial", "Regime", "Rate_type", "Meanrate", "Microbe", "Da", "Dabio", "Chem"]],
    on=["Trial", "Regime", "Chem"]
)


df["darec"] = 1/df["Da"]

for Reg in list(df.Regime.unique()):
    for c in list(df.Chem.unique()):
        baseDa = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].Da.iloc[0]
        baseDabio = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].Dabio.iloc[0]
        basedarec = df[(df.Regime == Reg) & (df.Trial == 'H') & (df.Chem == c)].darec.iloc[0]
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactDa'] = df[(df.Regime == Reg) & (df.Chem == c)].Da/baseDa
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactDabio'] = df[(df.Regime == Reg) & (df.Chem == c)].Dabio/baseDabio
        df.loc[(df.Regime == Reg) & (df.Chem == c), 'impactdarec'] = df[(df.Regime == Reg) & (df.Chem == c)].darec/basedarec
        print(df.shape)

df.to_csv("Z:/Saturated_flow/diffusion_transient/mean_da_biomassconc_impact.csv", sep = "\t")

finaldata = pd.read_csv("Z:/Saturated_flow/diffusion_transient/mean_da_biomassconc_impact.csv", sep = "\t") 
import plots.saturated_steady_state as sssp
import matplotlib.pyplot as plt
#Show variation of Da with spatial heterogeneity
ax1 = plt.subplot(321)
ax1 = sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "Dabio", "Normalized Da", "Logy")
plt.ylim(bottom = 0.01)
plt.title("Da# (normalized by biomass) with spatial heterogeneity")

ax2 = plt.subplot(322)
ax2 = sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "impactDabio", "Impact on normalized Da", "None")
plt.xlim(left = 0.01)
plt.title("Impact on Da (normalized by that in base case) with spatial heterogeneity")

ax3 = plt.subplot(323)
ax3 = sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "%del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

ax4 = plt.subplot(324)
ax4 = sssp.scatter_chem_regime (finaldata, "Dabio", "Normalized Da", "%del2massflux", "Impact on O, N removal", "Logx")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with normalized Da")

ax5 = plt.subplot(325)
ax5 = sssp.scatter_chem_regime (finaldata, "impactDabio", "Impact on normalized Da", "%del2massflux", "Impact on O, N removal (%)", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with impact on normalized Da")
  
sssp.scatter_chem_regime (finaldata, "Dabio", "Normalized Da", "%del2massflux", "Normalized impact on removal (%)", "None")
plt.xlim(left = 0.01)

import plots.general_plots as gp

gp.scatter(finaldata, "fraction", "Spatial heterogeneity", "%del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

gp.scatter (finaldata, "fraction", "Spatial heterogeneity", "Dabio", "Normalized Da", "Logy")
plt.ylim(bottom = 0.01)
plt.title("Da# (normalized by biomass) with spatial heterogeneity")

gp.scatter (finaldata, "fraction", "Spatial heterogeneity", "impactDabio", "Impact on normalized Da", "None")
plt.xlim(left = 0.01)
plt.title("Impact on Da (normalized by that in base case) with spatial heterogeneity")

gp.scatter (finaldata, "Dabio", "Normalized Da", "%del2massflux", "Impact on O, N removal", "Logx")
plt.xlim([0.001, 1000])
plt.title("Impact on chemical removal (normalized by that in base case) with normalized Da")

gp.scatter (finaldata, "impactDabio", "Impact on normalized Da", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with impact on normalized Da")

#histogram

import seaborn as sns
subset = df[["fraction", "%del2massflux", "Dabio", "impactDabio", "Regime", "Chem"]]
fig = sns.pairplot(subset, kind="scatter", hue="Regime", markers=["o", "s", "D"], palette="Set2")
fig.savefig("Z:/Saturated_flow/diffusion_transient/cycling-da.png", pad_inches = 0.1)

subset["DaLog"] = np.log(subset["Dabio"])
subset["%del2massfluxLog"] = np.log(subset["%del2massflux"])
subset2 = subset[["fraction", "%del2massfluxLog", "DaLog", "impactDabio","Regime", "Chem"]]
# Basic correlogram
fig = sns.pairplot(subset2, kind="scatter", hue="Chem", markers=["o", "s", "D"], palette="Set2")
fig.savefig("Z:/Saturated_flow/diffusion_transient/cycling-da-log.png", pad_inches = 0.1)

df['bin'] = pd.cut(df['darec'], bins)
subsetgroup = df[['bin','%del2massflux']].groupby('bin').median()
subsetgroup.plot(kind='bar')
dagroup = df.groupby(['darec'])['%del2massflux'].describe()

#Descriptive statistice
group.Chem.describe