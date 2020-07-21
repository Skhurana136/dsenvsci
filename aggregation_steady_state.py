# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plots.saturated_steady_state as sssp
import plots.general_plots as gp

path_da_data = "Z:/Saturated_flow/diffusion_transient/Da_mean_ss.csv"
da = pd.read_csv(path_da_data, sep = "\t")
da.columns
da.shape

Pelist = [2, 11, 22]
Reglist = ["Slow", "Medium", "Fast"]

for i in range(len(Reglist)):
    da.loc[da.Regime == Reglist[i], 'Pe'] = Pelist[i]

da["%del2massflux"] = da["del2massflux"]*100
da['PeDa'] = da['V_Da'] / da['Pe'] #this transformation is promising

finaldata = da

#Show distribution of Da numbers differentiated with Pe numbers investigated
#Kerneldensityplot

sns.distplot(finaldata[finaldata["Regime"] == "Slow"]["V_Da"], color = "indianred", label = "Slow flow", kde = False, bins = 10)
sns.distplot(finaldata[finaldata["Regime"] == "Medium"]["V_Da"], color = "g", label = "Medium flow", kde = False, bins = 10)
sns.distplot(finaldata[finaldata["Regime"] == "Fast"]["V_Da"], color = "steelblue", label = "Fast flow", kde = False, bins = 10)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.xlabel ("Damkohler number", fontsize = 15)
plt.ylabel ("Number of scenarios", fontsize = 15)
plt.title ("Distribution of Damkohler numbers", fontsize = 15)
#plt.xscale ("log")
plt.legend(fontsize = 12)

#Show variation of Da with spatial heterogeneity

sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "V_Da", "Da", "Logy")
plt.ylim(bottom = 0.01)
plt.title("Da# with spatial heterogeneity")

sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

sssp.scatter_chem_regime (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

#Generic plots without distinguishing between chemicals and flow regimes
gp.scatter(finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

gp.scatter (finaldata, "fraction", "Spatial heterogeneity", "V_Da", "Da", "Logy")
plt.ylim(bottom = 0.05)
plt.title("Da# (normalized by biomass) with spatial heterogeneity")

gp.scatter (finaldata, "V_Da", "Da", "del2massflux", "Impact on O, N removal", "Logx")
plt.xlim(left = 0.1)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

gp.scatter (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with impact on normalized Da")