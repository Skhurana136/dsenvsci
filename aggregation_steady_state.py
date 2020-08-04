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
my_pal = {2:"indianred", 11:"g", 22:"steelblue", "DO":"indianred", "Nitrate":"g", "Ammonium":"steelblue",
          "Slow":"indianred", "Medium":"g", "Fast":"steelblue"}

path_da_data = "Z:/Saturated_flow/diffusion_transient/Complex_Da_mean_ss.csv"
da = pd.read_csv(path_da_data, sep = "\t")
da.columns
da.shape

Pelist = [2, 11, 22]
Reglist = ["Slow", "Medium", "Fast"]

for i in range(len(Reglist)):
    da.loc[da.Regime == Reglist[i], 'Pe'] = Pelist[i]

da["%reldelmassflux"] = da["reldelmassflux"]*100
da["%del2massflux"] = da["del2massflux"]*100
da['PeDa'] = da['V_Da'] * da['Pe'] #this transformation is promising

finaldata = da

#Show distribution of Da numbers differentiated with Pe numbers investigated
#Kerneldensityplot

plotdistribution("MF/V_R")
plotdistribution("A/V_R")
plotboxplot(finaldata, "MF/V_R", "Regime")
plotboxplot(finaldata, "C/V_R", "Regime")

def plotboxplot(finaldata, variable, category):
    my_pal = {"Slow":"indianred", "Medium":"g", "Fast":"steelblue", "DO":"indianred", "Nitrate":"g", "Ammonium":"steelblue"}
    
    dummy = sns.boxplot( x=finaldata[category], y=finaldata[variable], palette = my_pal)

    medians = finaldata.groupby([category])[variable].median().values
    nobs = finaldata[category].value_counts().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]
 
    # Add it to the plot
    pos = range(len(nobs))
    for tick,label in zip(pos,dummy.get_xticklabels()):
        dummy.text(pos[tick], medians[tick] + 0.03, nobs[tick], horizontalalignment='center', size=12, color='w', weight='semibold')

    plt.ylabel("Ratio", fontsize = 15)
    plt.xlabel ("Reactive species", fontsize = 15)

def plotdistribution(variable):
    sns.distplot(finaldata[finaldata["Regime"] == "Slow"][variable], color = "indianred", label = "Slow flow", kde = False, bins = 10)
    sns.distplot(finaldata[finaldata["Regime"] == "Medium"][variable], color = "g", label = "Medium flow", kde = False, bins = 10)
    sns.distplot(finaldata[finaldata["Regime"] == "Fast"][variable], color = "steelblue", label = "Fast flow", kde = False, bins = 10)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.xlabel ("Damkohler number", fontsize = 15)
    plt.ylabel ("Number of scenarios", fontsize = 15)
    plt.title ("Distribution of Damkohler numbers", fontsize = 15)
    #plt.xscale ("log")
    plt.legend(fontsize = 12)

#Show variation of Da with spatial heterogeneity

sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "MF_Da", "Da", "Logy")
#plt.ylim(bottom = 0.05)
plt.ylim(bottom = 10)
plt.title("Da# (normalized by biomass) with spatial heterogeneity")

sssp.scatter_chem_regime (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

sssp.scatter_chem_regime (finaldata, "MF_Da", "Da", "del2massflux", "Impact on O, N removal", "Logx")
#plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

sssp.scatter_chem_regime (finaldata, "PeDa", "V_Daratio", "MF_Da", "MF_Da", "True")
plt.xlim(left = 1)
plt.ylim(bottom = 10)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "Conc_Da", "C_Da", "Logy")
#plt.ylim(bottom = 0.05)
#plt.ylim(bottom = 10)
plt.title("C_Da# (normalized by biomass) with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "V_Da", "V_Da", "Conc_Da", "C_Da", "True")
plt.xlim(left = 1)
plt.ylim(bottom = 0.1)
plt.title("Conncentration derived Da with volumetrically averaged Da")          

sns.boxplot(data = finaldata, x = finaldata['Chem'], y = finaldata['Conc_Da']/finaldata['V_Da'], hue = "Regime", palette = my_pal)
plt.xlim(left = 1)
plt.ylim(bottom = 0.1)
plt.title("Conncentration derived Da with volumetrically averaged Da")          

sssp.scatter_chem_regime (finaldata, "Conc_Da", "Concentration derived Da", "%reldelmassflux", "O, N removal", "Logx")
plt.xlim(left = 0.1)
plt.title("Reactive species removal with Damkohler numbers")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "Conc_Da", "C_Da", "Logy")
plt.ylim(bottom = 0.1)
plt.title("Concentration derived Da# with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "Conc_Da", "Conc_Da", "%reldelmassflux", "O,N removal", "Logy")
plt.ylim(bottom = 0.1)
plt.title("Concentration derived Da# with spatial heterogeneity")

#Generic plots without distinguishing between chemicals and flow regimes
gp.scatter(finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

gp.scatter (finaldata, "fraction", "Spatial heterogeneity", "MF_Da", "Da", "Logy")
plt.ylim(bottom = 0.05)
plt.title("Da# (normalized by biomass) with spatial heterogeneity")

gp.scatter (finaldata, "V_Da", "Da", "del2massflux", "Impact on O, N removal", "Logx")
plt.xlim(left = 0.1)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

gp.scatter (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with impact on normalized Da")

gp.scatter (finaldata, "Conc_Da", "Concentration derived Da", "%reldelmassflux", "O, N removal", "Logx")
plt.xlim(left = 0.1)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")
