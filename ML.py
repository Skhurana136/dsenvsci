# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:07:39 2020

@author: khurana
"""

import pandas as pd
import numpy as np

file = "Y:/Home/khurana/4. Publications/Paper4/N_cycling_microbes_activity_MK_2_SK.xlsx"

data = pd.read_excel(file)

data.head()
data.tail()
data.isnull().sum()
data = data.replace(["below_detection"], 1000)

filt = data[data['Rates_nmol_NOx_L_d'] != "not_measured"]
filt.info()

import matplotlib.pyplot as plt
plt.scatter("NH4_µmol_L","Rates_nmol_NOx_L_d", data = filt)
plt.scatter("NO2_µmol_L","Rates_nmol_NOx_L_d", data = filt)
plt.scatter("NO3_µmol_L","Rates_nmol_NOx_L_d", data = filt)

def pairplot (data, biomassindex):
    import seaborn as sns
    cols = data.columns.tolist()
    touse = [cols[biomassindex]] + cols[-7:-2]
    filt = data[touse]
    filt = filt[filt != "no_data"]
    filt[cols[biomassindex]] = np.log10(filt[cols[biomassindex]])
#    graph = filt.drop(axis = 1, columns = cols[biomassindex])
    sns.pairplot(data = filt)
    
    return None

pairplot (data, 3)
pairplot (data, 4)
pairplot (data, 5)
pairplot (data, 6)

plt.scatter(cols[3],"Rates_nmol_NOx_L_d", data = filt, label = cols[3])
plt.scatter(cols[4],"Rates_nmol_NOx_L_d", data = filt, label = cols[4])
plt.scatter(cols[5],"Rates_nmol_NOx_L_d", data = filt, label = cols[5])
#plt.scatter(cols[6],"Rates_nmol_NOx_L_d", data = filt, label = cols[6])
plt.legend()
plt.xscale("log")
plt.xlabel("Biomass")
plt.ylabel("Rate of NOx evolution_nmol/L/day")