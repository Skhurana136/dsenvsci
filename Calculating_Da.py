# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:24:47 2020

@author: khurana
"""
import pandas as pd
import data_reader.data_processing as proc

# Saturated flow regime
Regimes = ["Slow", "Medium", "Fast"]
fpre = "NS-A"
fsuf = r"/"
gw = 1

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
ratenames = proc.masterrates("saturated")

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

# Calculate bulk Damkohler numbers in the domain

chem_path_data = r"Y:\Home\khurana\4. Publications\Restructuring\Paper1\Figurecodes\massflux_withbreakthrough_forMartin_v4_complete.csv"
chemdata = pd.read_csv(chem_path_data, sep = '\t')
chemdata.columns
chemdata.Chem.unique()
gvarnames = ["DO", "Ammonium", "Nitrate"]
chemsubset = chemdata[chemdata['Chem'].isin (gvarnames)]
chemsubset.shape

rate_path_data = r"Z:\Saturated_flow\diffusion_transient\Mean_rate_ss.csv"
ratedata = pd.read_csv(rate_path_data, sep = '\t')
ratedata.columns
ratedata.Rate_type.unique()
gratenames = ['Immobile aerobic respiration',
       'Immobile nitrate respiration',
       'Immobile ammonium respration', 
       ]

ratesubset = ratedata[ratedata['Rate_type'].isin (gratenames)]

for Reg in Regimes:
    for c, r, in zip(gvarnames, gratenames):
        ratesubset.loc[(ratesubset.Regime == Reg) & (ratesubset.Rate_type == r), 'Chem'] = c

comb = pd.merge(chemsubset, ratesubset[["Regime", "Trial", "Rate_type", "Volumetric_Meanrate", "Arithmetic_Meanrate","Chem"]], on = ["Regime", "Trial", "Chem"])

comb['V_Da'] = comb['Breakthroughtime']*comb['Volumetric_Meanrate']
comb['A_Da'] = comb['Breakthroughtime']*comb['Arithmetic_Meanrate']

comb.to_csv(r"Z:\Saturated_flow\diffusion_transient\Da_mean_ss.csv", sep = '\t')

#Potential further data transformations not implemented right now: #No need to implement this i think so commenting it out
#for Reg in list(comb.Regime.unique()):
#    for c in list(comb.Chem.unique()):
#        baseVDa = comb[(comb.Regime == Reg) & (comb.Trial == 'H') & (comb.Chem == c)].V_Da.iloc[0]
#        comb.loc[(comb.Regime == Reg) & (comb.Chem == c), 'impactVDa'] = comb[(comb.Regime == Reg) & (comb.Chem == c)].V_Da/baseVDa
#        baseADa = comb[(comb.Regime == Reg) & (comb.Trial == 'H') & (comb.Chem == c)].A_Da.iloc[0]
#        comb.loc[(comb.Regime == Reg) & (comb.Chem == c), 'impactADa'] = comb[(comb.Regime == Reg) & (comb.Chem == c)].A_Da/baseADa
#        print(comb.shape)