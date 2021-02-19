# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:24:47 2020

@author: khurana
"""
import pandas as pd
import numpy as np
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

velocity = [0.00038, 0.0038, 0.038]
Pelist = [2, 11, 22]

for pe, vel, Reg in zip(Pelist, velocity,Regimes):
        chemdata.loc[(chemdata.Regime == Reg), 'Velocity'] = vel
        chemdata.loc[chemdata.Regime == Reg, 'Pe'] = pe
        for c in chemdata.Chem.unique().tolist():
            base = chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c) & (chemdata.Trial == 'H')]['reldelmassflux'].values[0]
#            print(base)
            chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c), 'base'] = base

chemdata['fraction_rel_delmf'] = chemdata.reldelmassflux/chemdata.base
chemdata['avgconc_in'] = chemdata.Inlet_total_mass_flux*0.2/(chemdata.Velocity*0.3)
chemdata['avgconc_out'] = chemdata.Outlet_mass_flux*0.2/(chemdata.Velocity*0.3)
chemdata['removal'] = (chemdata['avgconc_in'] - chemdata['avgconc_out']).abs()
chemdata['removal_rate'] = 100*chemdata['removal']/(chemdata['avgconc_in'] * chemdata['Breakthroughtime'])
chemdata['conc_Da'] = chemdata['removal_rate']*chemdata['Breakthroughtime']
#chemdata['conc_Da'] = chemdata['removal_rate']/np.log(0.05) #For first order Da estimates
chemdata["%reldelmassflux"] = chemdata["reldelmassflux"]*100
chemdata["%del2massflux"] = chemdata["del2massflux"]*100
chemdata["%fraction_rel_delmf"] = chemdata["fraction_rel_delmf"]*100
chemdata['PeDa'] = chemdata['conc_Da'] / chemdata['Pe'] #this transformation is promising

for Reg in Regimes:
    for c in chemdata.Chem.unique().tolist():
        base1 = chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c) & (chemdata.Trial == 'H')]['conc_Da'].values[0]
        chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c), 'Dabase'] = base1
        base2 = chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c) & (chemdata.Trial == 'H')]['PeDa'].values[0]
        chemdata.loc[(chemdata.Regime == Reg) & (chemdata.Chem == c), 'PeDabase'] = base2

chemdata['fraction_da'] = chemdata.conc_Da/chemdata.Dabase
chemdata['fraction_PeDa'] = chemdata.PeDa/chemdata.PeDabase
chemdata["%fraction_da"] = chemdata["fraction_da"]*100
chemdata['%fraction_PeDa'] = chemdata.fraction_PeDa*100
chemdata.to_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Norm_Conc_da_ss.csv", sep = "\t")

gvarnames = ["DO", "Ammonium", "Nitrate"]
chemsubset = chemdata[chemdata['Chem'].isin (gvarnames)]
chemsubset.shape

rate_path_data = r"Z:\Saturated_flow\diffusion_transient\Complex_Mean_rate_ss.csv"
ratedata = pd.read_csv(rate_path_data, sep = '\t')
ratedata.columns
ratedata.Rate_type.unique()
gratenames = ["DO removal", "Ammonium removal", "Nitrate removal"]

ratesubset = ratedata[ratedata['Rate_type'].isin (gratenames)]
velocity = [0.00038, 0.0038, 0.038]

for vel, Reg in zip(velocity,Regimes):
        ratesubset.loc[(ratesubset.Regime == Reg), 'Velocity'] = vel

comb = pd.merge(chemsubset, ratesubset[["Regime", "Trial", "Rate_type", "Volumetric_Meanrate", "Arithmetic_Meanrate","Chem", "Velocity"]], on = ["Regime", "Trial", "Chem"])

comb['V_Da'] = comb['Breakthroughtime']*comb['Volumetric_Meanrate']
comb['A_Da'] = comb['Breakthroughtime']*comb['Arithmetic_Meanrate']

avgconcdata = pd.read_csv("Z:/Saturated_flow/diffusion_transient/avgconc.csv", sep = ',')
avgconcdata['Conc_Rate'] = -1*avgconcdata['influx']*(avgconcdata['Inlet_conc'] - avgconcdata['Outlet_conc'])/(0.3*0.5*0.2)

comb = pd.merge(comb, avgconcdata[["Regime", "Trial", "Chem", "Inlet_conc", "Outlet_conc", "Conc_Rate"]], on = ["Regime", "Trial", "Chem"])
comb['Massflux_MeanRate'] = comb['delmassflux']/(0.3*0.5*0.2)
comb['MF_Da'] = comb['delmassflux']/(comb['Velocity']*0.3)
comb['Conc_Da'] = comb['Conc_Rate']*comb['Breakthroughtime']

comb['MF/V_Da'] = comb['MF_Da']/comb['V_Da']
comb['A/V_Da'] = comb['A_Da']/comb['V_Da']
comb['MF/V_R'] = comb['Massflux_MeanRate']/comb['Volumetric_Meanrate']
comb['A/V_R'] = comb['Arithmetic_Meanrate']/comb['Volumetric_Meanrate']
comb['C/V_Da'] = comb['Conc_Da']/comb['V_Da']
comb['C/V_R'] = comb['Conc_Rate']/comb['Volumetric_Meanrate']

comb.to_csv(r"Z:\Saturated_flow\diffusion_transient\Complex_Da_mean_ss.csv", sep = '\t')

#Potential further data transformations not implemented right now: #No need to implement this i think so commenting it out
for Reg in list(comb.Regime.unique()):
    for c in list(comb.Chem.unique()):
        baseVDa = comb[(comb.Regime == Reg) & (comb.Trial == 'H') & (comb.Chem == c)].V_Da.iloc[0]
        comb.loc[(comb.Regime == Reg) & (comb.Chem == c), 'impactVDa'] = comb[(comb.Regime == Reg) & (comb.Chem == c)].V_Da/baseVDa
        print(comb.shape)