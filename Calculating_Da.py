# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:24:47 2020

@author: khurana
"""
import pandas as pd
import data_reader.data_processing as proc
import numpy as np
import data_reader.reader as rdr

# Saturated flow regime
Regimes = ["Slow", "Equal", "Fast"]
fpre = "NS-A"
fsuf = r"/"
gw = 1
filename = "ratesAtFinish.dat"

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

comb = pd.merge(chemsubset, ratesubset[["Regime", "Trial", "Rate_type", "Meanrate", "Chem"]], on = ["Regime", "Trial", "Chem"])

comb['Da'] = comb['Meanrate']/comb['Breakthroughtime']

#row = []
#for c, r, in zip(gvarnames, gratenames)):
#    chemrow = chemsubset.loc[(chemsubset.Regime == Reg) & (chemsubset.Trial == t) & (chemsubset.Chem == bioname)]
#    traveltime = biorow.Breakthroughtime.iloc[0]
#    da = traveltime*meanrate
#    dabio = traveltime*meanrate/meanbiomass
#    ka = traveltime*meanrate
#    kabio = traveltime*meanrate/meanbiomass
#    row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate, bioname, meanbiomass, meanrate/meanbiomass, da, ka, dabio, kabio, c])
#
#row = []
#for Reg in Regimes:
#    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
#    if Reg == "Equal":
#        Reg = "Medium"
#    for t in Trial:
#        print(t)
#        filepath = directory + fpre + str(t) + fsuf + filename
#        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=16 + respindx)
#        df = np.load(directory + fpre + t + fsuf + fpre + t + "_df.npy")
#        for bioindx, bioname, i,c in zip(microvars, microbes, range(len(respindx)), gvarnames):
#            biorow = biomass.loc[(biomass.Regime == Reg) & (biomass.Trial == t) & (biomass.Chem == bioname)]
#            ratesdf = rdr.Converttoarray_1581(M, "rates")
#            meanbiomass = np.mean(df[bioindx, - 1, :, :])
#            traveltime = biorow.Breakthroughtime.iloc[0]
#            da = traveltime*meanrate
#            dabio = traveltime*meanrate/meanbiomass
#            ka = traveltime*meanrate
#            kabio = traveltime*meanrate/meanbiomass
#            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate, bioname, meanbiomass, meanrate/meanbiomass, da, ka, dabio, kabio, c])

#row = []
#for Reg in Regimes:
#    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
#    if Reg == "Equal":
#        Reg = "Medium"
#    for t in Trial:
#        print(t)
#        filepath = directory + fpre + str(t) + fsuf + filename
#        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=16 + respindx)
#        df = np.load(directory + fpre + t + fsuf + fpre + t + "_df.npy")
#        for bioindx, bioname, i,c in zip(microvars, microbes, range(len(respindx)), gvarnames):
##            biorow = biomass.loc[(biomass.Regime == Reg) & (biomass.Trial == t) & (biomass.Chem == bioname)]
#            meanrate = np.mean(M[:, i])
#            meanbiomass = np.mean(df[bioindx, - 1, :, :])
#            traveltime = biorow.Breakthroughtime.iloc[0]
#            da = traveltime*meanrate
#            dabio = traveltime*meanrate/meanbiomass
#            ka = traveltime*meanrate
#            kabio = traveltime*meanrate/meanbiomass
#            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate, bioname, meanbiomass, meanrate/meanbiomass, da, ka, dabio, kabio, c])

#df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Rate_type', 'Meanrate', 'Microbe', 'Meanbiomassconc', 'Rateperbio', 'Da', 'Ka', 'Dabio', 'Kabio', 'Chem'])


#df.to_csv(r"Z:\Saturated_flow\diffusion_transient\Da_mean_ss.csv", sep = '\t')