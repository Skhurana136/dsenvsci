# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:18:49 2020

@author: khurana
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""
import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import analyses.saturated_transient as sta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
domains = ["Original", "Half", "Double", "Big"]
domainodes = {"Original": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes
domaininvest = list(domainodes.keys())[:1]

vardict = proc.speciesdict("Saturated")
gvarnames = list(t for t in vardict.keys() if vardict[t]["Location"]=="Mobile") + ["Nitrogen", "TOC"]

#Sensitivity
row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domadd = domain + "_"
        else:
            domadd = ""
        basedata = np.load("D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_0/NS-AH/NS-AH_df.npy")
        for t in ["1", "2", "5"]:
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    amplitude, ampmax, amplitudebase, basemax = sta.conc_norm_amplitude(data, basedata, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], amplitudebase[gvarnames.index(g)], basemax[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Sensitivity", "Timloc_max", "Sensitivitybase", "Timloc_maxbase"])

tracerdata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Normalized_RMSamplitude_chem_withloc.csv", sep="\t")

import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import analyses.saturated_transient as sta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
domains = ["Original", "Half", "Double", "Big"]
domainodes = {"Original": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes
domaininvest = list(domainodes.keys())[:1]

vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))

#Sensitivity
row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domadd = domain + "_"
        else:
            domadd = ""
        benchmark = np.load("D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_0/NS-AH/NS-AH_df.npy")
        for t in ["1", "2", "5"]:
#            directory = "//tsclient/D/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    amplitude, ampmax, baseamp, basemax = sta.mass_norm_amplitude(data, benchmark, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], baseamp[gvarnames.index(g)], baseamp[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Sensitivity", "Timloc_max", "Sensitivitybase", "Timloc_maxbase"])

tracerdata = pd.read_csv("Z:/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv("Z:/Normalized_RMSamplitude_biomass_withtimeloc.csv", sep="\t")