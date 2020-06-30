# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:04:50 2020

@author: khurana
"""
import pandas as pd
import numpy as np


def processdataframe(data, variablenames):
    data["Regime"] = data["Regime"].astype(int)
    data["Regime"] = data["Regime"].replace([0, 1, 2], ["Slow", "Medium", "Fast"])
    data["Regime"] = data["Regime"].astype(str)
    data["Trial"] = data["Trial"].astype(int)
    data["Trial"] = data["Trial"].replace([0], "H")
    data["Trial"] = data["Trial"].astype(str)
    data["Chem"] = data["Chem"].astype(int)
    for k in range(len(variablenames)):
        data["Chem"] = data["Chem"].replace([k], variablenames[k])
    data["Chem"] = data["Chem"].astype(str)

    return data


def processchembiomassfiles(biomassfile, regime):
    biomass = pd.read_csv(biomassfile, delimiter="\t")
    biomass["Chemb"] = biomass["Chem"]
    biomass["VA"] = biomass["Variance"] * biomass["Anisotropy"]
    #    chemfile = pd.read_csv(chemfiles, delimiter = "\t")
    #    oxy = pd.merge(biomass[biomass['Chemb']=='Aerobes'], chemfile[chemfile['Chem']=='DO'], on = 'Trial', suffixes=("_b","_c"))
    #    amm = pd.merge(biomass[biomass['Chemb']=='Ammonia oxidizers'], chemfile[chemfile['Chem']=='Ammonium'], on = 'Trial', suffixes=("_b","_c"))
    #    nitra = pd.merge(biomass[biomass['Chemb']=='Nitrate reducers'], chemfile[chemfile['Chem']=='Nitrate'], on = 'Trial', suffixes=("_b","_c"))
    #    allbiomass = pd.concat([oxy,amm,nitra], axis=0, ignore_index = True)
    #    biomass['delbiomass%']=biomass['TotalRatio']*100
    #    biomass['del2biomass%']=biomass['TotalRatio_Time']*100
    #    biomass['delspbiomass%']=biomass['SpeciesRatio']*100
    #    biomass['del2spbiomass%']=biomass['SpeciesRatio_Time']*100
    biomass["Regime"] = regime

    return biomass


def processchemfiles(chemfile, regime):
    chemfile = pd.read_csv(chemfile, delimiter="\t")
    chemfile["Regime"] = regime
    chemfile["VA"] = chemfile["Variance"] * chemfile["Anisotropy"]
    #    chemfile['RemovalRatio%']=chemfile['RemovalRatio']*100
    #    chemfile['RemovalRatio_Time%']=chemfile['RemovalRatio_Time']*100
    return chemfile


def tracerstudies(filename):
    combined_tracer = pd.read_csv(filename,
        delimiter="\t",
    )
    combined_tracer.loc[combined_tracer.Regime == "Equal", "Regime"] = "Medium"
    combined_tracer["fraction_withslow"] = (
        combined_tracer["Time"] / combined_tracer["Time"][0]
    )
    l = []
    for i in range(len(combined_tracer)):
        if combined_tracer["Variance"][i] == 0.1:
            l.append(
                str(combined_tracer["Variance"][i])
                + ":"
                + str(combined_tracer["Anisotropy"][i])
            )
        else:
            l.append(
                str(int(combined_tracer["Variance"][i]))
                + ":"
                + str(combined_tracer["Anisotropy"][i])
            )

    combined_tracer["Xlabels"] = l
    combined_tracer = combined_tracer.sort_values(by=["Variance", "Anisotropy"])

    return combined_tracer


def localmaxmin(y1, y2, init):
    ymax = np.max(np.max(np.append(y1, y2), axis=0))
    ymin = np.min(np.min(np.append(y1, y2), axis=0))
    xposmax = (np.argmax(y1) + init) * 5
    xposmin = (np.argmin(y1) + init) * 5
    print(ymax, ymin, xposmax, xposmin)
    return ymax, ymin, xposmax, xposmin

#master dictionary for the scenarios that have been run
def masterscenarios():
    #keys in this master dictionary refer to scenario names/subfolder names
    #Het refers to the variance in permeability enforced while generating the SRF
    #Anis refer to Anisotropy enforced while generating the SRF
    #Truncatedpowerlaw was used as the variogram model
    masterdictionary = {
            'H': {"Het": 0, "Anis": 1},
            '37': {"Het": 0.1, "Anis": 2},
            '38': {"Het": 0.1, "Anis": 5},
            '39': {"Het": 0.1, "Anis": 10},
            '40': {"Het": 1, "Anis": 2},
            '41': {"Het": 1, "Anis": 5},
            '42': {"Het": 1, "Anis": 10},
            '43': {"Het": 10, "Anis": 2},
            '44': {"Het": 10, "Anis": 5},
            '45': {"Het": 10, "Anis": 10},
            '46': {"Het": 0.1, "Anis": 2},
            '47': {"Het": 0.1, "Anis": 5},
            '48': {"Het": 0.1, "Anis": 10},
            '49': {"Het": 1, "Anis": 2},
            '50': {"Het": 1, "Anis": 5},
            '51': {"Het": 1, "Anis": 10},
            '52': {"Het": 10, "Anis": 2},
            '53': {"Het": 10, "Anis": 5},
            '54': {"Het": 10, "Anis": 10},
            '55': {"Het": 0.1, "Anis": 2},
            '56': {"Het": 0.1, "Anis": 5},
            '57': {"Het": 0.1, "Anis": 10},
            '58': {"Het": 1, "Anis": 2},
            '59': {"Het": 1, "Anis": 5},
            '60': {"Het": 1, "Anis": 10},
            '61': {"Het": 10, "Anis": 2},
            '62': {"Het": 10, "Anis": 5},
            '63': {"Het": 10, "Anis": 10},
            '64': {"Het": 0.1, "Anis": 2},
            '65': {"Het": 0.1, "Anis": 5},
            '66': {"Het": 0.1, "Anis": 10},
            '67': {"Het": 1, "Anis": 2},
            '68': {"Het": 1, "Anis": 5},
            '69': {"Het": 1, "Anis": 10},
            '70': {"Het": 5, "Anis": 2},
            '71': {"Het": 5, "Anis": 5},
            '72': {"Het": 5, "Anis": 10},
            '73': {"Het": 10, "Anis": 2},
            '74': {"Het": 10, "Anis": 5},
            '75': {"Het": 10, "Anis": 10},
            '76': {"Het": 5, "Anis": 2},
            '77': {"Het": 5, "Anis": 5},
            '78': {"Het": 5, "Anis": 10},
            '79': {"Het": 5, "Anis": 2},
            '80': {"Het": 5, "Anis": 5},
            '81': {"Het": 5, "Anis": 10},
            '82': {"Het": 5, "Anis": 2},
            '83': {"Het": 5, "Anis": 5},
            '84': {"Het": 5, "Anis": 10},
            }
    
    return masterdictionary