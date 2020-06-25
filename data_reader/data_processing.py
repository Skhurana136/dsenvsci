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


def tracerstudies():
    di = r"Z:\Saturated_flow\Steady_state\Tracer_studies\/"
    filename = "NS-ATracerH_84_MF_1081220191522_equal.csv"
    equal = pd.read_csv(di + filename, delimiter="\t")
    equal["Regime"] = "Medium"

    filename = "NS-ATracerH_84_MF_1121220192228_fast_0.00002.csv"
    fast = pd.read_csv(di + filename, delimiter="\t")
    fast["Regime"] = "Fast"

    filename = "NS-ATracerH_84_MF_1081220191539_slowsk.csv"
    slow = pd.read_csv(di + filename, delimiter=",")
    slow["Regime"] = "Slow"

    breakthrough = pd.concat([fast, slow, equal], axis=0, ignore_index=True)
    breakthrough["Heterogeneity"] = breakthrough["Variance"]
    breakthrough["VA"] = breakthrough["Heterogeneity"] * breakthrough["Anisotropy"]

    breakthrough["%ofhomogeneous"] = breakthrough["del"] * 100
    #    bth = breakthrough.rename(columns={'Scenario':'Trial'}).astype(str)

    l = []
    for i in range(len(breakthrough)):
        if breakthrough["Variance"][i] == 0.1:
            l.append(
                str(breakthrough["Variance"][i])
                + ":"
                + str(breakthrough["Anisotropy"][i])
            )
        else:
            l.append(
                str(int(breakthrough["Variance"][i]))
                + ":"
                + str(breakthrough["Anisotropy"][i])
            )

    breakthrough["Xlabels"] = l
    breakthrough = breakthrough.sort_values(by=["Variance", "Anisotropy"])

    combined_tracer = pd.read_csv(
        "X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv",
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

    return breakthrough, combined_tracer


def localmaxmin(y1, y2, init):
    ymax = np.max(np.max(np.append(y1, y2), axis=0))
    ymin = np.min(np.min(np.append(y1, y2), axis=0))
    xposmax = (np.argmax(y1) + init) * 5
    xposmin = (np.argmin(y1) + init) * 5
    print(ymax, ymin, xposmax, xposmin)
    return ymax, ymin, xposmax, xposmin
