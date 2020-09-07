# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 16:32:33 2020

@author: khurana
"""

#Comparison between different domains
import pandas as pd
import data_reader.data_processing as proc
import numpy as np
import analyses.saturated_transient as sta

# Saturated flow regime
domainodes = {"Base": {'yout' : 50},
              "Big" : {'yout' : 125},
              "Double" : {'yout' : 100},
              "Half" : {'yout' : 25}}
regsteps = {"Slow": {'step' : 0.01, 'freq' : 100}, 
            "Equal" : {'step' : 0.005, 'freq' : 2},
            "Fast" : {'step' : 0.0002, 'freq' : 5}}
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

horiznodes = 31
fpre = "NS-A"
fsuf = r"/"
gw = 1

# Default:
yin = 0
xleft = 0
xright = 30

#Investigating scenarios:
trialist = ['H','37', '38', '39', '40', '41', '42', '43', '44', '45']

trialist = ['37', '45']
reginvest = list(regsteps.keys())[:1]
domaininvest = list(domainodes.keys())[1:2]

#Exploring velocity field of tracer studies in different domains
import seaborn as sns
import matplotlib.pyplot as plt

for Reg in reginvest:
    step = regsteps[Reg]['freq']*regsteps[Reg]['step']
    for domain in domaininvest:
        directory = "X:/Saturated_flow/Steady_state/Tracer_studies/" + domain + "_" + Reg + "AR/" + fpre #change directory as per flow regime
        print (Reg, domain)
        yout = domainodes[domain]['yout']
        for j in trialist:
            print(j)
            df = np.load(directory + j + fsuf + fpre + j + "_df.npy")
            print(np.mean(df[2, -1, :, :]))
            plt.figure()
            sns.heatmap(df[2,-1,:,:])
            plt.title(j)
            print (Reg, domain, j, np.mean(df[2, -1, :, :]))

#Tracer studies for different domains:
tr1 = 7 - gw
vely = 5
velx = 4
vars = [tr1]
gvarnames = ["Tracer"]

for Reg in reginvest:
    step = regsteps[Reg]['freq']*regsteps[Reg]['step']
    for domain in domaininvest:
        directory = "X:/Saturated_flow/Steady_state/Tracer_studies/" + domain + "_" + Reg + "AR/" #change directory as per flow regime
        print (Reg, domain)
        yout = domainodes[domain]['yout']
        for j in trialist:
            df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                'H', scdict['H']['Het'], scdict['H']['Anis'], gw, directory, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames
                )
            print(j)
            plt.figure()
            plt.plot(conctime[-1,:,0])
            plt.title(j + ": Concentration profile at the end time step")
            plt.figure()
            plt.plot(conctime[:,-1,0])
            plt.title(j + ": Concentration profile at the outlet")
            
idx = 1
row = []
for Reg in reginvest:
    s = regsteps[Reg]['freq']*regsteps[Reg]['step']
    for domain in domaininvest:
        d = "X:/Saturated_flow/Steady_state/Tracer_studies" + "/" + domain+ "_" + Reg + "AR/"
        print (Reg, domain)
        yout = domainodes[domain]['yout']
        df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                'H', scdict['H']['Het'], scdict['H']['Anis'], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames
                )
        Time = np.where(np.round(conctime[:, -1, 0], 3) > 9)
        initial = s * Time[0][0]
        for j in trialist:
            df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                    j, scdict[j]['Het'], scdict[j]['Anis'], gw, d, fpre, fsuf, yin,yout, xleft, xright, vars, gvarnames
                    )
            print(conctime[-1,-1,0])
            Time = np.where(np.round(conctime[:, -1, 0], 3) > 9)
            Time2 = np.where(np.round(df[-1, :, -1, :], 3) > 9)
            print(s * Time[0][0], s * Time2[0][0], initial, (s * Time[0][0]) / initial)
            row.append([idx, j,scdict[j]['Het'], scdict[j]['Anis'],"Tracer",s * Time[0][0],(s * Time[0][0]) / initial,Reg,domain])
            idx = idx + 1

df = pd.DataFrame.from_records(row, columns = ["Sno", "Trial", "Variance", "Anisotropy", "Chem", "Time", "fraction", "Regime", "Domain"])
df.to_csv("X:/Saturated_flow/Steady_state/Tracer_studies/tracer_half_double.csv", sep = "\t")

#Explore change in mass flux in the different domains and compare with the original domain
vardict = proc.masterdissolvedspecies()
gvarnames = list(vardict.keys())
vars = list(vardict[t]["TecIndex"]-1 for t in gvarnames)
gvarnames.append("Nitrogen")
gvarnames.append("TOC")


row = []
for Reg in reginvest:
    for domain in domaininvest:
        yout = domainodes[domain]['yout']
        if domain != "":
            domadd = domain + "_"
        else:
            domadd = ""
        d = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_0/"#change directory as per flow regime
        print (Reg, domain)
        df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                'H', scdict['H']['Het'], scdict['H']['Anis'], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames
                )
        basecasein = conctime[-1, 0, :]
        basecaseout = conctime[-1, -1, :]
        basedelconc = basecasein - basecaseout
        basereldelconc = 100*basedelconc/basecasein
        for j in trialist:
            df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                    j, scdict[j]['Het'], scdict[j]['Anis'], gw, d, fpre, fsuf, yin,yout, xleft, xright, vars, gvarnames
                    )
            concin = conctime[-1, 0, :]
            concout = conctime[-1, -1, :]
            delconc = concin - concout
            reldelconc = 100*delconc/concin
            fractiondelconc = delconc/basedelconc
            fractionreldelconc = reldelconc/basereldelconc
            for g in gvarnames:
                row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, g, delconc[gvarnames.index(g)], reldelconc[gvarnames.index(g)], fractiondelconc[gvarnames.index(g)], fractionreldelconc[gvarnames.index(g)]])

df = pd.DataFrame.from_records(row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Chem", "delconc", "reldelconc", "fraction_delconc", "fraction_reldelconc"])
df.to_csv("X:/Saturated_flow/changedkindox_transient/" + domain + "_massflux.csv", sep = "\t")
