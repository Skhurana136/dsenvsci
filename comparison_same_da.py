# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:02:57 2020

@author: khurana
"""

#Comparison between different domains
import csv
import analyses.saturated_steady_state as sssa
import data_reader.data_processing as proc
import numpy as np

# Saturated flow regime
Regimes = ["Equal", "Fast"]
domains = ["", "Double_", "Half_"]
domains = ["", "Half_"]
youtlist = [50, 100, 25]
fpre = "NS-A"
fsuf = r"/"
gw = 1

original = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Default:
yin = 0
xleft = 0
xright = 30
# Assign index to Variable
doc1 = 10 - gw
dox1 = 11 - gw
Amm1 = 12 - gw
nitra1 = 17 - gw
sulpha1 = 22 - gw
tr1 = 29 - gw
Bfo1 = 8 - gw
Bfn1 = 15 - gw
Bfs1 = 20 - gw
Bfa1 = 25 - gw
Bmo1 = 9 - gw
Bmn1 = 16 - gw
Bms1 = 21 - gw
Bma1 = 26 - gw
Bifo1 = 13 - gw
Bifn1 = 18 - gw
Bifs1 = 23 - gw
Bifa1 = 27 - gw
Bimo1 = 14 - gw
Bimn1 = 19 - gw
Bims1 = 24 - gw
Bima1 = 28 - gw
POM1 = 30 - gw
vely = 5
velx = 4
vars = [doc1, dox1, Amm1, nitra1, sulpha1, Bmo1, Bma1, Bmn1, Bimo1, Bima1, Bimn1]
gvarnames = [
    "DOC",
    "DO",
    "Ammonium",
    "Nitrate",
    "Sulphate",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
    "Nitrogen",
    "TOC",
]
AFbiomassvars = [
    Bfo1,
    Bfa1,
    Bfn1,
    Bmo1,
    Bma1,
    Bmn1,
    Bifo1,
    Bifa1,
    Bifn1,
    Bimo1,
    Bima1,
    Bimn1,
]
AFbiomassgvarnames = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive fixed Aerobes",
    "Inactive fixed Ammonia oxidizers",
    "Inactive fixed Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
]

fnamemf = ("Z:/Saturated_flow/diffusion_transient/comparison_sameDa.csv")
csvfile = open(fnamemf, "w")
writer = csv.writer(
                csvfile,
                delimiter="\t",
                quotechar="\t",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
                )
writer.writerow(
                [
                        "Sno",
                        "Trial",
                        "Variance",
                        "Anisotropy",
                        "Chem",
                        "Inlet_total_mass_flux",
                        "Outlet_mass_flux",
                        "delmassflux",
                        "del2massflux",
                        "fdelmassflux",
                        "domain",
                        "Regime"
                        ]
                )
for Reg in Regimes:
    for domain, yout in zip(domains, youtlist):
        scdict = original
        Trial = ['H']
        Het = [0]
        Anis = [1]
        trialist = ['37', '38', '39', '40', '41', '42', '43', '44', '45']
        Trial.extend(list(int(t) for t,values in scdict.items() if t in trialist))
        Het.extend(list(values['Het'] for t,values in scdict.items() if t in trialist))
        Anis.extend(list(values['Anis'] for t,values in scdict.items() if t in trialist))
        directory = "Z:/Saturated_flow/diffusion_transient" +  "/" + domain + Reg + "AR_0/"
        mf = sssa.calcmassflux(Trial,Het, Anis,gw, directory, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
        biomasssum = sssa.calcsum(
                Trial,
                Het,
                Anis,
                gw,
                directory,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                AFbiomassvars,
                AFbiomassgvarnames,
                )
        # Writing the results into csv files
        if domain == "":
            domain = "Base"
        if Reg == "Equal":
            r = "Medium"
        else:
            r = Reg            
        for j in Trial:
            for i in gvarnames:
                idx = Trial.index(j) * len(gvarnames) + gvarnames.index(i)
                writer.writerow([idx + 1,j,Het[Trial.index(j)],Anis[Trial.index(j)],i,mf[idx, 4],mf[idx, 5],mf[idx, 6],mf[idx, 7],mf[idx, 8], domain, Reg])
            df = np.load(directory + fpre + str(j) + fsuf + fpre + str(j) + "_df.npy")
            vel = np.mean(df[vely - 3, -1, :, :])
            idx = idx + 1
            writer.writerow([idx + 1,j,Het[Trial.index(j)],Anis[Trial.index(j)],"Velocity",vel,vel,vel,vel,vel, domain, Reg])
csvfile.close()

import seaborn as sns
import matplotlib.pyplot as plt
dom = "Half"
vely = 5

for Regime in ["Equal", "Fast"]:
    directory = "X:/Saturated_flow/changedkindox_transient/" + dom + "_" + Regime + "AR_0/" + fpre #change directory as per flow regime
    for j in Trial:
        print(j)
        df = np.load(directory + j + fsuf + fpre + j + "_df.npy")
        print(np.mean(df[2, -1, :, :]))
        plt.figure()
        sns.heatmap(df[2,-1,:,:])
        print (Regime, dom, j, np.mean(df[vely - 3, -1, :, :]))
            
            
            
            
            