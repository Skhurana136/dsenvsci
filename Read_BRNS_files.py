# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:25:44 2020

@author: khurana
"""
import numpy as np
import pandas as pd
import data_reader.data_processing as proc

#variables of interest:

chemdict = proc.masterdissolvedspecies() #master dictionary of dissolved species
microdict = proc.mastermicrobialspecies() #master dictionary of microbial species
ratenames = proc.masterrates("saturated")

#look for file names
Chems = list(values['Var'] for t,values in chemdict.items())
Chemnames = list(values['Graphname'] for t,values in chemdict.items())
chemindx = list(values['TecIndex']-7 for t,values in chemdict.items())
microbes = list(values['Var'] for t,values in microdict.items())
microbenames = list(values['Graphname'] for t,values in microdict.items())
microbeindx = list(values['TecIndex']-5 for t,values in microdict.items())

gvarnames = ["DO", "Ammonium", "Nitrate"]
indices = list(Chemnames.index(i) for i in gvarnames)
gvarindex = list(chemindx[i] for i in indices)

rates = [
    "Fixedaeroresp",
    "Fixedammresp",
    "Fixednitraresp"
    ]
gratenames = ["Immobile aerobic respiration", "Immobile ammonium respiration", "Immobile nitrate respiration"]
respindx = list(ratenames.index(i)+1 for i in rates)
ratefiles = list("xrate"+str(i) for i in respindx)

scdict = {
        1: {'Pe': 2, 'time': 1316},
        2: {'Pe': 11, 'time': 132},
        3: {'Pe': 22, 'time': 13},
        4: {'Pe': 50, 'time': 13},
        5: {'Pe': 22, 'time': 132},
        7: {'Pe': 2, 'time': 132},
        8: {'Pe': 11, 'time': 13},
        9: {'Pe': 2, 'time': 1.3},
        10: {'Pe': 11, 'time': 6.6},
        12: {'Pe': 50, 'time': 26.3},
        13: {'Pe': 450, 'time': 263},
        15: {'Pe': 2, 'time': 26},
        16: {'Pe': 22, 'time': 263},
        17: {'Pe': 50, 'time': 526},
        19: {'Pe': 11, 'time': 7368},
        20: {'Pe': 22, 'time': 14474},
        22: {'Pe': 2, 'time': 13},
        23: {'Pe': 11, 'time': 1316},
        24: {'Pe': 22, 'time': 1316},
        25: {'Pe': 50, 'time': 1316},
        26: {'Pe': 2, 'time': 13},
        27: {'Pe': 50, 'time': 132},
        28: {'Pe': 2, 'time': 132},
        29: {'Pe': 22, 'time': 132},
        }

scdict = {
        1: {'Pe': 2, 'time': 1316},
        2: {'Pe': 12, 'time': 132},
        3: {'Pe': 22, 'time': 13},
        4: {'Pe': 47, 'time': 13},
        5: {'Pe': 32, 'time': 132},
        7: {'Pe': 4, 'time': 132},
        8: {'Pe': 10, 'time': 13},
        9: {'Pe': 2, 'time': 1.3},
        10: {'Pe': 11, 'time': 6.6},
        12: {'Pe': 50, 'time': 26.3},
        13: {'Pe': 450, 'time': 263},
        15: {'Pe': 2, 'time': 26},
        16: {'Pe': 22, 'time': 263},
        17: {'Pe': 50, 'time': 526},
        19: {'Pe': 11, 'time': 7368},
        20: {'Pe': 22, 'time': 14474},
        22: {'Pe': 3, 'time': 13},
        23: {'Pe': 11, 'time': 1316},
        24: {'Pe': 22, 'time': 1316},
        25: {'Pe': 46, 'time': 1316},
        }

fileindex = list(t for t,values in scdict.items())
restime = list(values['time'] for t,values in scdict.items())
Pelist = list(values['Pe'] for t,values in scdict.items())

row = []
#Path to file

for f,t, p in zip(fileindex, restime, Pelist):
    path_concdata = "X:\Pap1_discussion\BRNS_PeDa_"+ str(f) + "\conc.dat"
    concdata = np.loadtxt(path_concdata, skiprows = 1)
    np.shape(concdata)
    delta = (concdata[0,1:] - concdata[-1,1:])/concdata[0,1:]
    for c, gvar, rfile, r, rvar in zip(gvarindex, gvarnames, ratefiles, rates, gratenames):
        concdata = np.loadtxt(path_concdata, skiprows = 1)    
        path_ratedata = "X:\Pap1_discussion\BRNS_PeDa_" + str(f) +  r"/" + rfile + ".dat"
        ratedata = np.loadtxt(path_ratedata)
        firstnode = ratedata[-100,0]
        lastnode = ratedata[-1,0]
        sidelength = ratedata[-1,1] - ratedata[-2,1]
        domainvol = ratedata[-1,1] - ratedata[-100,1]
        volmeanrate = ((firstnode + lastnode)*sidelength/2 + sum(ratedata[-101:-1,0])*sidelength)/domainvol
        V_Da = t*volmeanrate
        row.append([f, t, gvar, delta[c-1], rvar,  volmeanrate, V_Da, p])

df = pd.DataFrame.from_records(row, columns = ["Scenario", "Residence_time", "Chem", "reldelmassflux", "Rate_type", "Volumetric_MeanRate", "Da", "Pe"])

df.to_csv("X:\Pap1_discussion\PE_Da_chem_summary.csv", sep = '\t')
subset = df.sort_values(by=['Pe','Da'])

import matplotlib.pyplot as plt
marklist = ['o', '^', 's']
colorlist = ['indianred', 'g', 'steelblue', 'orange']
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
red_patch = mpatches.Patch(color="indianred", alpha = 0.5, label="2")
green_patch = mpatches.Patch(color="g", alpha = 0.5, label="11")
blue_patch = mpatches.Patch(color="steelblue", alpha = 0.5, label="22")
orange_patch = mpatches.Patch(color="orange", alpha = 0.5, label="50")
do = mlines.Line2D([], [], Linestyle = 'None', marker = 'o', alpha = 0.5, label="Ammonium")
amm = mlines.Line2D([], [], Linestyle = 'None', marker = '^', alpha = 0.5, label="DO")
nitra = mlines.Line2D([], [], Linestyle = 'None', marker = 's', alpha = 0.5, label="Nitrate")
patchlist = [red_patch, green_patch, blue_patch, orange_patch, do, amm, nitra]

subset = pd.read_csv("X:\Pap1_discussion\PE_Da_chem_summary.csv", sep = '\t')
subset["%reldelmassflux"] = subset["reldelmassflux"]*100
subset = subset[subset['Pe']<100]
subset = subset[subset['Da']<300]

subset = subset[subset['Scenario'].isin([1,2,3,4,5,7,8,22,23,24,25, 26, 27, 28, 29])]
subset = subset[subset['Scenario'].isin([1,2,3])]

plt.figure()
for c in ["Ammonium", "DO", "Nitrate"]:
    for p in [2,11,22, 50]:
        data = subset[(subset['Pe']==p) & (subset['Chem']==c)]
        plt.scatter(data['Da'], data['%reldelmassflux'], color = colorlist[[2,11,22,50].index(p)], marker = marklist[["Ammonium", "DO", "Nitrate"].index(c)], alpha = 0.5, label = data["Scenario"])
plt.xscale('log')
#plt.yscale('log')
plt.xlim(left = 0.00001)
#plt.ylim(bottom = 0.1)
plt.xlabel ("Dammkohler number")
plt.ylabel ("Normalized removal")
plt.legend(title = "Chem")
plt.legend(handles = patchlist, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=len(patchlist), mode="expand", borderaxespad=0., title = "Pe and Chem")
plt.savefig("Z:\Saturated_flow\diffusion_transient\PE_sameDa_Chem_Removal.png", dpi = 300, pad_inches = 0.1)

plt.figure()
for f in fileindex[:3]:
    subset2 = subset[subset['Scenario']==f]
    for p in subset2['Pe'].unique().tolist():
        data = subset2[(subset2['Pe']==p) ]
        plt.scatter(data['Da'], data['reldelmassflux'], marker = marklist[subset['Pe'].unique().tolist().index(p)], alpha = 0.5, label = p)
plt.xscale('log')
plt.xlim(left = 0.00001)
plt.xlabel ("Dammkohler number")
plt.ylabel ("Normalized removal")
plt.legend(title = "Pe")