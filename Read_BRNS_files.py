# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:25:44 2020

@author: khurana
"""
import numpy as np
import pandas as pd
import data_reader.data_processing as proc

#Required functions:
def calcvolmeanrate (ratedata, domain):
    firstnode = ratedata.iloc[-100]
    lastnode = ratedata.iloc[-1]        
    sidelength = domain.iloc[-1] - domain.iloc[-2]
    domainvol = domain.iloc[-1] - domain.iloc[-100]
    volmeanrate = ((firstnode + lastnode)*sidelength/2 + sum(ratedata.iloc[-101:-1])*sidelength)/domainvol
    
    return volmeanrate

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

gratenames = ["Aerobic respiration", "Ammonium respiration", "Nitrate respiration"]
respindx = range(66)[1:]

#list(ratenames.index(i) for i in rates)
#ratefiles = = list("xrate" + str(i) for i in respindx)
#r1 = np.loadtxt()

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

fileindex = list(t for t,values in scdict.items())
restime = list(values['time'] for t,values in scdict.items())
Pelist = list(values['Pe'] for t,values in scdict.items())

row = []
for f,t, p in zip(fileindex, restime, Pelist):
    path_concdata = "X:\Pap1_discussion\BRNS_PeDa_"+ str(f) + "\conc.dat"
    concdata = np.loadtxt(path_concdata, skiprows = 1)
    np.shape(concdata)
    delta = (concdata[0,1:] - concdata[-1,1:])
    reldelta = (concdata[0,1:] - concdata[-1,1:])/concdata[0,1:]
    ratedata = []
    for i in respindx:
        data = np.loadtxt("X:/Pap1_discussion/BRNS_PeDa_" + str(f) + "/" + "xrate" + str(i) + ".dat")
        ratedata.append(data[:,0])
    ratedf = pd.DataFrame.from_records(ratedata).T
    ratedf["X"] = data[:,1]
    np.shape(ratedf)
    r1 = calcvolmeanrate(ratedf[0], ratedf["X"])
    r2 = calcvolmeanrate(ratedf[1], ratedf["X"])
    r3 = calcvolmeanrate(ratedf[2], ratedf["X"])
    r4 = calcvolmeanrate(ratedf[3], ratedf["X"])
    r17 = calcvolmeanrate(ratedf[16], ratedf["X"])
    r18 = calcvolmeanrate(ratedf[17], ratedf["X"])
    r19 = calcvolmeanrate(ratedf[18], ratedf["X"])
    r20 = calcvolmeanrate(ratedf[19], ratedf["X"])
    r35 = calcvolmeanrate(ratedf[34], ratedf["X"])
    r36 = calcvolmeanrate(ratedf[35], ratedf["X"])
    r49 = calcvolmeanrate(ratedf[48], ratedf["X"])
    r50 = calcvolmeanrate(ratedf[49], ratedf["X"])
    r51 = calcvolmeanrate(ratedf[50], ratedf["X"])
    r52 = calcvolmeanrate(ratedf[51], ratedf["X"])
    r65 = calcvolmeanrate(ratedf[64], ratedf["X"])
    meanrates= [-1*(r1 + r2) - 1*(r49 + r50), 0.1*r65 - 0.1*(r3 + r4 + r19 + r20 + r35 + r36 + r51 + r52) - 0.5*(r49 + r50), -0.8*(r17 + r18)]
    meanrateswocoef= [-(r1 + r2) - (r49 + r50), r65 - (r3 + r4 + r19 + r20 + r35 + r36 + r51 + r52) - (r49 + r50), -(r17 + r18)]
    v_da = [abs(t*m) for m in meanrates]
    v_dawocoef = [abs(t*m) for m in meanrateswocoef]
    for c, gvar, rvar in zip(gvarindex, gvarnames, gratenames):
        mfrate = delta[c-1]/(0.2*t) #From Meile and Jung, 2010
        row.append([f, t, gvar, delta[c-1], reldelta[c-1], rvar,  abs(meanrates[gvarnames.index(gvar)]), v_da[gvarnames.index(gvar)], abs(meanrateswocoef[gvarnames.index(gvar)]), v_dawocoef[gvarnames.index(gvar)], mfrate, p])
    
df = pd.DataFrame.from_records(row, columns = ["Scenario", "Residence_time", "Chem", "delmassflux", "reldelmassflux", "Rate_type", "Volumetric_MeanRate", "Da", "WOCF_VMR", "WOCF_Da", "Massflux_MeanRate","Pe"])

df.to_csv("X:\Pap1_discussion\PE_Da_chem_summary.csv", sep = '\t')

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

alldata = pd.read_csv("X:\Pap1_discussion\PE_Da_chem_summary.csv", sep = '\t')
alldata["%reldelmassflux"] = alldata["reldelmassflux"]*100
subset = alldata[alldata['Scenario'].isin([1,2,3,4,5,7,8,22,23,24,25, 26, 27, 28, 29])]
subset = alldata[alldata['Scenario'].isin([1,2,3])]

plt.figure()
for c in ["Ammonium", "DO", "Nitrate"]:
    for p in [2,11,22, 50]:
        data = subset[(subset['Pe']==p) & (subset['Chem']==c)]
        plt.scatter(abs(data['Da']*data['Pe']), abs(data['Massflux_MeanRate']/data['Volumetric_MeanRate']), color = colorlist[[2,11,22,50].index(p)], marker = marklist[["Ammonium", "DO", "Nitrate"].index(c)], alpha = 0.5, label = data["Scenario"])
plt.xscale('log')
plt.yscale('log')
#plt.xlim ([0.0001, 10000])
#plt.ylim ([0.0001, 10000])
plt.xlabel ("Da/Pe")
plt.ylabel ("Concentration derived rate/Volumetric mean rate")
plt.legend(title = "Chem")
plt.legend(handles = patchlist, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=len(patchlist), mode="expand", borderaxespad=0., title = "Pe and Chem")
plt.savefig("Z:\Saturated_flow\diffusion_transient\Rates_DaPe_ratio.png", dpi = 300, pad_inches = 0.1)

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

#Figure S6 Comparison of Da numbers - mass flux derived numbers and volumetrically averaged rates
#Explore how mass flux derived rate can differ from volunetrically derived rates for a discussion on Da in both spatial heterogeneity and temporal heterogeneity
alldata = pd.read_csv("X:\Pap1_discussion\PE_Da_chem_summary.csv", sep = '\t')
alldata["%reldelmassflux"] = alldata["reldelmassflux"]*100
alldata["Ratio_MF_VRate"] = alldata['Massflux_MeanRate']/alldata['Volumetric_MeanRate']

data = alldata[alldata['Scenario'].isin([1,2,3,4,5,7,8,22,23,24,25, 26, 27, 28, 29])]

dummy = sns.boxplot( x=data["Chem"], y=data["Ratio_MF_VRate"], hue = data["Pe"], palette=my_pal)

medians = data.groupby(['Chem', 'Pe'])["Ratio_MF_VRate"].median().values
nobs = data['Chem'].value_counts().values
nobs = [str(x) for x in nobs.tolist()]
nobs = ["n: " + i for i in nobs]
 
# Add it to the plot
pos = range(len(nobs))
for tick,label in zip(pos,dummy.get_xticklabels()):
    dummy.text(pos[tick], medians[tick] + 0.03, nobs[tick], horizontalalignment='center', size=12, color='k', weight='semibold')

plt.ylabel("Ratio", fontsize = 15)
plt.xlabel ("Reactive species", fontsize = 15)

plt.title("Ratio of concentration derived rate and\n volumetrically averaged rate", fontsize = 15)
pic = dummy.get_figure()
pic.savefig("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/FigureS6_Ratio_MF_Vol_Rate.png",pi = 300, bbox_inches = 'tight', pad_inches = 0.1)