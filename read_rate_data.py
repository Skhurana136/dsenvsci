# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:53:39 2020

@author: khurana
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:24:47 2020

@author: khurana
"""
import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import data_reader.reader as rdr

# Saturated flow regime
Regimes = ["Slow", "Equal", "Fast"]
fpre = "NS-A"
fsuf = r"/"
filename = "ratesAtFinish.dat"

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
ratenames = proc.masterrates("saturated")

#Domain description
xleft = 0
xright = 30
yin = 0
yout = 50
vedge = 0.005
velem = 0.01
domainvol = 4*vedge**2 + 2*vedge*velem*(29+49) + (1581 - 2*49 - 2*29)*velem**2

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

rates = [
    "Fixedaeroresp",
    "Mobaeroresp",
    "Fixedaerogwth",
    "Mobaerogwth",
    "Fixednitraresp",
    "Mobnitraresp",
    "Fixednitragwth",
    "Mobnitragwth",
    "Fixedammresp",
    "Mobammresp",
    "Fixedammgwth",
    "Mobammgwth",
    "Hydrolysis"
    ]

respindx = np.array(list(ratenames.index(i) for i in rates))

gratenames = [
    "Immobile aerobic respiration",
    "Mobile aerobic respiration",
    "Immobile aerobic growth",
    "Mobile aerobic growth",
    "Immobile nitrate respiration",
    "Mobile nitrate respiration",
    "Immobile nitrate growth",
    "Mobile nitrate growth",
    "Immobile ammonium respration",
    "Mobile ammonium respiration",
    "Immobile ammonium growth",
    "Mobile ammonium growth",
    "Hydrolysis"
    ]
    
#Reading rates of interest, saving in numpy array, calculating volume averaged rates
row = []
for Reg in Regimes:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    if (Reg == "Equal"):
        Reg = "Medium"
    for t in Trial:
        print(t)
        filepath = directory + fpre + str(t) + fsuf + filename
        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=16 + respindx)
        ratesarray = rdr.Converttoarray_1581(M, "rates")
        np.save(directory + fpre + str(t) + fsuf + fpre + str(t) + "_resp_growth_hydrolysis_rates.npy", ratesarray)
        for i in range(len(respindx)):
            sumratecorner = (ratesarray[i,yin,xleft] + ratesarray[i,yout,xleft] + ratesarray[i,yin,xright] + ratesarray[i,yout,xright])*vedge**2
            sumratebound = (sum(ratesarray[i,yin,xleft+1:xright]) + sum(ratesarray[i,yout,xleft+1:xright]) + sum(ratesarray[i,yin+1:yout,xleft]) + sum(ratesarray[i,yin+1:yout,xright]))*vedge*velem
            sumrateelem = sum(ratesarray[i,yin+1:yout, xleft+1:xright])*velem**2
            meanrate = sum(sumratecorner+sumratebound+sumrateelem)/domainvol
            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate])

df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Rate_type', 'Meanrate'])

df.to_csv(r"Z:\Saturated_flow\diffusion_transient\Mean_rate_ss.csv", sep = '\t')