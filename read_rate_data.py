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

def calcvolmeanrate (ratedata, domain):
    firstnode = ratedata.iloc[-100]
    lastnode = ratedata.iloc[-1]        
    sidelength = domain.iloc[-1] - domain.iloc[-2]
    domainvol = domain.iloc[-1] - domain.iloc[-100]
    volmeanrate = ((firstnode + lastnode)*sidelength/2 + sum(ratedata.iloc[-101:-1])*sidelength)/domainvol
    
    return volmeanrate

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

ratesofinterest = [1,2,3,4,17,18,19,20,35,36,49,50,51,52,65]

respindx = np.array(range(66)[1:]) #List all rates

gvarnames = ["DO", "Ammonium", "Nitrate"]
gratenames = ["DO removal", "Ammonium removal", "Nitrate removal"]

#Reading rates of interest and saving in numpy array

for Reg in Regimes:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    for t in Trial:
        print(t)
        filepath = directory + fpre + str(t) + fsuf + filename
        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=18 + respindx)
        ratesarray = rdr.Converttoarray_1581(M, "rates")
        np.save(directory + fpre + str(t) + fsuf + fpre + str(t) + "_all_rates.npy", ratesarray)
  
#Calculating mean rates
row = []        
for Reg in Regimes:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    if (Reg == "Equal"):
        Reg = "Medium"
    for t in Trial:
        print(t)
        ratesarray = np.load(directory + fpre + str(t) + fsuf + fpre + str(t) + "_all_rates.npy")
        r = []
        mr = []
        for i in list(k-1 for k in respindx):
            sumratecorner = (ratesarray[i,yin,xleft] + ratesarray[i,yout,xleft] + ratesarray[i,yin,xright] + ratesarray[i,yout,xright])*vedge**2
            sumratebound = (sum(ratesarray[i,yin,xleft+1:xright]) + sum(ratesarray[i,yout,xleft+1:xright]) + sum(ratesarray[i,yin+1:yout,xleft]) + sum(ratesarray[i,yin+1:yout,xright]))*vedge*velem
            sumrateelem = sum(ratesarray[i,yin+1:yout, xleft+1:xright])*velem**2
            volmeanrate = sum(sumratecorner+sumratebound+sumrateelem)/domainvol
            meanrate = np.mean(ratesarray[i,:,:])
            r.append(volmeanrate)
            mr.append(meanrate)
        volmeanrates= [-1*(r[0] + r[1]) - 1*(r[48] + r[49]), 0.1*r[64] - 0.1*(r[2] + r[3] + r[18] + r[19] + r[34] + r[35] + r[50] + r[51]) - 0.5*(r[48] + r[49]), -0.8*(r[16] + r[17])]        
        meanrates= [-1*(mr[0] + mr[1]) - 1*(mr[48] + mr[49]), 0.1*mr[64] - 0.1*(mr[2] + mr[3] + mr[18] + mr[19] + mr[34] + mr[35] + mr[50] + mr[51]) - 0.5*(mr[48] + mr[49]), -0.8*(mr[16] + mr[17])]                
        for gvar, rvar in zip(gvarnames, gratenames):
            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gvar, rvar, abs(volmeanrates[gvarnames.index(gvar)]), abs(meanrates[gvarnames.index(gvar)])])

#        for i in range(len(respindx)):
#            sumratecorner = (ratesarray[i,yin,xleft] + ratesarray[i,yout,xleft] + ratesarray[i,yin,xright] + ratesarray[i,yout,xright])*vedge**2
#            sumratebound = (sum(ratesarray[i,yin,xleft+1:xright]) + sum(ratesarray[i,yout,xleft+1:xright]) + sum(ratesarray[i,yin+1:yout,xleft]) + sum(ratesarray[i,yin+1:yout,xright]))*vedge*velem
#            sumrateelem = sum(ratesarray[i,yin+1:yout, xleft+1:xright])*velem**2
#            volmeanrate = sum(sumratecorner+sumratebound+sumrateelem)/domainvol
#            meanrate = np.mean(ratesarray[i,:,:])
#            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], volmeanrate, meanrate])

df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Chem', 'Rate_type', 'Volumetric_Meanrate', 'Arithmetic_Meanrate'])

df.to_csv(r"Z:\Saturated_flow\diffusion_transient\Complex_Mean_rate_ss.csv", sep = '\t')