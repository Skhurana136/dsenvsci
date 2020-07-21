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

fileindex = [1, 2, 3, 9, 10, 12, 13]
Pelist = [2, 11, 22, 2, 11, 45, 450]
restime = [1316, 132, 13, 1.3, 6.6, 13.2, 26.3, 263]
row = []
#Path to file

for f,t, p in zip(fileindex,restime, Pelist):
    path_concdata = "X:\Pap1_discussion\BRNS_PeDa_"+ str(f) + "\conc.dat"
    concdata = np.loadtxt(path_concdata, skiprows = 1)
    np.shape(concdata)
    delta = (concdata[-1,1:])/concdata[0,1:]
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
        row.append([f, t, gvar, delta[c-1], rvar, V_Da, p])

df = pd.DataFrame.from_records(row, columns = ["Scenario", "Residence_time", "Chem", "reldelmassflux", "Rate_type", "Da", "Pe"])