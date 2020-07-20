# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:27:25 2020

@author: khurana
"""

import numpy as np
import data_reader.data_processing as proc

# Saturated flow regime
Regimes = ["Slow", "Equal", "Fast"]
fpre = "NS-A"
fsuf = r"/"
filename = "ratesAtFinish.dat"

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
ratenames = proc.masterrates("saturated")

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

rates = [
    "Fixedaeroresp",
    "Mobaeroresp",
    "Fixednitraresp",
    "Mobnitraresp",
    "Fixedammresp",
    "Mobammresp",
    ]

respindx = np.array(list(ratenames.index(i) for i in rates))
gvarnames = ["DO", "DO", "Nitrate", "Nitrate", "Ammonium", "Ammonium"]

#Calculating volume averaged rates

for Reg in Regimes:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    for t in Trial:
        print(t)
        filepath = directory + fpre + str(t) + fsuf + fpre + str(t) + "_rates.npy"
        ratesarray = np.load(filepath)
        for i in range(len(respindx)):
            np.save(directory + fpre + str(t) + fsuf + fpre + str(t) + "_rates.npy", ratesarray)
            
            