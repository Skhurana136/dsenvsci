# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:25:44 2020

@author: khurana
"""
import numpy as np
import data_reader.data_processing as proc

#variables of interest:

chemdict = proc.masterdissolvedspecies() #master dictionary of dissolved species
microdict = proc.mastermicrobialspecies() #master dictionary of microbial species
ratenames = proc.masterrates("saturated")

#look for file names
Chems = list(values['Var'] for t,values in chemdict.items())
chemindx = list(values['TecIndex']-5 for t,values in chemdict.items())
microbes = list(values['Var'] for t,values in microdict.items())
microbeindx = list(values['TecIndex']-5 for t,values in microdict.items())

rates = [
    "Fixedaeroresp",
    "Mobaeroresp",
    "Fixednitraresp",
    "Mobnitraresp",
    "Fixedammresp",
    "Mobammresp",
    ]

respindx = list(ratenames.index(i)+1 for i in ratenames)
ratefiles = list("xrate"+str(i) for i in respindx)

#Path to file
path_data = "X:\Pap1_discussion\BRNS_PeDa_12\conc.dat"
columnstoread = chemindx.append(microbeindx) #WIP not working yet

data = np.loadtxt(path_data, skiprows = 1)#, usecols = chemindx)
np.shape(data)

delta = (data[0,1:] - data[-1,1:])/data[0,1:]
delta[chemindx[1]]
