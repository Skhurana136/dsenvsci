# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 14:26:55 2020

@author: khurana
"""

#Reading tec files and storing them as numpy array - possible options to read 51x31 matrices or 101x61 matrices in data_reader library
import data_reader.reader as rdr
import data_reader.data_processing as proc
import numpy as np

#set up basic constants 
fpre = "NS-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
fsuf = r"/"
filename = "model_domain_quad.tec" #same filename in each subfolder

# Default:
trialist = ['H','37', '38', '39', '40', '41', '42', '43', '44', '45']
Trial = list(t for t,values in scdict.items() if t in trialist)
Het = list(values['Het'] for t,values in scdict.items() if t in Trial)
Anis = list(values['Anis'] for t,values in scdict.items() if t in Trial)

#Variations:
#notlist = [43,54]
#Trial = list(t for t,values in scdict.items() if t not in notlist)
#Het = list(values['Het'] for t,values in scdict.items() if t not in notlist)
#Anis = list(values['Anis'] for t,values in scdict.items() if t not in notlist)

#domdict = {"Half": 'verticalnodes' = 26, "Double": 'verticalnodes' = 101}
dom = "Big"
horiznodes = 31
vely = 5

# Reading and storing in numpy array
Regime = "Equal"
directory = "X:/Saturated_flow/Steady_state/Tracer_studies/" + dom + "_" + Regime + "AR/" #change directory as per flow regime
j = '45'
fwithd = directory + j + fsuf + filename #complete path to file

def file():
    print(fwithd)
    print("Reading tech file....")
    size, steps, Headers, D = rdr.readTecfile(fwithd) #read the tec file
    print("Converting to array....")
    df = rdr.Convertetoarray(D, "tec", 126, horiznodes) #Convert the coarse grid in to 51x31 array
    print(np.mean(df[vely - 3, -1, :, :]))
    print(np.mean(df[0, -1, 0, :]))