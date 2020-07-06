# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:42:48 2020

@author: khurana
"""
#Reading tec files and storing them as numpy array - possible options to read 51x31 matrices or 101x61 matrices in data_reader library
import numpy as np
import data_reader.reader as rdr
import data_reader.data_processing as proc

#set up basic constants 
fpre = "RF-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
fsuf = r"/"
filename = "model_domain_quad.tec" #same filename in each subfolder

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

#Variations:
#notlist = [43,54]
#Trial = list(t for t,values in scdict.items() if t not in notlist)
#Het = list(values['Het'] for t,values in scdict.items() if t not in notlist)
#Anis = list(values['Anis'] for t,values in scdict.items() if t not in notlist)

# Reading and storing in numpy array
for Regime in ["Slow"]:
    directory = r"X:/Richards_flow/Tracer_studies/" + Regime + "AR_0/" + fpre #change directory as per flow regime
    for j in range(len(Trial)):
        print(str(Trial[j]))
        fwithd = directory + str(Trial[j]) + fsuf + filename #complete path to file
        print("Reading tech file....")
        size, steps, Headers, D = rdr.readTecfile(fwithd) #read the tec file
        print("Converting to array....")
        df = rdr.Converttomarr_1581(D) #Convert the coarse grid in to 51x31 array
        print("Saving numpy array....") 
        np.save(directory + str(Trial[j]) + "_df", df)
        # Test for correct orientation of the data
        for i in range(np.shape(df)[0]):
            print(Headers[i + 3], np.mean(df[i, steps - 1, 0, :]))