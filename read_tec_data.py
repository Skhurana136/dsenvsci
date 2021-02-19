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
Regimes = ["Slow", "Equal", "Fast"]
domains = ["", "Half", "Double", "Big"]
domainodes = {"": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
fpre = "RF-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
fsuf = r"/"
filename = "model_domain_quad.tec" #same filename in each subfolder
horiznodes = 31

# Default:
trialist = ['H','37', '38', '39', '40', '41', '42', '43', '44', '45']
#trialist = ['37','45']
reginvest = Regimes [1:]
domaininvest = domains [1:2]
# Reading and storing in numpy array
for Regime in reginvest:
    for dom in domaininvest:
        directory = "X:/Saturated_flow/changedkindox_transient/" + dom + "_" + Regime + "AR_0/" + fpre #change directory as per flow regime
        #directory = "X:/Saturated_flow/Steady_state/Tracer_studies/" + dom + "_" + Regime + "AR/" + fpre #change directory as per flow regime
        for j in trialist:
            print(j)
            fwithd = directory + j + fsuf + filename #complete path to file
            print("Reading tech file....")
            size, steps, Headers, D = rdr.readTecfile(fwithd) #read the tec file
            print("Converting to array....")
            df = rdr.Convertetoarray(D, "tec", domainodes[dom]['ynodes'], horiznodes) #Convert the coarse grid in to 51x31 array
            print("Saving numpy array....") 
            np.save(directory + j + fsuf + fpre + j +  "_df", df)
            # Test for correct orientation of the data
            for i in range(np.shape(df)[0]):
                print(Headers[i + 3], np.mean(df[i, -1, 0, :]))
            print(np.mean(df[2, -1, :, :]))

for Regime in Regimes:
    directory = "X:/Richards_flow/Tracer_studies/" + Regime + "AR/" + fpre #change directory as per flow regime
    for j in list(scdict.keys()):
        print(j)
        fwithd = directory + j + fsuf + filename #complete path to file
        print("Reading tech file....")
        size, steps, Headers, D = rdr.readTecfile(fwithd) #read the tec file
        print("Converting to array....")
        df = rdr.Convertetoarray(D, "tec", 51, horiznodes) #Convert the coarse grid in to 51x31 array
        print("Saving numpy array....") 
        np.save(directory + j + fsuf + fpre + j +  "_df", df)
        # Test for correct orientation of the data
        for i in range(np.shape(df)[0]):
            print(Headers[i + 3], np.mean(df[i, -1, 0, :]))
            print(np.mean(df[2, -1, :, :]))


def file(j):
    fwithd = directory + j + fsuf + filename #complete path to file
    print("Reading tech file....")
    size, steps, Headers, D = rdr.readTecfile(fwithd) #read the tec file
    print("Converting to array....")
    df = rdr.Convertetoarray(D, "tec", domainodes[dom]['ynodes'], horiznodes) #Convert the coarse grid in to 51x31 array
    print("Saving numpy array....") 
    # Test for correct orientation of the data
    for i in range(np.shape(df)[0]):
        print(Headers[i + 3], np.mean(df[i, -1, 0, :]))
    print(np.mean(df[2, -1, :, :]))