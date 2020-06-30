# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:42:48 2020

@author: khurana
"""
#Reading tec files and storing them as numpy array - possible options to read 51x31 matrices or 101x61 matrices in data_reader library
import numpy as np
import data_reader.reader as rdr

#set up basic constants 
fpre = "RF-A"
masterTrial, masterAnis, masterHet = proc.masterScenarios()
masterTrial = [
    "H",
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
]
masterHet = [
    0,
    0.1,
    0.1,
    0.1,
    1,
    1,
    1,
    10,
    10,
    10,
    0.1,
    0.1,
    0.1,
    1,
    1,
    1,
    10,
    10,
    10,
    0.1,
    0.1,
    0.1,
    1,
    1,
    1,
    10,
    10,
    10,
    0.1,
    0.1,
    0.1,
    1,
    1,
    1,
    5,
    5,
    5,
    10,
    10,
    10,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
]
masterAnis = [
    1,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
    2,
    5,
    10,
]
fsuf = r"/"

filename = "model_domain_quad.tec" #same filename in each subfolder

# setup what we really want to read
# Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

# Reading and storing in numpy array
for Regime in ["Equal", "Fast"]:
    directory = r"X:/Richards_flow/Tracer_studies/" + Regime + "AR/" + fpre #change directory as per flow regime
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