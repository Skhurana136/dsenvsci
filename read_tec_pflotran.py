# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:42:48 2020

@author: khurana
"""
#Reading tec files and storing them as numpy array - possible options to read 51x31 matrices or 101x61 matrices in data_reader library
import numpy as np

#set up basic constants 
fsuf = r"/"
filename = "pflotran_Swamini_2909-001.tec" #same filename in each subfolder
horiznodes = 31

# Default:
directory = "Z:/pflotran_het"
fwithd = directory + fsuf + filename #complete path to file
print("Reading tech file....")

fid = open(fwithd, "r")
tline_old = fid.readline()
n = 1
while len(tline_old) > 0:
    n = n + 1
    tline_old = fid.readline()

fid.seek(0)
line_ex = (fid.readline())  # check if readline includes or excludes new line character. chose function that removes new line character
Timepoint = line_ex.split(" ")[4]
line_ex = fid.readline()
Headers = line_ex[11:].split(",")
fid.close()
with open(fwithd) as fid:
    M = np.loadtxt(fid,dtype=str,delimiter="  ",skiprows=3)

D = M[:,:-2]

vertnodes = 51
horiznodes = 31

Headers = np.shape(D)[1]
size = np.shape(D)[0]
df2 = np.ndarray([Headers - 3, vertnodes, horiznodes])
for i in range(Headers - 3):
    counter = 0
    for k in range(vertnodes):
        for l in range(horiznodes):
            df2[i, k, l] = D[counter, i+3]
            counter = counter + 1

df = np.ndarray([Headers - 3, vertnodes, horiznodes])
for j in range(vertnodes):
    df[:, vertnodes - 1 - j, :] = df2[:, j, :]

import seaborn as sns

sns.heatmap(df[-2, :, :])