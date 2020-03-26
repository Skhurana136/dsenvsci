# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:41:04 2019

@author: khurana
"""

import math
import numpy  as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import pandas as pd
from matplotlib.colors import LogNorm
import linecache

#reading Tec files
#filename = 'model_domain_quad.tec'

def readTecfile (filename):
    numberOfHeaderLines = 3

    # pre-analysis to extract total number of lines in the data file which is 
    # used further on to calculate number of time steps
    fid   = open(filename, "r")

    #with open("t149_domain_quad.tec", "r") as fid
    tline_old = fid.readline()
    n = 1
    while len(tline_old)>0:
        n = n+1
        tline_old = fid.readline()

    #go to beginning of file;
    fid.seek(0)
    line_ex = fid.readline() #check if readline includes or excludes new line character. chose function that removes new line character
    line_ex_a = line_ex.lstrip('VARIABLES  = "X","Y","Z",') #deleting this section of the variable string
    #Alternatively, slice string:
    line_ex_final = line_ex[12:len(line_ex)].strip()
    Headers = line_ex_final.split(',') #unpacking the variables
#    Headers_new = [] REDUNDANT AFTER LINE 36
#    for i in Headers:
#        i = i.strip().replace('"','')
#        Headers_new.append(i)
        
    line_ex = fid.readline()
    line_ex_final = line_ex.split(',') #unpacking the variables
    fid.close() #close the file or terminate the connection with the file system
        
    #Extracting number of data points
    line_ex_final = line_ex[line_ex.find('N='):].partition(',')
    dataPoints = int(line_ex_final[0][2:len(line_ex_final[0])])
    
    # Extracting useless data
    uslessDataLines = line_ex[line_ex.find('E='):].partition(',')
    uslessDataLines = int(uslessDataLines[0][2:len(uslessDataLines[0])])
     
    numberOfLinesInEachTimeStep = numberOfHeaderLines + dataPoints + uslessDataLines
    numberOfTimeSteps = int(math.floor(n/numberOfLinesInEachTimeStep))
    print ("Number of time steps: ",numberOfTimeSteps)
    D = np.zeros([dataPoints, len(Headers), numberOfTimeSteps]) #Initializing np.darray
      
    #Reading select lines from the file
    count = 0
    for ii in range(numberOfTimeSteps):
        print(ii)
        with open(filename) as fid:
#            c = 0
#            for i in range(dataPoints):
#                line = linecache.getline(filename, (ii+1)*numberOfHeaderLines+ii*(uslessDataLines+dataPoints)+i+1)
#                linex = line.split(' ')
#                if (len(linex)!=len(Headers)+1):
#                    c = c + 1
#            if (c > 0):
#                count = count+1
#                continue
            M =  np.loadtxt(fid, dtype = str,
                             delimiter = ' ',
                             skiprows=(ii+1)*numberOfHeaderLines+
                             ii*(uslessDataLines+dataPoints),
                             max_rows = dataPoints)
            if (np.shape(M)[1]!=len(Headers)+1):
                count = count+1
                continue
            D[0:dataPoints,0:len(Headers),ii] = M[:,0:len(M[1])-1]
    print ("Total corrupted time steps: ", count)
    return dataPoints, numberOfTimeSteps, Headers, D

#Generating infiltration values based on precipitation, degree of saturation, porosity, depth of zone
def Qin (precip, ws, n, Zr):
     if (n*Zr*(1 - ws) > precip):
         return precip
     elif (n*Zr*(1-ws) < precip):
         return n*Zr*(1-ws)  

#Generating recharge values based on precipitation, degree of saturation, porosity, depth of zone 
def Qout(ws, wsmin, wsmax, Ks, c, ETmax, n, Zr, A):
    if  (ws <= wsmin):
        return 0
    elif (ws >= wsmax):
        return Ks*(ws**c) + ETmax
    elif (ws < wsmax):
        return Ks*(ws**c) + ws/(n*Zr*A)

#Converting tec data to datafram for heat maps
def Converttomarr (D):
    steps = np.shape(D)[2]
    Headers = np.shape(D)[1]
    size = np.shape(D)[0]
    df2 = np.ndarray([Headers-3, steps, 51,31])
    for i in range(Headers-3):
        counter = 0
        for j in range(steps):
#           print (Headers[i+3]) 
            while counter < (j+1)*size:
                for k in range(51):
                    for l in range(31):
                        df2[i, j, k, l] = D[counter-(j)*size, i+3, j]
                        counter = counter+1
    df = np.ndarray([Headers-3, steps, 51,31])
    for j in range(51):
        df [:,:,50-j,:] = df2 [:,:,j,:]
    return df

#Converting rate data to dataframe for heat maps
def RateConverttomarr (D):
    Trial = np.shape(D)[0]
    Headers = np.shape(D)[2]
    size = np.shape(D)[1]
    df2 = np.ndarray([Headers, Trial, 51,31])
    for i in range(Headers):
        counter = 0
        for j in range(Trial):
#           print (Headers[i+3]) 
            while counter < (j+1)*size:
                for k in range(51):
                    for l in range(31):
                        df2[i,j,k,l] = D[j, counter-(j)*size, i]
                        counter = counter+1
    df = np.ndarray([Headers, Trial, 51,31])
    for j in range(51):
        df [:,:,50-j,:] = df2 [:,:,j,:]
    return df

def calcmassflux (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    mf = np.zeros([len(Trial)*len(gvarnames),8])
    for j in range(len(Trial)):
        di =  d+fpre+str(Trial[j])+fsuf
        print (str(Trial[j]))
#        di+fpre+str(Tforfpre[k])+str(Trial[j])+'_df'
        df = np.load(di+fpre+str(Trial[j])+"_df.npy")
        veliredg = df[2, np.shape(df)[1]-1,yin,xright]
        veliledg = df[2,np.shape(df)[1]-1,yin,xleft]
        veloredg = df[2, np.shape(df)[1]-1,yout,xright]
        veloledg = df[2,np.shape(df)[1]-1,yout,xleft]
        veloelem = df[2,np.shape(df)[1]-1,yout,xleft+1:xright]
        velielem = df[2,np.shape(df)[1]-1,yin,xleft+1:xright]
        velelem = df[2,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]
        vellelem = df[2,np.shape(df)[1]-1,yin+1:yout,xleft]
        velrelem = df[2,np.shape(df)[1]-1,yin+1:yout,xright]
        if (gw ==1):
            satielem = 1
            satoelem = 1
            satlelem = 1
            satrelem = 1
            satiredg = 1
            satiledg = 1
            satoledg = 1
            satoredg = 1
            satelem = 1
        else:
            satielem = df[4,np.shape(df)[1]-1,yin,xleft+1:xright]
            satoelem = df[4,np.shape(df)[1]-1,yout,xleft+1:xright]
            satiredg = df[4,np.shape(df)[1]-1,yin,xright]
            satiledg = df[4,np.shape(df)[1]-1,yin,xright]
            satoledg = df[4,np.shape(df)[1]-1,yout,xleft]
            satoredg = df[4,np.shape(df)[1]-1,yout,xright]
            satlelem = df[4,np.shape(df)[1]-1,yin+1:yout,xleft]
            satrelem = df[4,np.shape(df)[1]-1,yin+1:yout,xright]
            satelem = df[4,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]
        for i in range(len(gvarnames)):
            idx = j*len(gvarnames)+i
            mf[idx,0] = j
            mf[idx,1] = Het[j]
            mf[idx,2] = Anis[j]
            mf[idx,3] = i
            if (gvarnames[i]=="Nitrogen"):
                mf[idx,4] = ((df[vars[3]-3,np.shape(df)[1]-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[3]-3,np.shape(df)[1]-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[3]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem*velielem*velem)) + df[vars[2]-3,np.shape(df)[1]-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[2]-3,np.shape(df)[1]-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[2]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem*velielem*velem))/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx,5] = ((df[vars[3]-3,np.shape(df)[1]-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[3]-3,np.shape(df)[1]-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[3]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem*veloelem*velem)) + df[vars[2]-3,np.shape(df)[1]-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[2]-3,np.shape(df)[1]-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[2]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            else:
                mf[idx,4] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem*velielem*velem))/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx,5] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
#            mf[idx,4] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem*velielem*velem))/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
#            mf[idx,5] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            mf[idx,6] = (mf[idx,4]-mf[idx,5])/mf[idx,4]
            mf[idx,7] = (mf[idx,6] - mf[i,6])/mf[i,6]
    return mf

def calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    sumall = np.zeros([len(Trial)*len(vars),6])
    for j in range(len(Trial)):
        di =  d+fpre+str(Trial[j])+fsuf
        print (str(Trial[j]))
#        di+fpre+str(Tforfpre[k])+str(Trial[j])+'_df'
        df = np.load(di+fpre+str(Trial[j])+"_df.npy")
#        veliredg = df[2, np.shape(df)[1]-1,yin,xright]
#        veliledg = df[2,np.shape(df)[1]-1,yin,xleft]
#        veloredg = df[2, np.shape(df)[1]-1,yout,xright]
#        veloledg = df[2,np.shape(df)[1]-1,yout,xleft]
#        veloelem = df[2,np.shape(df)[1]-1,yout,xleft+1:xright]
#        velielem = df[2,np.shape(df)[1]-1,yin,xleft+1:xright]
#        velelem = df[2,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]
#        vellelem = df[2,np.shape(df)[1]-1,yin+1:yout,xleft]
#        velrelem = df[2,np.shape(df)[1]-1,yin+1:yout,xright]
        if (gw ==1):
            satielem = 1
            satoelem = 1
            satlelem = 1
            satrelem = 1
            satiredg = 1
            satiledg = 1
            satoledg = 1
            satoredg = 1
            satelem = 1
        else:
            satielem = df[4,np.shape(df)[1]-1,yin,xleft+1:xright]
            satoelem = df[4,np.shape(df)[1]-1,yout,xleft+1:xright]
            satiredg = df[4,np.shape(df)[1]-1,yin,xright]
            satiledg = df[4,np.shape(df)[1]-1,yin,xright]
            satoledg = df[4,np.shape(df)[1]-1,yout,xleft]
            satoredg = df[4,np.shape(df)[1]-1,yout,xright]
            satlelem = df[4,np.shape(df)[1]-1,yin+1:yout,xleft]
            satrelem = df[4,np.shape(df)[1]-1,yin+1:yout,xright]
            satelem = df[4,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]
        for i in range(len(gvarnames)):
            idx = j*len(gvarnames)+i
            sumall[idx,0] = Trial[1]+j
            sumall[idx,1] = Het[j]
            sumall[idx,2] = Anis[j]
            sumall[idx,3] = i
            sumall[idx,4] = ((df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*vedge*vedge + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*vedge*velem) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]*satelem))*velem*velem
            sumall[idx,5] = (sumall[idx,4] - sumall[i,4])/sumall[i,4]
    return sumall


def calcmassfluxtime (Trialhead, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    mft = np.ones([len(Trial)*len(vars),8])
    dire = d + fpre + str(Trialhead)
    for j in range(len(Trial)):
        print (str(Trial[j]))
        di =  dire+str(Trial[j])+fsuf
        df = np.load(di+fpre+str(Trialhead)+str(Trial[j])+"_df.npy")
        veliredg = df[2, 1:,yin,xright]
        veliledg = df[2,1:,yin,xleft]
        veloredg = df[2, 1:,yout,xright]
        veloledg = df[2,1:,yout,xleft]
        veloelem = df[2,1:,yout,xleft+1:xright]
        velielem = df[2,1:,yin,xleft+1:xright]
        velelem = df[2,1:,yin+1:yout,xleft+1:xright]
        vellelem = df[2,1:,yin+1:yout,xleft]
        velrelem = df[2,1:,yin+1:yout,xright]
        if (gw ==1):
            satielem = 1
            satoelem = 1
            satlelem = 1
            satrelem = 1
            satiredg = 1
            satiledg = 1
            satoledg = 1
            satoredg = 1
            satelem = 1
        else:
            satiredg = df[4, 1:,yin,xright]
            satiledg = df[4,1:,yin,xleft]
            satoredg = df[4,1:,yout,xright]
            satoledg = df[4,1:,yout,xleft]
            satoelem = df[4,1:,yout,xleft+1:xright]
            satielem = df[4,1:,yin,xleft+1:xright]
            satlelem = df[4,1:,yin+1:yout,xleft]
            satrelem = df[4,1:,yin+1:yout,xright]
            satelem = df[4,1:,yin+1:yout,xleft+1:xright]
        for i in range(len(vars)):
            idx = j*len(vars)+i
            if (type(Trial[1]) == str):
                mft[idx,0] = j
            else:
               mft[idx,0] = Trial[1]+j
            mft[idx,1] = Het[j]
            mft[idx,2] = Anis[j]
            mft[idx,3] = i
            mft[idx,4] = (sum(df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[vars[i]-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)))/sum(np.sum(velielem*velem, axis = -1) + (veliledg+veliredg)*vedge)
            mft[idx,5] = (sum(df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[vars[i]-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)))/sum((np.sum(veloelem*velem, axis = -1) + (veloledg+veloredg)*vedge))
            mft[idx,6] = (mft[idx,4]-mft[idx,5])/mft[idx,4]
            mft[idx,7] = (mft[idx,6] - mft[i,6])/mft[i,6]
    return mft

def calcaggres (mf):
    H = sorted(np.unique(mf[:,1,...]))
    A = sorted(np.unique(mf[...,2]))
    C = sorted(np.unique(mf[...,3]))
    meanresults = []
    stdresults = []
    covresults = []
    for i in H:
        remH = mf[np.where(mf[...,1]==i)]
        for j in A:
            remHA = remH[np.where(remH[...,2]==j)]
            for k in C:
                meanresults.append([np.mean(remHA[np.where(remHA[...,3] == k)][...,6]),
                                np.mean(remHA[np.where(remHA[...,3] == k)][...,7]), i, j, k])
                stdresults.append([np.std(remHA[np.where(remHA[...,3] == k)][...,6]),
                                np.std(remHA[np.where(remHA[...,3] == k)][...,7]), i, j, k])
                covresults.append([np.cov(remHA[np.where(remHA[...,3] == k)][...,7]),
                                np.cov(remHA[np.where(remHA[...,3] == k)][...,7]), i, j, k])
    cleanmresults = [x for x in meanresults if str(x[0]) != 'nan']
    cleansresults = [x for x in stdresults if str(x[0]) != 'nan']
    cleancresults = [x for x in covresults if str(x[0]) != 'nan']
    meanresultsarr = np.array(cleanmresults)
    stdresultsarr = np.array(cleansresults)
    covresultsarr = np.array(cleancresults)
    return meanresultsarr, stdresultsarr, covresultsarr

def calcaggrestime (mf):
    H = sorted(np.unique(mf[:,2,...]))
    A = sorted(np.unique(mf[...,3]))
    C = sorted(np.unique(mf[...,6]))
    Ch = sorted(np.unique(mf[...,11]))
    meanresults = []
    stdresults = []
    covresults = []
    for i in H:
        remH = mf[np.where(mf[...,2]==i)]
        for j in A:
            remHA = remH[np.where(remH[...,3]==j)]
            for k in C:
                remHAC = remHA[np.where(remHA[...,6]==k)]
                for l in Ch:
                    meanresults.append([np.mean(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-3]),
                                np.mean(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-2]), i, j, k,l])
                    stdresults.append([np.std(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-3]),
                                np.std(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-2]), i, j, k,l])
                    covresults.append([np.cov(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-3]),
                                np.cov(remHAC[np.where(remHAC[...,np.shape(mf[1])[0]-1] == l)][...,np.shape(mf[1])[0]-2]), i, j, k,l])
    cleanmresults = [x for x in meanresults if str(x[0]) != 'nan']
    cleansresults = [x for x in stdresults if str(x[0]) != 'nan']
    cleancresults = [x for x in covresults if str(x[0]) != 'nan']
    meanresultsarr = np.array(cleanmresults)
    stdresultsarr = np.array(cleansresults)
    covresultsarr = np.array(cleancresults)
    return meanresultsarr, stdresultsarr, covresultsarr

def calcconcmasstime (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    massendtime = np.zeros([len(vars)])
    massendtimey = np.zeros([51,len(vars)+1])
    massendtimey[:,0] = range(51)
    di =  d+fpre+str(Trial)+fsuf
    print (str(Trial))
    df = np.load(di+fpre+str(Trial)+"_df.npy")
    masstime = np.zeros([np.shape(df)[1], 51, len(vars)+1])
    conctime = np.zeros([np.shape(df)[1], 51, len(vars)+1])
    veliredg = df[2,1:,yin,xright]
    veliledg = df[2,1:,yin,xleft]
    veloredg = df[2,1:,yout,xright]
    veloledg = df[2,1:,yout,xleft]
    veloelem = df[2,1:,yout,xleft+1:xright]
    velielem = df[2,1:,yin,xleft+1:xright]
    velelem = df[2,1:,yin+1:yout,xleft+1:xright]
    vellelem = df[2,1:,yin+1:yout,xleft]
    velrelem = df[2,1:,yin+1:yout,xright]
    if (gw ==1):
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
        for i in range(len(vars)):
#            massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
#            massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
#            massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
#            massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
            masstime[1:,yin,i+1] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright-1]*satielem, axis=-1))*velem*vedge
            masstime[1:,yout,i+1] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright-1]*satoelem, axis = -1))*velem*vedge
            masstime[1:,yin+1:yout, i+1] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
            conctime[1:,yin,i+1] = ((df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg + df[vars[i]-3,1:,yin,xright]*satiredg*veliredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem, axis=-1))*velem)/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
            conctime[1:,yout,i+1] = ((df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg + df[vars[i]-3,1:,yout,xright]*satoredg*veloredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem, axis = -1))*velem)/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
            conctime[1:,yin+1:yout, i+1] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*velem*velelem,axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem*vellelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem*velrelem)*vedge)/(vedge*(vellelem+velrelem) + np.sum(velem*velelem, axis = -1))
    else:
       satielem = df[4,1:,yin,xleft+1:xright]
       satoelem = df[4,1:,yout,xleft+1:xright]
       satlelem = df[4,1:,yin+1:yout,xleft]
       satrelem = df[4,1:,yin+1:yout,xright]
       satiredg = df[4,1:,yin,xright]
       satiledg = df[4,1:,yin,xleft]
       satoledg = df[4,1:,yout,xleft]
       satoredg = df[4,1:,yout,xright]
       satelem = df[4,1:,yin+1:yout,xleft+1:xright]
       for i in range(len(vars)):
           massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[int(np.shape(df)[1])-2])*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[int(np.shape(df)[1])-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
           massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
           massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           masstime[1:,yin,i+1] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
           masstime[1:,yout,i+1] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
           masstime[1:,yin+1:yout, i+1] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
           conctime[1:,yin,i+1] = ((df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1)*velem*vedge)/(vbc*vedge)
           conctime[1:,yout,i+1] = ((df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
           conctime[1:,yin+1:yout, i+1] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge)/(vbc*velem)
    TotalFlow = (veliledg + veloledg + veliredg + veloredg)*vedge + (np.sum(vellelem) + np.sum(velrelem) + np.sum(velelem) + np.sum(velielem) + np.sum(veloelem))*velem
#    Velocity = np.mean ([InVelocity, OutVelocity, MidVelocity])
    Velocity = df[2,np.shape(df)[1]-1,:,:]
    return df, massendtime, masstime, conctime, np.mean(Velocity)

def biomasstime (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    biomassendtime = np.zeros([len(vars)])
    biomassendtimey = np.zeros([51,len(vars)+1])
    biomassendtimey[:,0] = range(51)
    di =  d+fpre+str(Trial)+fsuf
    print (str(Trial))
    df = np.load(di+fpre+str(Trial)+"_df.npy")
    biomasstime = np.zeros([np.shape(df)[1]-1, 51, len(vars)+1])
    bioconctime = np.zeros([np.shape(df)[1]-1, 51, len(vars)+1])
    if (gw ==1):
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
        for i in range(len(vars)):
            biomassendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
            biomassendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
            biomassendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
            biomassendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
            biomasstime[:,yin,i+1] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright-1]*satielem, axis=-1))*velem*vedge
            biomasstime[:,yout,i+1] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright-1]*satoelem, axis = -1))*velem*vedge
            biomasstime[:,yin+1:yout, i+1] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
            bioconctime[:,yin,i+1] = ((df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright-1]*satielem, axis=-1))*velem*vedge)/(vbc*vedge)
            bioconctime[:,yout,i+1] = ((df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
            bioconctime[:,yin+1:yout, i+1] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge)/(vbc*velem)
    else:
       satielem = df[4,1:,yin,xleft+1:xright]
       satoelem = df[4,1:,yout,xleft+1:xright]
       satlelem = df[4,1:,yin+1:yout,xleft]
       satrelem = df[4,1:,yin+1:yout,xright]
       satiredg = df[4,1:,yin,xright]
       satiledg = df[4,1:,yin,xleft]
       satoledg = df[4,1:,yout,xleft]
       satoredg = df[4,1:,yout,xright]
       satelem = df[4,1:,yin+1:yout,xleft+1:xright]
       for i in range(len(vars)):
           biomassendtime[i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[int(np.shape(df)[1])-2])*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[int(np.shape(df)[1])-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomasstime[:,yin,i+1] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
           biomasstime[:,yout,i+1] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
           biomasstime[:,yin+1:yout, i+1] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
           bioconctime[:,yin,i+1] = ((df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1)*velem*vedge)/(vbc*vedge)
           bioconctime[:,yout,i+1] = ((df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
           bioconctime[:,yin+1:yout, i+1] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge)/(vbc*velem)
    return df, biomassendtime, biomasstime, bioconctime

def heatmapconcdist (df, vars, Trial, gvarnames, d, fpre):
    #Heatmaps
    from matplotlib.ticker import FuncFormatter
    fmt = lambda x,pos: '{:4.0f}'.format(x)
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    vel = df[2,np.shape(df)[1]-1,:,:]*-1
    if ('DOC' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(10,10), sharex = True, sharey = True)
        plt.suptitle ("Concentration distribution at steady state (uM)_A-" + str(Trial))
        for i in range(len(vars)):
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt)})
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(vel, square = False, cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[5],
                           cbar_kws={'format': FuncFormatter(fmtvel)})
        ax.flat[5].set_title("Velocity")
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_chem.png"
    elif ('Aerobes' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(7,10), sharex = True, sharey = True)
        plt.suptitle ("Biomass distribution at steady state (uM)_A-" + str(Trial))
        for i in range(len(vars)):
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt)})
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(vel, square = False, cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[3],
                           cbar_kws={'format': FuncFormatter(fmtvel)})
        ax.flat[3].set_title("Velocity")
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass.png"
    pictosave = name.get_figure()
    pictosave.savefig(picname)    
    return (plt)

def heatmapconcdistLOG (df, vars, Trial, gvarnames, d, fpre):
    #Heatmaps
    from matplotlib.ticker import FuncFormatter
    fmt = lambda x,pos: '{:4.0f}'.format(x)
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    vel = df[2,np.shape(df)[1]-1,:,:]*-1
    if ('DOC' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(10,10), sharex = True, sharey = True)
        plt.suptitle ("Concentration distribution at steady state (uM)_A-" + str(Trial))
        for i in range(len(vars)):
            log_norm = LogNorm(vmin=df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min(), vmax=df[vars[i]-3,np.shape(df)[1]-1,:,:].max().max())
            cbar_ticks = [math.pow(10, x) for x in range(int(math.floor(math.log10(df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min()))), 1+int(math.ceil(math.log10(df[vars[i]-3,np.shape(df)[1]-1,:,:].max().max()))))]
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], norm = log_norm, square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt), 'ticks':cbar_ticks}, vmin = df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min(), vmax = df[vars[i]-3,np.shape(df)[1]-1,:,:].min().max().max())
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(vel, square = False, norm = LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()), cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[5],
                           cbar_kws={'format': FuncFormatter(fmtvel), 'ticks':[math.pow(10, x) for x in range(int(math.floor(math.log10(np.min(abs(vel))))), 1+int(math.ceil(math.log10(np.max(abs(vel))))))]},
                           vmin = abs(vel).min().min(), vmax = abs(vel).max().max())
        ax.flat[5].set_title("Velocity")
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_chem.png"
    elif ('Aerobes' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(7,10), sharex = True, sharey = True)
        plt.suptitle ("Biomass distribution at steady state (uM)_A-" + str(Trial))
        for i in range(len(vars)):
            log_norm = LogNorm(vmin=df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min(), vmax=df[vars[i]-3,np.shape(df)[1]-1,:,:].max().max())
            cbar_ticks = [math.pow(10, x) for x in range(int(math.floor(math.log10(df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min()))), 1+int(math.ceil(math.log10(df[vars[i]-3,np.shape(df)[1]-1,:,:].max().max()))))]
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], norm = log_norm, square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt), 'ticks':cbar_ticks}, vmin = df[vars[i]-3,np.shape(df)[1]-1,:,:].min().min(), vmax = df[vars[i]-3,np.shape(df)[1]-1,:,:].min().max().max())
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(vel, square = False, norm = LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()), cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[3],
                            cbar_kws={'format': FuncFormatter(fmtvel), 'ticks':[math.pow(10, x) for x in range(int(math.floor(math.log10(np.min(abs(vel))))), 1+int(math.ceil(math.log10(np.max(abs(vel))))))]},
                           vmin = abs(vel).min().min(), vmax = abs(vel).max().max())
        ax.flat[3].set_title("Velocity")
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass.png"
    pictosave = name.get_figure()
    pictosave.savefig(picname)    
    return (plt)

def heatmapratedist (df, vars, ratenames, nrows, ncols, figsize):
    #Heatmaps
    fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize, sharex = True, sharey = True)
    plt.suptitle ("Impact on respiration rate distribution at steady state (uM)")
    for i in range(len(vars)):
        #    plt.title(gvarnames[i], ax = ax.flat[i])
        name = sns.heatmap(df[i,:,:], square = False, cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[i])
        ax.flat[i].set_title(ratenames[i])
    return name

def heatmapratedistLOG (df, vars, ratenames, nrows, ncols, figsize):
    #Heatmaps
    fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize, sharex = True, sharey = True)
    plt.suptitle ("Impact on respiration rate distribution at steady state (uM)")
    for i in range(len(vars)):
        log_normrate = LogNorm(vmin=abs(df[i,:,:].min().min()), vmax=abs(df[i,:,:].max().max()))
        cbar_ticks = [math.pow(10, x) for x in range(int(math.floor(math.log10(abs(df[i,:,:].min().min())))), 1+int(math.ceil(math.log10(abs(df[i,:,:].max().max())))))]
        name = sns.heatmap(df[i,:,:], square = False, cmap = "YlGnBu",
                           norm = log_normrate,
                           cbar_kws = {"ticks":cbar_ticks}, vmin =abs(df[i,:,:].min().min()), vmax = abs(df[i,:,:].max().max()),
                          xticklabels=False, yticklabels=False, ax=ax.flat[i])
        ax.flat[i].set_title(ratenames[i])
    return name

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plottimeseries (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, inttime, colorstime, nrows, ncols, figsize, chem, lylim, uylim):
    col=0
    #figbig = plt.figure()
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(inttime)):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[inttime[k],0:51,gvarnames.index(chem)+1], label = inttime[k]*50, color = colorstime[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
            ax[int(math.floor(col/ncols))][col%ncols].legend()
            ax[int(math.floor(col/ncols))][col%ncols].set_ylim(lylim, uylim)
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+ "_" + chem + "_" + datetime.now().strftime("%d%m%Y%H%M")+".png"
    plt.savefig(picname)
    return(plt)

def concprofileatss (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors):
    for i in intsce:
        if ('DOC' in gvarnames):
            df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            #    print (np.mean(df[6, np.shape(df)[1]-1, yin, :]), np.mean(df[6, np.shape(df)[1]-1, yout-1, :]), np.mean(df[6, np.shape(df)[1]-1, yout, :]))
            fig, host = plt.subplots()
            fig.subplots_adjust(right=0.75)
            
            par1 = host.twinx()
            par2 = host.twinx()
            
            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["right"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["right"].set_visible(True)
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,1], label = gvarnames[0], color = colors[0])
            p2, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,2], label = gvarnames[1], color = colors[1])
            p3, = par1.plot(conctime[np.shape(conctime)[0]-1,0:51,3], label = gvarnames[2], color = colors[2])
            p4, = par2.plot(conctime[np.shape(conctime)[0]-1,0:51,4], label = gvarnames[3], color = colors[3])
            
            host.set_xlim(0, 51)
            host.set_ylim(0, 800)
            par1.set_ylim(30, 60)
            par2.set_ylim(50, 260)
            
            host.set_xlabel("Y")
            host.set_ylabel("DOC, DO (uM)")
            par1.set_ylabel(str(gvarnames[2])+" (uM)")
            par2.set_ylabel(str(gvarnames[3])+" (uM)")
            
            host.yaxis.label.set_color(p1.get_color())
            par1.yaxis.label.set_color(p3.get_color())
            par2.yaxis.label.set_color(p4.get_color())
            
            tkw = dict(size=4, width=1.5)
            host.tick_params(axis='y', colors=p1.get_color(), **tkw)
            par1.tick_params(axis='y', colors=p3.get_color(), **tkw)
            par2.tick_params(axis='y', colors=p4.get_color(), **tkw)
            host.tick_params(axis='x', **tkw)
            
            lines = [p1, p2, p3, p4]
            
            host.legend(lines, [l.get_label() for l in lines])
            plt.title (Trial[Trial.index(i)])
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_chem_concatss.png"
        
        elif ('Aerobes' in gvarnames):
            df,massendtime, masstime, conctime = biomasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            fig, host = plt.subplots()
            fig.subplots_adjust(right=0.75)
            
            par1 = host.twinx()
            par2 = host.twinx()
            
            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["right"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["right"].set_visible(True)
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,1], label = gvarnames[0], color = colors[0])
            p2, = par1.plot(conctime[np.shape(conctime)[0]-1,0:51,2], label = gvarnames[1], color = colors[2])
            p3, = par2.plot(conctime[np.shape(conctime)[0]-1,0:51,3], label = gvarnames[2], color = colors[3])
            
            host.set_xlim(0, 51)
            host.set_ylim(0, 500)
            par1.set_ylim(0, 30)
            par2.set_ylim(0, 30)
            
            host.set_xlabel("Y")
            host.set_ylabel(str(gvarnames[0])+" (uM)")
            par1.set_ylabel(str(gvarnames[1])+" (uM)")
            par2.set_ylabel(str(gvarnames[2])+" (uM)")
            
            host.yaxis.label.set_color(p1.get_color())
            par1.yaxis.label.set_color(p2.get_color())
            par2.yaxis.label.set_color(p3.get_color())
            
            tkw = dict(size=4, width=1.5)
            host.tick_params(axis='y', colors=p1.get_color(), **tkw)
            par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
            par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
            host.tick_params(axis='x', **tkw)
            
            lines = [p1, p2, p3]
            
            host.legend(lines, [l.get_label() for l in lines])
            plt.title (Trial[Trial.index(i)])
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass_concatss.png"
        plt.savefig(picname)
    return None

def plotconcallatss(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors, nrows, ncols, figsize):
    col=0
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(vars)-2):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[np.shape(conctime)[0]-1,0:51,k+1], label = gvarnames[k], color = colors[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+"_" +datetime.now().strftime("%d%m%Y%H%M")+"_concatss.png"
    plt.savefig(picname)
    
    return None

def plotconcwithhet (Regimes, intsce, chem, colorfamily, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    colors = sns.color_palette(colorfamily, len(intsce))
    fig, axes = plt.subplots(nrows = 1, ncols = len(Regimes), figsize = [16,5])
    for r in Regimes:
        d = r"X:/Saturated_flow/Anaerobic/" + r+ "AR_changedkindox/"
        lines = []
        for i in intsce:
            df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            p1, = axes.flat[Regimes.index(r)].plot(conctime[np.shape(conctime)[0]-1,0:51,gvarnames.index(chem)+1],
                           label = str(Het[Trial.index(i)]) + ":" + str(Anis[Trial.index(i)]), color = colors[intsce.index(i)])
            axes.flat[Regimes.index(r)].set_ylim(0, 250)
            axes.flat[Regimes.index(r)].set_xlabel("Y (cm)")
            axes.flat[Regimes.index(r)].set_ylabel("DO (uM)")
            lines.append(p1,)
                
            axes.flat[Regimes.index(r)].legend(lines, [l.get_label() for l in lines])
            if (r == "Equal"):
                axes.flat[Regimes.index(r)].set_title ("Medium")
            else:
                axes.flat[Regimes.index(r)].set_title (r)
    return fig

def plotoutflowtimeseries (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, inttime, colorstime, nrows, ncols, figsize, chem, lylim, uylim, start, end):
    col=0
    #figbig = plt.figure()
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[start:end,yout,gvarnames.index(chem)+1], label = "Tracer", color = "black")
        ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
        #        ax[int(math.floor(col/ncols))][i%ncols].set_xlabel("Y")
        #        ax[int(math.floor(col/ncols))][i%ncols].set_ylabel("Concentration (uM)")
        ax[int(math.floor(col/ncols))][col%ncols].set_ylim(lylim, uylim)
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+ "_" + chem + "_time.png"
    plt.savefig(picname)
    return(plt)
    
def plotconcallatsstime(Trialhead, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors, nrows, ncols, figsize):
    col=0
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        Trial2 = str(Trialhead) + str(i)
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial2[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(vars)-2):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[np.shape(conctime)[0]-1,0:51,k+1], label = gvarnames[k], color = colors[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+"_" +datetime.now().strftime("%d%m%Y%H%M")+"_concatss.png"
    plt.savefig(picname)
    
    return None

def processchembiomassfiles(biomassfile, chemfiles, regime):
    biomass = pd.read_csv(biomassfile, delimiter = "\t")
    biomass['Chemb'] = biomass['Chem']
    biomass['VA']=biomass['Variance']*biomass['Anisotropy']
    chemfile = pd.read_csv(chemfiles, delimiter = "\t")
    oxy = pd.merge(biomass[biomass['Chemb']=='Aerobes'], chemfile[chemfile['Chem']=='DO'], on = 'Trial', suffixes=("_b","_c"))
    amm = pd.merge(biomass[biomass['Chemb']=='Ammonia oxidizers'], chemfile[chemfile['Chem']=='Ammonium'], on = 'Trial', suffixes=("_b","_c"))
    nitra = pd.merge(biomass[biomass['Chemb']=='Nitrate reducers'], chemfile[chemfile['Chem']=='Nitrate'], on = 'Trial', suffixes=("_b","_c"))
    allbiomass = pd.concat([oxy,amm,nitra], axis=0)
    allbiomass['delbiomass%']=allbiomass['Change_umoles']*100
    allbiomass['del2massflux%']=allbiomass['del2massflux']*10
    allbiomass['Regime']=regime
    
    return allbiomass

def processchemfiles(chemfile, regime):
    chemfile = pd.read_csv(chemfile, delimiter = "\t")
    chemfile['Regime']=regime
    chemfile['VA']=chemfile['Variance']*chemfile['Anisotropy']
    chemfile['delmassflux%']=chemfile['delmassflux']*100
    chemfile['del2massflux%']=chemfile['del2massflux']*100
    
    return chemfile

def scatterbioflux (dataset1, dataset2, dataset3):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [16,10])
    plt.suptitle("Change in removal of carbon and nitrogen with respect to change in biomass at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.scatterplot(x = "delbiomass%", y = "del2massflux%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 270,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax, row in zip(axes[:,0], Chems):
        if (row == "Aerobes"):
            ax.annotate("DO", xy = (0, 0.5), xytext = (-ax.yaxis.labelpad - 70,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'right', va = 'center')
        if (row == "Ammonia oxidizers"):
            ax.annotate("Ammonium", xy = (0, 0.5), xytext = (-ax.yaxis.labelpad - 70,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'right', va = 'center')
        if (row == "Nitrate reducers"):
            ax.annotate("Nitrate", xy = (0, 0.5), xytext = (-ax.yaxis.labelpad - 70,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'right', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in biomass (%)")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxVAflux(dataset1, dataset2, dataset3, chemseries):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [16,10])
    plt.suptitle("Change in removal of carbon and nitrogen at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "VA", y = "del2massflux%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 270,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Variance x Anisotropy (-)")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxVAbio(dataset1, dataset2, dataset3):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [16,10])
    plt.suptitle("Change in biomass at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "VA", y = "delbiomass%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 270,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Variance x Anisotropy (-)")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxV_Aflux(dataset1, dataset2, dataset3, chemseries, imgsize):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = imgsize)
    plt.suptitle("Change in removal of carbon and nitrogen at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "del2massflux%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 400,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[4]:
        ax.set_xlabel("Variance : Anisotropy")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxV_Abio(dataset1, dataset2, dataset3, imgsize):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
#    dfall = pd.concat([equalbc, slowbc, fastbc], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance_b'][i]) + ":" + str(dfall['Anisotropy_b'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance_b','Anisotropy_b'])
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = imgsize)
    plt.suptitle("Change in total biomass at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "delbiomass%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 400,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Variance : Anisotropy")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def tracerstudies ():
    di = r"X:\Saturated_flow\Steady_state\Tracer_studies\/"
    filename = "NS-ATracerH_84_MF_1081220191522_equal.csv"
    equal = pd.read_csv(di+filename, delimiter = "\t")
    equal['Regime'] = "Medium"

    filename = "NS-ATracerH_84_MF_1121220192228_fast_0.00002.csv"
    fast = pd.read_csv(di+filename, delimiter = "\t")
    fast['Regime'] = "Fast"
    
    filename = "NS-ATracerH_84_MF_1081220191539_slowsk.csv"
    slow = pd.read_csv(di+filename, delimiter = ",")
    slow['Regime'] = "Slow"
    
    breakthrough = pd.concat([fast, slow, equal], axis=0, ignore_index = True)
    breakthrough ['Heterogeneity'] = breakthrough ['Variance'] 
    breakthrough['VA']=breakthrough['Heterogeneity']*breakthrough['Anisotropy']
    
    breakthrough['%ofhomogeneous']=breakthrough['del']*100
#    bth = breakthrough.rename(columns={'Scenario':'Trial'}).astype(str)
    
    l = []
    for i in range(len(breakthrough)):
        l.append(str(breakthrough['Variance'][i]) + ":" + str(breakthrough['Anisotropy'][i]))
        
    breakthrough['Xlabels'] = l
    breakthrough = breakthrough.sort_values(by=['Variance','Anisotropy'])
    
    return breakthrough

def boxrestime_flux(dataset1, dataset2, dataset3, chemseries):
    
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    bins = [-60, -40, -20, 0, 20]
    dfall2['binned'] = pd.cut(dfall2['%ofhomogeneous'].astype(float), bins)
    
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
 
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10])
    plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state")
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "binned", y = "del2massflux%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in residence time (%)")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxrestime_biomass(dataset1, dataset2, dataset3):
    
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance_b'][i]) + ":" + str(dfall['Anisotropy_b'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance_b','Anisotropy_b'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    bins = [-60, -40, -20, 0, 20]
    dfall2['binned'] = pd.cut(dfall2['%ofhomogeneous'].astype(float), bins)
    
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
 
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10])
    plt.suptitle("Change in biomass with respect to homogeneous scenario at steady state")
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "binned", y = "delbiomass%", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in residence time (%)")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def scatterrestime_flux(dataset1, dataset2, dataset3, chemseries):
    
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    bins = [-60, -40, -20, 0, 20]
#    dfall2['binned'] = pd.cut(dfall2['%ofhomogeneous'].astype(float), bins)
    
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)
 
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10])
    plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state")
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            axes[colidx1][colidx2].scatter(x = "%ofhomogeneous", y = "del2massflux%", c = colseries[Regimes.index(i)], data = dfc)
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            axes[colidx1][colidx2].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
            axes[colidx1][colidx2].set_xticklabels([])
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[-1]:
        ax.set_xlabel("Relative difference in residence time (%)")
        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10))
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def scatterrestime_biomass(dataset1, dataset2, dataset3):
    
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance_b'][i]) + ":" + str(dfall['Anisotropy_b'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance_b','Anisotropy_b'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    bins = [-60, -40, -20, 0, 20]
#    dfall2['binned'] = pd.cut(dfall2['%ofhomogeneous'].astype(float), bins)
    
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)
 
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10], sharex='col')
    plt.suptitle("Change in biomass with respect to homogeneous scenario at steady state")
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            axes[colidx1][colidx2].scatter("%ofhomogeneous", "delbiomass%", c = colseries[Regimes.index(i)], data = dfc,
                cmap = colseries[Regimes.index(i)])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            axes[colidx1][colidx2].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
            axes[colidx1][colidx2].set_xticklabels([])
#            axes[colidx1][colidx2].set_xticklabels([20, 0, -20, -40, -60, -80])
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[-1]:
        ax.set_xlabel("Relative difference in residence time (%)")
        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10))
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def calcoxiccells(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    oxiccells= np.zeros([len(Trial), 51,1])
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        c = []
        for k in range(51):
            c.append(np.count_nonzero(df[vars[1]-3,np.shape(df)[1]-1,k,:]>limit))
#        print (Trial[j], c)
        print ("Calculating number of oxic cells")
        oxiccells [j,:,0] = c
    
    return oxiccells

def plotoxiccellssdo(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    oxiccells = calcoxiccells(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    nrows = 7
    ncols = 7
    figsize = [21,28]
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True, sharey = True)
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity = calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        axes.flat[j].plot(conctime[-1,:,2], 'r-')
#        axes.flat[j].set_yticks([])
        ax2 = axes.flat[j].twinx()
        ax2.plot(oxiccells[j,:,0],'b-')
#        ax2.set_yticks([])
#        ax2.set_yticklabels([])
        axes.flat[j].set_title(Trial[j])
    for ax in axes[:,0]:
        ax.set_ylabel("DO (uM)", color='r')
#        ax.set_yticks([])
    for ax in axes[:,-1]:
        ax.twinx().set_ylabel("Number of oxic cells (-)", color='b')
#        ax.twinx().set_yticklabels(labels = [0,10,20,30],color = 'b')
    for ax in axes[-1]:
        ax.set_xlabel("Y (cm)")
    return fig
        
