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
    D = np.zeros([dataPoints, len(Headers), numberOfTimeSteps]) #Initializing np.darray
      
    #Reading select lines from the file
    count = 0
    for ii in range(numberOfTimeSteps):
#        c = 0
#        for i in range(750):
#            line = linecache.getline(filename, (ii+1)*numberOfHeaderLines+ii*(uslessDataLines+dataPoints)+i+1)
#            linex = line.split(' ')
#            if (len(linex)!=len(Headers)+1):
#                c = c + 1
#        if (c > 0):
#            print(ii, ":", c)
#            count = count+1
#            continue
        with open(filename) as fid:
            M =  np.loadtxt(fid, dtype = str,
                            delimiter = ' ',
                            skiprows=(ii+1)*numberOfHeaderLines+
                            ii*(uslessDataLines+dataPoints),
                            max_rows = dataPoints)
            if (np.shape(M)[1]!=len(Headers)+1):
                count = count+1
                continue
            elif (np.shape(M)[0]!=dataPoints):
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
    por = 0.2
    doc1=10-gw
    Bmo1 = 9 - gw
    Bmn1 = 16 - gw
    Bms1 = 21 - gw
    Bma1 = 26 - gw
    Bimo1 = 14 - gw
    Bimn1 = 19 - gw
    Bims1 = 24 - gw
    Bima1 = 28 - gw
    POM1 = 30 - gw
    Amm1=12-gw
    nitra1=17-gw
    Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    mf = np.zeros([len(Trial)*len(gvarnames),9])
    for j in range(len(Trial)):
        di =  d+fpre+str(Trial[j])+fsuf
        print (str(Trial[j]))
#        di+fpre+str(Tforfpre[k])+str(Trial[j])+'_df'
        df = np.load(di+fpre+str(Trial[j])+"_df.npy")
        veliredg = df[2,-1,yin,xright]
        veliledg = df[2,-1,yin,xleft]
        veloredg = df[2,-1,yout,xright]
        veloledg = df[2,-1,yout,xleft]
        veloelem = df[2,-1,yout,xleft+1:xright]
        velielem = df[2,-1,yin,xleft+1:xright]
        velelem = df[2,-1,yin+1:yout,xleft+1:xright]
        vellelem = df[2,-1,yin+1:yout,xleft]
        velrelem = df[2,-1,yin+1:yout,xright]
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
            satielem = df[4,-1,yin,xleft+1:xright]
            satoelem = df[4,-1,yout,xleft+1:xright]
            satiredg = df[4,-1,yin,xright]
            satiledg = df[4,-1,yin,xright]
            satoledg = df[4,-1,yout,xleft]
            satoredg = df[4,-1,yout,xright]
            satlelem = df[4,-1,yin+1:yout,xleft]
            satrelem = df[4,-1,yin+1:yout,xright]
            satelem = df[4,-1,yin+1:yout,xleft+1:xright]
        for i in range(len(gvarnames)):
            idx = j*len(gvarnames)+i
            mf[idx,0] = j
            mf[idx,1] = Het[j]
            mf[idx,2] = Anis[j]
            mf[idx,3] = i
            if (gvarnames[i]=="Nitrogen"):
                ninlet = 0
                noutlet = 0
                for n in Nspecies:
                    ninlet = ninlet + (df[n-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[n-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[n-3,-1,yin,xleft+1:xright]*satielem*velielem*velem))/por
                    noutlet = noutlet + (df[n-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[n-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[n-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/por
                ninlet = ninlet/10 + (df[Amm1-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[Amm1-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[Amm1-3,-1,yin,xleft+1:xright]*satielem*velielem*velem) + df[nitra1-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[nitra1-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[nitra1-3,-1,yin,xleft+1:xright]*satielem*velielem*velem))/por
                noutlet = noutlet/10 + (df[Amm1-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[Amm1-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[Amm1-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem) + df[nitra1-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[nitra1-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[nitra1-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/por
                mf[idx,4] = abs(ninlet)#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx,5] = abs(noutlet)#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            elif (gvarnames[i]=="TOC"):
                cinlet = 0
                coutlet = 0
                for c in Cspecies:
                    cinlet = cinlet + (df[c-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[c-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[c-3,-1,yin,xleft+1:xright]*satielem*velielem*velem))/por
                    coutlet = coutlet + (df[c-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[c-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[c-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem))/por
                mf[idx,4] = abs(cinlet)#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx,5] = abs(coutlet)#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            else:
                mf[idx,4] = abs(((df[vars[i]-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[i]-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[i]-3,-1,yin,xleft+1:xright]*satielem*velielem*velem)))/por)#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx,5] = abs(((df[vars[i]-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[i]-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[i]-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem)))/por)#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            #calculating removal in mass
            mf[idx,6] = (mf[idx,4]-mf[idx,5])
            #normalizing removal in mass with incoming mass
            mf[idx,7] = (mf[idx,4]-mf[idx,5])/mf[idx,4]
            #comparing removal with homogeneous case
            mf[idx,8] = mf[idx,6]/mf[i,6]
    return mf

def calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    por = 0.2
    sumall = np.zeros([len(Trial)*len(vars),9])
    for j in range(len(Trial)):
        di =  d+fpre+str(Trial[j])+fsuf
        print (str(Trial[j]))
        df = np.load(di+fpre+str(Trial[j])+"_df.npy")
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
            sumall[idx,4] = (((df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*vedge*vedge + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*vedge*velem) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]*satelem))*velem*velem)*por
            sumall[idx,5] = (sumall[idx,4])/sumall[i,4]
        total = np.sum(sumall[idx-11:idx+1,4])
        for k in range(len(gvarnames)):
            idxk = j*len(gvarnames)+k
            sumall[idxk,6] = sumall[idxk,4]/total
            #comparing the fraction of species with homogeneous case in the uniform flow rate scenario
            sumall[idxk,7] = sumall[idxk,6]/sumall[i,6]
            #comparing the fraction of species with homogeneous case in the varying flow rate scenario
            sumall[idxk,8] = sumall[idxk,6]/sumall[Trial.indexd(k)*len(gvarnames)+k,6]
    return sumall

def calcmft_temp(Tforfpre, Trial, vars, gvarnames, biomassvars, biomassgvarnames, d, fpre, fsuf, Het, Anis, gw):
    yout = 50
    yin = 0
    xleft = 0
    xright = 30
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    doc1=10-gw
    por = 0.2
    Bmo1 = 9 - gw
    Bmn1 = 16 - gw
    Bms1 = 21 - gw
    Bma1 = 26 - gw
    Bimo1 = 14 - gw
    Bimn1 = 19 - gw
    Bims1 = 24 - gw
    Bima1 = 28 - gw
    POM1 = 30 - gw
    Amm1=12-gw
    nitra1=17-gw
    Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    mft = np.zeros([(len(Trial)*len(Tforfpre))*len(gvarnames),11])
    #dire = d + fpre + str(Trialhead)
    sumall = np.zeros([(len(Trial)*len(Tforfpre))*len(biomassgvarnames),13])
    for j in Tforfpre:
        print (j)
        di =  d+j
        for k in Trial: 
            df = np.load(di+fpre+str(k)+fsuf+fpre+str(k)+"_df.npy")
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
            for i in range(len(gvarnames)):
                idx = (Tforfpre.index(j)*len(Trial)+Trial.index(k))*len(gvarnames)+i
                #            print (j,k, ": ", idx)
                if (type(k) == str):
                    mft[idx,0] = 0
                else:
                    mft[idx,0] = Trial[Trial.index(k)]
                mft[idx,1] = Het[Trial.index(k)]
                mft[idx,2] = Anis[Trial.index(k)]
                mft[idx,3] = i
                if (gvarnames[i]=="Nitrogen"):
                    ninlet = 0
                    noutlet = 0
                    for n in Nspecies:
                        ninlet = ninlet + (sum(df[n-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[n-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[n-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)))/por
                        noutlet = noutlet + (sum(df[n-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[n-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[n-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)))/por
                    ninlet = ninlet/10 + (sum(df[Amm1-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[Amm1-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[Amm1-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)) + sum(df[nitra1-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[nitra1-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[nitra1-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)))/por
                    noutlet = noutlet/10 + (sum(df[Amm1-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[Amm1-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[Amm1-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)) + sum(df[nitra1-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[nitra1-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[nitra1-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)))/por
                    mft[idx,4] = ninlet#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                    mft[idx,5] = noutlet#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
                elif (gvarnames[i]=="TOC"):
                    cinlet = 0
                    coutlet = 0
                    for c in Cspecies:
                        cinlet = cinlet + (sum(df[c-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[c-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[c-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)))/por
                        coutlet = coutlet + (sum(df[c-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[c-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[c-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)))/por
                    mft[idx,4] = cinlet#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                    mft[idx,5] = coutlet#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)      
                else:
                    mft[idx,4] = (sum(df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg*vedge) + sum(df[vars[i]-3,1:,yin,xright]*satiredg*veliredg*vedge) + sum(np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1)))/por
                    mft[idx,5] = (sum(df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg*vedge) + sum(df[vars[i]-3,1:,yout,xright]*satoredg*veloredg*vedge) + sum(np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1)))/por
                #calculating removal in mass
                mft[idx,6] = (mft[idx,4]-mft[idx,5])
                #comparing removal with homogeneous case in the uniform flow rate scenario
                mft[idx,7] = (mft[idx,6])/mft[i,6]
                #comparing removal with heterogeneous case in the uniform flow rate scenario
                mft[idx,8] = (mft[idx,6])/mft[Trial.index(k)*len(gvarnames)+i,6]
                #time series index
                mft[idx,9] = Tforfpre.index(j)
                #comparing removal with homogeneous case in the varying flow rate scenario
                mft[idx,10] = (mft[idx,6])/mft[Tforfpre.index(j)*len(Trial)*len(gvarnames)+i,6]
            total=0
            for b in range(len(biomassgvarnames)):
                idxb = (Tforfpre.index(j)*len(Trial)+Trial.index(k))*len(biomassgvarnames)+b
                if (type(k) == str):
                    sumall[idxb,0] = 0
                else:
                    sumall[idxb,0] = Trial[Trial.index(k)]
                sumall[idxb,1] = Het[Trial.index(k)]
                sumall[idxb,2] = Anis[Trial.index(k)]
                sumall[idxb,3] = b
                #total biomass in the domain
                sumall[idxb, 4] = (por*(df[biomassvars[b]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[biomassvars[b]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[biomassvars[b]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[biomassvars[b]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*vedge*vedge + (sum(df[biomassvars[b]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem) + sum(df[biomassvars[b]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem) + sum(df[biomassvars[b]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[biomassvars[b]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*vedge*velem) + por*(sum(sum(df[biomassvars[b]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]*satelem))*velem*velem)
                total = total+sumall[idxb,b]
                #comparing biomass with homogeneous case in the uniform flow rate scenario
                sumall[idxb,5] = (sumall[idxb,4])/sumall[b,4]
                #comparing biomass with heterogeneous case in the uniform flow rate scenario
                sumall[idxb,6] = (sumall[idxb,4])/sumall[Trial.index(k)*len(biomassgvarnames)+b,4]
                #comparing biomass with homogenous case in the varying flow rate scenario
                sumall[idxb,7] = (sumall[idxb,4])/sumall[Tforfpre.index(j)*len(Trial)*len(biomassgvarnames)+b,4]
                #time series index
                sumall[idxb,8] = Tforfpre.index(j)
#            total = np.sum(sumall[idxb-11:idxb+1,4])
#            print (k,j)
            for g in range(len(biomassgvarnames)):
                idxk = (Tforfpre.index(j)*len(Trial)+Trial.index(k))*len(biomassgvarnames)+g
                #calculating fraction of species with total biomass in the domain
                sumall[idxk,9] = sumall[idxk,4]/total
                #comparing the fraction of species with homogeneous case in the uniform flow rate scenario
                sumall[idxk,10] = sumall[idxk,9]/sumall[g,9]
                #comparing the fraction of species with heterogeneous case in the uniform flow rate scenario
                sumall[idxk,11] = sumall[idxk,9]/sumall[Trial.index(k)*len(biomassgvarnames)+g,9]
                #comparing the fraction of species with homogeneous case in the varying flow rate scenario
                sumall[idxk,12] = sumall[idxk,9]/sumall[Tforfpre.index(j)*len(Trial)*len(biomassgvarnames)+g,9]
    return mft, sumall 

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

def calcconcmasstime (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    por = 0.2
    doc1=10-gw
    Bmo1 = 9 - gw
    Bmn1 = 16 - gw
    Bms1 = 21 - gw
    Bma1 = 26 - gw
    Bimo1 = 14 - gw
    Bimn1 = 19 - gw
    Bims1 = 24 - gw
    Bima1 = 28 - gw
    POM1 = 30 - gw
    Amm1=12-gw
    nitra1=17-gw
    Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    massendtime = np.zeros([len(gvarnames)])
    massendtimey = np.zeros([51,len(gvarnames)])
    massendtimey[:,0] = range(51)
    di =  d+fpre+str(Trial)+fsuf
    print (str(Trial))
    df = np.load(di+fpre+str(Trial)+"_df.npy")
    masstime = np.zeros([np.shape(df)[1], 51, len(gvarnames)])
    conctime = np.zeros([np.shape(df)[1], 51, len(gvarnames)])
    veliredg = df[2,1:,yin,xright]
    veliledg = df[2,1:,yin,xleft]
    veloredg = df[2,1:,yout,xright]
    veloledg = df[2,1:,yout,xleft]
    veloelem = df[2,1:,yout,xleft+1:xright]
    velielem = df[2,1:,yin,xleft+1:xright]
    velelem = df[2,1:,yin+1:yout,xleft+1:xright]
    vellelem = df[2,1:,yin+1:yout,xleft]
    velrelem = df[2,1:,yin+1:yout,xright]
    if (gw == 1):
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
        for i in range(len(gvarnames)):
            if (gvarnames[i]=="Nitrogen"):
                ninlet = 0
                noutlet = 0
                for n in Nspecies:
                    ninlet = ninlet + (df[n-3,1:,yin,xleft]*satiledg*veliledg*vedge + df[n-3,1:,yin,xright]*satiredg*veliredg*vedge + np.sum(df[n-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1))/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
                    noutlet = noutlet + (df[n-3,1:,yout,xleft]*satoledg*veloledg*vedge + df[n-3,1:,yout,xright]*satoredg*veloredg*vedge + np.sum(df[n-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1))/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
                conctime[1:,yin,i] = ninlet/10 + (df[Amm1-3,1:,yin,xleft]*satiledg*veliledg*vedge + df[Amm1-3,1:,yin,xright]*satiredg*veliredg*vedge + np.sum(df[Amm1-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1) + df[nitra1-3,1:,yin,xleft]*satiledg*veliledg*vedge + df[nitra1-3,1:,yin,xright]*satiredg*veliredg*vedge + np.sum(df[nitra1-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1))/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
                conctime[1:,yout,i] = noutlet/10 + (df[Amm1-3,1:,yout,xleft]*satoledg*veloledg*vedge + df[Amm1-3,1:,yout,xright]*satoredg*veloredg*vedge + np.sum(df[Amm1-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1) + df[nitra1-3,1:,yout,xleft]*satoledg*veloledg*vedge + df[nitra1-3,1:,yout,xright]*satoredg*veloredg*vedge + np.sum(df[nitra1-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1))/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
            elif (gvarnames[i]=="TOC"):
                cinlet = 0
                coutlet = 0
                for c in Cspecies:
                    cinlet = cinlet + (df[c-3,1:,yin,xleft]*satiledg*veliledg*vedge + df[c-3,1:,yin,xright]*satiredg*veliredg*vedge + np.sum(df[c-3,1:,yin,xleft+1:xright]*satielem*velielem*velem, axis = -1))/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
                    coutlet = coutlet + (df[c-3,1:,yout,xleft]*satoledg*veloledg*vedge + df[c-3,1:,yout,xright]*satoredg*veloredg*vedge + np.sum(df[c-3,1:,yout,xleft+1:xright]*satoelem*veloelem*velem, axis = -1))/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
                conctime[1:,yin,i] = cinlet#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                conctime[1:,yout,i] = coutlet#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)      
            else:
#                massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
#                massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
#                massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
#                massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
                masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
                masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
                masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
                conctime[1:,yin,i] = ((df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg + df[vars[i]-3,1:,yin,xright]*satiredg*veliredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem, axis=-1))*velem)/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
                conctime[1:,yout,i] = ((df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg + df[vars[i]-3,1:,yout,xright]*satoredg*veloredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem, axis = -1))*velem)/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
#                conctime[1:,yin+1:yout, i] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*velem*velelem,axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem*vellelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem*velrelem)*vedge)/(vedge*(vellelem+velrelem) + np.sum(velem*velelem, axis = -1))
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
           massendtimey[yin,i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
           massendtimey[yout,i] = (df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
           massendtimey[yin+1:yout, i] = sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
           masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
           masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
           conctime[1:,yin,i] = ((df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1)*velem*vedge)/(vbc*vedge)
           conctime[1:,yout,i] = ((df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
           conctime[1:,yin+1:yout, i] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge)/(vbc*velem)
    TotalFlow = (veliledg + veloledg + veliredg + veloredg)*vedge + (np.sum(vellelem) + np.sum(velrelem) + np.sum(velelem) + np.sum(velielem) + np.sum(veloelem))*velem
#    Velocity = np.mean ([InVelocity, OutVelocity, MidVelocity])
    Velocity = df[2,np.shape(df)[1]-1,:,:]
    Headinlettime = np.mean(df[2,1:,yin,:], axis = -1)*-1
    return df, massendtime, masstime, conctime, np.mean(Velocity), Headinlettime

def calcconcmasstimeX (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    massendtime = np.zeros([len(vars)])
    massendtimey = np.zeros([51,len(vars)+1])
    massendtimey[:,0] = range(51)
    di =  d+fpre+str(Trial)+fsuf
    print (str(Trial))
    df = np.load(di+fpre+str(Trial)+"_df.npy")
    masstime = np.zeros([np.shape(df)[1], 31, len(vars)+1])
    conctime = np.zeros([np.shape(df)[1], 31, len(vars)+1])
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
            masstime[1:,xleft,i+1] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yout,xleft]*satoledg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft]*satielem, axis=-1))*velem*vedge
            masstime[1:,xright,i+1] = (df[vars[i]-3,1:,yin,xright]*satiredg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin+1:yout,xright]*satoelem, axis = -1))*velem*vedge
            masstime[1:,xleft+1:xright, i+1] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-2) + ((df[vars[i]-3,1:,yin,xleft+1:xright]*satlelem + df[vars[i]-3,1:,yout,xleft+1:xright])*satrelem)*velem*vedge
            conctime[1:,xleft,i+1] = ((df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg + df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg)*(vedge) + (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem*vellelem, axis=-1))*velem)/(vedge*(veloledg+veliledg) + np.sum(velem*vellelem, axis = -1))
            conctime[1:,xright,i+1] = ((df[vars[i]-3,1:,yin,xright]*satiredg*veliredg+ df[vars[i]-3,1:,yout,xright]*satoredg*veloredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yin+1:yout,xright]*satrelem*velrelem, axis = -1))*velem)/(vedge*(veloredg+veliredg) + np.sum(velem*velrelem, axis = -1))
            conctime[1:,xleft+1:xright, i+1] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*velem*velelem,axis=-2) + (df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem + df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem)*vedge)/(vedge*(veloelem+velielem) + np.sum(velem*velelem, axis = -2))
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

def biomasstimefunc (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, biomassvars):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    biomassendtime = np.zeros([len(biomassvars)])
    biomassendtimey = np.zeros([51,len(biomassvars)])
    biomassendtimey[:,0] = range(51)
    di =  d+fpre+str(Trial)+fsuf
    print (str(Trial))
    df = np.load(di+fpre+str(Trial)+"_df.npy")
    biomasstime = np.zeros([np.shape(df)[1]-1, 51, len(biomassvars)])
    bioconctime = np.zeros([np.shape(df)[1]-1, 51, len(biomassvars)])
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
        for i in range(len(biomassvars)):
            biomassendtime[i] = (df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
            biomassendtimey[yin,i] = (df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright]*satielem))*velem*vedge
            biomassendtimey[yout,i] = (df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright]*satoelem))*velem*vedge
            biomassendtimey[yin+1:yout, i] = sum(sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright]*satelem*(velem**2))) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[biomassvars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
            biomasstime[:,yin,i] = (df[biomassvars[i]-3,1:,yin,xleft]*satiledg + df[biomassvars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
            biomasstime[:,yout,i] = (df[biomassvars[i]-3,1:,yout,xleft]*satoledg + df[biomassvars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
            biomasstime[:,yin+1:yout, i] = np.sum(df[biomassvars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[biomassvars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[biomassvars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
            bioconctime[:,yin,i] = ((df[biomassvars[i]-3,1:,yin,xleft]*satiledg + df[biomassvars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yin,xleft+1:xright-1]*satielem, axis=-1))*velem*vedge)/(vbc*vedge)
            bioconctime[:,yout,i] = ((df[biomassvars[i]-3,1:,yout,xleft]*satoledg + df[biomassvars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
            bioconctime[:,yin+1:yout, i] = (np.sum(df[biomassvars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[biomassvars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[biomassvars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge)/(vbc*velem)
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
       for i in range(len(biomassvars)):
           biomassendtime[i] = (df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[int(np.shape(df)[1])-2] + df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[int(np.shape(df)[1])-2] + df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[int(np.shape(df)[1])-2] + df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[int(np.shape(df)[1])-2])*(vedge**2) + sum(sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[int(np.shape(df)[1])-2,:,:]*(velem**2))) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]) + sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]) + sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yin,i] = (df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yout,i] = (df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomassendtimey[yin+1:yout, i] = sum(sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[biomassvars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
           biomasstime[:,yin,i] = (df[biomassvars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
           biomasstime[:,yout,i] = (df[biomassvars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
           biomasstime[:,yin+1:yout, i] = np.sum(df[biomassvars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[biomassvars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[biomassvars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
           bioconctime[:,yin,i] = ((df[biomassvars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + np.sum(df[biomassvars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1)*velem*vedge)/(vbc*vedge)
           bioconctime[:,yout,i] = ((df[biomassvars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[biomassvars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge)/(vbc*vedge)
           bioconctime[:,yin+1:yout, i] = (np.sum(df[biomassvars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[biomassvars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[biomassvars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge)/(vbc*velem)
    return df, biomassendtime, biomasstime, bioconctime

def heatmapconcdist (df, vars, Trial, gvarnames, d, fpre, title):
    #Heatmaps
    from matplotlib.ticker import FuncFormatter
    fmt = lambda x,pos: '{:4.0f}'.format(x)
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    titlesize = 20
    subtitlesize = 18
    vel = df[2,np.shape(df)[1]-1,:,:]*-1
    if ('DOC' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(8,8), sharex = True, sharey = True)
#        plt.tight_layout()
        plt.suptitle (title, ha='center', va='center', fontsize = titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt)})
            ax.flat[i].set_title(gvarnames[i], fontsize = subtitlesize)
        name = sns.heatmap(vel, square = False, cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[5],
                           cbar_kws={'format': FuncFormatter(fmtvel)})
        ax.flat[5].set_title("Velocity", fontsize = subtitlesize)
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_chem.png"
    elif ('Aerobes' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(7,9), sharex = True, sharey = True)
#        plt.tight_layout()
        plt.suptitle (title, ha='center', va='center', fontsize = titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(df[vars[i]-3,np.shape(df)[1]-1,:,:], square = False, cmap = "YlGnBu", xticklabels=10, yticklabels=10, ax=ax.flat[i],
                                  cbar_kws={'format': FuncFormatter(fmt)})
            ax.flat[i].set_title(gvarnames[i], fontsize = subtitlesize)
        name = sns.heatmap(vel, square = False, cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=ax.flat[3],
                           cbar_kws={'format': FuncFormatter(fmtvel)})
        ax.flat[3].set_title("Velocity", fontsize = subtitlesize)
        picname = d+fpre+str(Trial)+ "_" + "heatmap" + "_biomass.png"
    pictosave = name.get_figure()
    pictosave.savefig(picname, dpi = 300)    
    return (plt)

def heatmapconcdistLOG (df, vars, Trial, gvarnames, d, fpre, title):
    #Heatmaps
    from matplotlib.ticker import FuncFormatter
    fmt = lambda x,pos: '{:4.0f}'.format(x)
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    vel = df[2,np.shape(df)[1]-1,:,:]*-1
    if ('DOC' in gvarnames):
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize=(10,10), sharex = True, sharey = True)
        plt.suptitle (title)
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
        plt.suptitle (title)
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
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(inttime)):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[inttime[k],0:51,gvarnames.index(chem)], label = inttime[k]*50, color = colorstime[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
            ax[int(math.floor(col/ncols))][col%ncols].legend()
            ax[int(math.floor(col/ncols))][col%ncols].set_ylim(lylim, uylim)
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+ "_" + chem + "_" + datetime.now().strftime("%d%m%Y%H%M")+".png"
    plt.savefig(picname)
    return(plt)

def concprofileatss (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors):
    legendsize = 16
    axissize = 16
    ticksize = 14
    titlesize = 20
    for i in intsce:
        if ('DOC' in gvarnames):
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            #    print (np.mean(df[6, np.shape(df)[1]-1, yin, :]), np.mean(df[6, np.shape(df)[1]-1, yout-1, :]), np.mean(df[6, np.shape(df)[1]-1, yout, :]))
            yindex = list(range(51))
            fig, host = plt.subplots()
            fig.subplots_adjust(top=0.8)
        
            par1 = host.twiny()
            par2 = host.twiny()
            
            # Offset the top spine of par2.  The ticks and label have already been
            # placed on the top by twiny above.
            par2.spines["top"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["top"].set_visible(True)
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,0], yindex, label = gvarnames[0], color = colors[0])
            p2, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,1], yindex, label = gvarnames[1], color = colors[1])
            p3, = par1.plot(conctime[np.shape(conctime)[0]-1,0:51,2], yindex, label = gvarnames[2], color = colors[2])
            p4, = par2.plot(conctime[np.shape(conctime)[0]-1,0:51,3], yindex, label = gvarnames[3], color = colors[3])
            
            host.set_ylim(0,51)
            host.set_xlim(0, 800)
            par1.set_xlim(30, 60)
            par2.set_xlim(50, 260)
            
            plt.gca().invert_yaxis()
            
            host.set_ylabel("Y (cm)",fontsize = axissize)
            host.set_xlabel("DOC, DO (uM)",fontsize = axissize)
            par1.set_xlabel(str(gvarnames[2])+" (uM)",fontsize = axissize)
            par2.set_xlabel(str(gvarnames[3])+" (uM)",fontsize = axissize)
            
            host.xaxis.label.set_color(p1.get_color())
            par1.xaxis.label.set_color(p3.get_color())
            par2.xaxis.label.set_color(p4.get_color())
            
            tkw = dict(size=4, width=1.5, labelsize = ticksize)
            host.tick_params(axis='x', colors=p1.get_color(), **tkw)
            par1.tick_params(axis='x', colors=p3.get_color(), **tkw)
            par2.tick_params(axis='x', colors=p4.get_color(), **tkw)
            host.tick_params(axis='y', **tkw)
            
            lines = [p1, p2, p3, p4]
            
            host.legend(lines, [l.get_label() for l in lines], fontsize = legendsize, loc = "lower right")
            plt.title ("Variance "+str(Het[Trial.index(i)])+" : Anisotropy "+str(Anis[Trial.index(i)]), fontsize = titlesize, pad = 100)
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_chem_concatss.png"
        
        elif ('Aerobes' in gvarnames):
            df,massendtime, masstime, conctime = biomasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            yindex = list(range(51))
            fig, host = plt.subplots()
            fig.subplots_adjust(top=0.8)
            
            par1 = host.twiny()
            par2 = host.twiny()
            
            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["top"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["top"].set_visible(True)
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:51,0], yindex, label = gvarnames[0], color = colors[0])
            p2, = par1.plot(conctime[np.shape(conctime)[0]-1,0:51,1], yindex, label = gvarnames[1], color = colors[2])
            p3, = par2.plot(conctime[np.shape(conctime)[0]-1,0:51,2], yindex, label = gvarnames[2], color = colors[3])
            
            host.set_ylim(0, 51)
            host.set_xlim(0, 500)
            par1.set_xlim(0, 30)
            par2.set_xlim(0, 60)
            
            plt.gca().invert_yaxis()
            
            host.set_ylabel("Y (cm)",fontsize = axissize)
            host.set_xlabel(str(gvarnames[0])+" (uM)",fontsize = axissize)
            par1.set_xlabel(str(gvarnames[1])+" (uM)",fontsize = axissize)
            par2.set_xlabel(str(gvarnames[2])+" (uM)",fontsize = axissize)
            
            host.xaxis.label.set_color(p1.get_color())
            par1.xaxis.label.set_color(p2.get_color())
            par2.xaxis.label.set_color(p3.get_color())
            
            tkw = dict(size=4, width=1.5, labelsize = ticksize)
            host.tick_params(axis='x', colors=p1.get_color(), **tkw)
            par1.tick_params(axis='x', colors=p2.get_color(), **tkw)
            par2.tick_params(axis='x', colors=p3.get_color(), **tkw)
            host.tick_params(axis='y', **tkw)
            
            lines = [p1, p2, p3]
            
            host.legend(lines, [l.get_label() for l in lines], fontsize = legendsize, loc = "lower right")
            plt.title ("Variance "+str(Het[Trial.index(i)])+" : Anisotropy "+str(Anis[Trial.index(i)]), fontsize = titlesize, pad = 100)
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass_concatss.png"
        plt.savefig(picname, bbox_inches = 'tight')
    return None

def concprofilewithH (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars, gvarnames, intsce, colors,i):
    legendsize = 16
    axissize = 16
    ticksize = 14
    titlesize = 20
    if ('DOC' in gvarnames):
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars)
        dfh,massendtimeh, masstimeh, conctimeh, Velocityh, headh = calcconcmasstime (Trial[Trial.index('H')], Het[Trial.index('H')], Anis[Trial.index('H')], gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars)
        yindex = list(range(51))
        fig, host = plt.subplots()
        fig.subplots_adjust(top=0.8)
        
        par1 = host.twiny()
        par2 = host.twiny()
            
        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par2.spines["top"].set_visible(True)
    
        p1, = host.plot(conctimeh[-1,0:51,0], yindex, label = gvarnames[0], color = colors[0],
                                  linestyle = '-.')
        p1, = host.plot(conctime[-1,0:51,0], yindex, label = gvarnames[0], color = colors[0],
                                  linestyle = '-')
        p2, = host.plot(conctimeh[-1,0:51,1], yindex, label = gvarnames[1], color = colors[1],
                                  linestyle = '-.')
        p2, = host.plot(conctime[-1,0:51,1], yindex, label = gvarnames[1], color = colors[1],
                                  linestyle = '-')
        p3, = par1.plot(conctimeh[-1,0:51,2], yindex, label = gvarnames[2], color = colors[2],
                                  linestyle = '-.')
        p3, = par1.plot(conctime[-1,0:51,2], yindex, label = gvarnames[2], color = colors[2],
                                  linestyle = '-')
        p4, = par2.plot(conctimeh[-1,0:51,3], yindex, label = gvarnames[3], color = colors[3],
                              linestyle = '-.')
        p4, = par2.plot(conctime[-1,0:51,3], yindex, label = gvarnames[3], color = colors[3],
                                  linestyle = '-')
            
        host.set_ylim(0, 51)
        host.set_xlim(0, 800)
        par1.set_xlim(30, 60)
        par2.set_xlim(50, 260)
        
        host.set_ylabel("Y (cm)", fontsize = axissize)
        host.set_xlabel("DOC, DO (uM)", fontsize = axissize)
        par1.set_xlabel(str(gvarnames[2])+" (uM)", fontsize = axissize)
        par2.set_xlabel(str(gvarnames[3])+" (uM)", fontsize = axissize)
        
        host.xaxis.label.set_color(p1.get_color())
        par1.xaxis.label.set_color(p3.get_color())
        par2.xaxis.label.set_color(p4.get_color())
        
        tkw = dict(size=4, width=1.5, labelsize = ticksize)
        host.tick_params(axis='x', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='x', colors=p3.get_color(),**tkw)
        par2.tick_params(axis='x', colors=p4.get_color(),**tkw)
        host.tick_params(axis='y', **tkw)
        
        plt.gca().invert_yaxis()
        
        lines = [p1, p2, p3, p4]
        
        host.legend(lines, [l.get_label() for l in lines], fontsize = legendsize, loc = "lower right")
        plt.title ("Variance "+str(Het[Trial.index(i)])+" : Anisotropy "+str(Anis[Trial.index(i)]), fontsize = titlesize, pad = 100)
        
    elif ('Aerobes' in gvarnames):
        df,massendtime, ma, bioconctime = biomasstimefunc (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars)
        dfh,massendtimeh, mh, bioconctimeh = biomasstimefunc ('H', Het[Trial.index('H')], Anis[Trial.index('H')], gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars)
        yindex = list(range(51))
        fig, host = plt.subplots()
        fig.subplots_adjust(top=0.8)
        
        par1 = host.twiny()
        par2 = host.twiny()
        
        # Offset the right spine of par2.  The ticks and label have already been
        # placed on the right by twinx above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par2.spines["top"].set_visible(True)
        
        p1, = host.plot(bioconctimeh[-1,:,0], yindex, label = gvarnames[0], color = colors[0],
                                  linestyle = '-.')
        p1, = host.plot(bioconctime[-1,:,0], yindex, label = gvarnames[0], color = colors[0],
                                  linestyle = '-')
        p2, = par1.plot(bioconctimeh[-1,:,1], yindex, label = gvarnames[1], color = colors[2],
                                  linestyle = '-.')
        p2, = par1.plot(bioconctime[-1,:,1], yindex, label = gvarnames[1], color = colors[2],
                                  linestyle = '-')
        p3, = par2.plot(bioconctimeh[-1,:,2], yindex, label = gvarnames[2], color = colors[3],
                                  linestyle = '-.')
        p3, = par2.plot(bioconctime[-1,:,2], yindex, label = gvarnames[2], color = colors[3],
                                  linestyle = '-')

        
        host.set_ylim(0, 51)
        host.set_xlim(0, 500)
        par1.set_xlim(0, 30)
        par2.set_xlim(0, 200)
        
        lblw = dict(fontsize = axissize)
        host.set_ylabel("Y (cm)",**lblw)
        host.set_xlabel(str(gvarnames[0])+" (uM)", **lblw)
        par1.set_xlabel(str(gvarnames[1])+" (uM)", **lblw)
        par2.set_xlabel(str(gvarnames[2])+" (uM)", **lblw)
    
        host.xaxis.label.set_color(p1.get_color())
        par1.xaxis.label.set_color(p2.get_color())
        par2.xaxis.label.set_color(p3.get_color())
        
        tkw = dict(size=4, width=1.5, labelsize = ticksize)
        host.tick_params(axis='x', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='x', colors=p2.get_color(), **tkw)
        par2.tick_params(axis='x', colors=p3.get_color(), **tkw)
        host.tick_params(axis='y', **tkw)
        
        plt.gca().invert_yaxis()
        
        lines = [p1, p2, p3]
        
        host.legend(lines, [l.get_label() for l in lines], fontsize = legendsize, loc = "lower right")
        plt.title ("Variance "+str(Het[Trial.index(i)])+" : Anisotropy "+str(Anis[Trial.index(i)]), fontsize = titlesize, pad = 100)

    return plt


def concprofileatssX (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors):
    for i in intsce:
        if ('DOC' in gvarnames):
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstimeX (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
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
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:31,0], label = gvarnames[0], color = colors[0])
            p2, = host.plot(conctime[np.shape(conctime)[0]-1,0:31,1], label = gvarnames[1], color = colors[1])
            p3, = par1.plot(conctime[np.shape(conctime)[0]-1,0:31,2], label = gvarnames[2], color = colors[2])
            p4, = par2.plot(conctime[np.shape(conctime)[0]-1,0:31,3], label = gvarnames[3], color = colors[3])
            
            host.set_xlim(0, 31)
            host.set_ylim(0, 800)
            par1.set_ylim(30, 60)
            par2.set_ylim(50, 260)
            
            host.set_xlabel("X (cm)")
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
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_chem_concatss_X.png"
        
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
            
            p1, = host.plot(conctime[np.shape(conctime)[0]-1,0:31,0], label = gvarnames[0], color = colors[0])
            p2, = par1.plot(conctime[np.shape(conctime)[0]-1,0:31,1], label = gvarnames[1], color = colors[2])
            p3, = par2.plot(conctime[np.shape(conctime)[0]-1,0:31,2], label = gvarnames[2], color = colors[3])
            
            host.set_xlim(0, 31)
            host.set_ylim(0, 500)
            par1.set_ylim(0, 30)
            par2.set_ylim(0, 60)
            
            host.set_xlabel("X (cm)")
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
            picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass_concatss_X.png"
        plt.savefig(picname)
    return None

def plotconcallatss(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors, nrows, ncols, figsize):
    col=0
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(vars)-2):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[np.shape(conctime)[0]-1,0:51,k], label = gvarnames[k], color = colors[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (plt.title (str(Het[Trial.index(i)])+" : "+ str(Anis[Trial.index(i)])))
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+"_" +datetime.now().strftime("%d%m%Y%H%M")+"_concatss.png"
    plt.savefig(picname)
    
    return None

def plotconcwithhet (Regimes, intsce, chem, colorfamily, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15
    yindex = list(range(51))
    fig, axes = plt.subplots(nrows = 1, ncols = len(Regimes), figsize = [16,5])
    for r in Regimes:
        colors = sns.color_palette(colorfamily[Regimes.index(r)], len(intsce))
        d = r"X:/Saturated_flow/diffusion/" + r+ "AR_changedkindox/"
        lines = []
        for i in intsce:
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
            p1, = axes.flat[Regimes.index(r)].plot(conctime[np.shape(conctime)[0]-1,0:51,gvarnames.index(chem)], yindex,
                           label = str(Het[Trial.index(i)]) + ":" + str(Anis[Trial.index(i)]), color = colors[intsce.index(i)])
            lines.append(p1,)
            axes.flat[Regimes.index(r)].legend(lines, [l.get_label() for l in lines], fontsize = legendsize)

        if (r == "Equal"):
            axes.flat[Regimes.index(r)].set_title ("Medium flow", fontsize = titlesize)
        else:
            axes.flat[Regimes.index(r)].set_title (r+" flow", fontsize = titlesize)

        axes.flat[Regimes.index(r)].set_xlim(0, 250)
        axes.flat[Regimes.index(r)].set_ylabel("Y (cm)", fontsize = axissize)
        axes.flat[Regimes.index(r)].set_xlabel("DO (uM)", fontsize = axissize)
        axes.flat[Regimes.index(r)].tick_params(axis='y', labelsize = ticksize)
        axes.flat[Regimes.index(r)].tick_params(axis='x', labelsize = ticksize)
        axes.flat[Regimes.index(r)].invert_yaxis()  
    return fig

def plotoutflowtimeseries (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, inttime, colorstime, nrows, ncols, figsize, chem, lylim, uylim, start, end):
    col=0
    #figbig = plt.figure()
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[start:end,yout,gvarnames.index(chem)], label = "Tracer", color = "black")
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
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial2[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        for k in range(len(vars)-2):
            ax[int(math.floor(col/ncols))][col%ncols].plot(conctime[np.shape(conctime)[0]-1,0:51,k], label = gvarnames[k], color = colors[k])
            ax[int(math.floor(col/ncols))][col%ncols].set_title (Trial[Trial.index(i)])
        col = col+1
    plt.legend()
    picname = d+fpre+str(intsce[0])+"_"+str(intsce[-1])+"_" +datetime.now().strftime("%d%m%Y%H%M")+"_concatss.png"
    plt.savefig(picname)
    
    return None

def processchembiomassfiles(biomassfile, regime):
    biomass = pd.read_csv(biomassfile, delimiter = "\t")
    biomass['Chemb'] = biomass['Chem']
    biomass['VA']=biomass['Variance']*biomass['Anisotropy']
#    chemfile = pd.read_csv(chemfiles, delimiter = "\t")
#    oxy = pd.merge(biomass[biomass['Chemb']=='Aerobes'], chemfile[chemfile['Chem']=='DO'], on = 'Trial', suffixes=("_b","_c"))
#    amm = pd.merge(biomass[biomass['Chemb']=='Ammonia oxidizers'], chemfile[chemfile['Chem']=='Ammonium'], on = 'Trial', suffixes=("_b","_c"))
#    nitra = pd.merge(biomass[biomass['Chemb']=='Nitrate reducers'], chemfile[chemfile['Chem']=='Nitrate'], on = 'Trial', suffixes=("_b","_c"))
#    allbiomass = pd.concat([oxy,amm,nitra], axis=0, ignore_index = True)
#    biomass['delbiomass%']=biomass['TotalRatio']*100
#    biomass['del2biomass%']=biomass['TotalRatio_Time']*100
#    biomass['delspbiomass%']=biomass['SpeciesRatio']*100
#    biomass['del2spbiomass%']=biomass['SpeciesRatio_Time']*100
    biomass['Regime']=regime
    
    return biomass

def processchemfiles(chemfile, regime):
    chemfile = pd.read_csv(chemfile, delimiter = '\t')
    chemfile['Regime']=regime
    chemfile['VA']=chemfile['Variance']*chemfile['Anisotropy']
#    chemfile['RemovalRatio%']=chemfile['RemovalRatio']*100
#    chemfile['RemovalRatio_Time%']=chemfile['RemovalRatio_Time']*100    
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
            dum = sns.boxplot(x = "Xlabels", y = "fdelmassflux", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 400,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[-1]:
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
    di = r"Z:\Saturated_flow\Steady_state\Tracer_studies\/"
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
        if (breakthrough['Variance'][i] == 0.1):
            l.append(str(breakthrough['Variance'][i]) + ":" + str(breakthrough['Anisotropy'][i]))
        else:
            l.append(str(int(breakthrough['Variance'][i])) + ":" + str(breakthrough['Anisotropy'][i]))
        
    breakthrough['Xlabels'] = l
    breakthrough = breakthrough.sort_values(by=['Variance','Anisotropy'])
    
    combined_tracer = pd.read_csv("X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv", delimiter = "\t")
    combined_tracer.loc[combined_tracer.Regime=="Equal", 'Regime'] = "Medium"
    combined_tracer['fraction_withslow'] = combined_tracer["Time"]/combined_tracer["Time"][0] 
    l = []
    for i in range(len(combined_tracer)):
        if (combined_tracer['Variance'][i] == 0.1):
            l.append(str(combined_tracer['Variance'][i]) + ":" + str(combined_tracer['Anisotropy'][i]))
        else:
            l.append(str(int(combined_tracer['Variance'][i])) + ":" + str(combined_tracer['Anisotropy'][i]))
        
    combined_tracer['Xlabels'] = l
    combined_tracer = combined_tracer.sort_values(by=['Variance','Anisotropy'])

    return breakthrough, combined_tracer

def plot_tracer ():
    breakthrough, combined_tracer = tracerstudies()
    combined_tracer['%fraction'] = combined_tracer['fraction']*100
    combined_tracer['%fraction_withslow'] = combined_tracer['fraction_withslow']*100
    sns.set(rc = {'figure.figsize':(7,4)})
    sns.set_style("whitegrid")
    sns.boxplot(x = "Xlabels", y = "%fraction", hue = "Regime", data = combined_tracer,
                      hue_order = ["Slow", "Medium", "Fast"], palette = ["coral", "mediumseagreen", "steelblue"])
    plt.xlabel("Variance:Anisotropy")
    plt.ylabel("% of homogeneous scenario")
    plt.title("Time taken for tracer breakthrough")
    plt.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.png", dpi = 300)
    plt.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.pdf", dpi = 300)
    

    sns.set(rc = {'figure.figsize':(7,7)})
    sns.set_style("whitegrid")
    sns.boxplot(x = "Xlabels", y = "%fraction_withslow", hue = "Regime", data = combined_tracer,
                      hue_order = ["Slow", "Medium", "Fast"], palette = ["coral", "mediumseagreen", "steelblue"])
    plt.xlabel("Variance:Anisotropy")
    plt.ylabel("% of homogeneous scenario in slow flow regime")
    plt.yscale("log")
    plt.title("Time taken for tracer breakthrough")
    plt.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction_withslow.png", dpi = 300)
    plt.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction_withslow.pdf", dpi = 300)

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
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10])
    plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state", fontsize = 20)
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
#            pred = dfc["%ofhomogeneous"].to_numpy()
#            resp = dfc["del2massflux%"].to_numpy()
#            intercept , slope= np.polyfit(pred, resp,1)
            print (i, k, intercept, slope)
            dum = sns.lmplot(x = "%ofhomogeneous", y = "del2massflux%", data = dfall2, palette = colseries[Regimes.index(i)], col = "Regime", row = "Chem")
#            axes[colidx1][colidx2].scatter(pred, resp, c = colseries[Regimes.index(i)])
#            axes[colidx1][colidx2].plot(pred, slope + intercept*pred, '-', color = colseries[Regimes.index(i)])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            axes[colidx1][colidx2].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
            axes[colidx1][colidx2].set_xticklabels([])
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col+" flow", ha = 'center', fontsize = titlesize)
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 255,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    for ax in axes[-1]:
        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10), size = ticksize)
    plt.figtext(0.5, 0.08, 'Relative difference in breakthrough time (%)', ha='center', va='center', fontsize = 20)
    plt.figtext(0.08, 0.5, 'Relative difference (%)', ha='center', va='center', rotation='vertical', fontsize = 20)

    return dfall2, fig

def scatterrestime_biomass(dataset1, dataset2, dataset3):
    
    dfall = pd.concat([dataset1, dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    
    bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)
 
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10], sharex='col')
    plt.suptitle("Change in biomass with respect to homogeneous scenario at steady state", fontsize = 20)
    for i in Regimes:
    #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            pred = dfc["%ofhomogeneous"].to_numpy()
            resp = dfc["delbiomass%"].to_numpy()
            intercept , slope= np.polyfit(pred, resp,1)
            print (i, k, intercept, slope)
            axes[colidx1][colidx2].scatter(pred, resp, c = colseries[Regimes.index(i)], cmap = colseries[Regimes.index(i)])
            axes[colidx1][colidx2].plot(pred, slope + intercept*pred, '-', color = colseries[Regimes.index(i)])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            axes[colidx1][colidx2].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
            axes[colidx1][colidx2].set_xticklabels([])
#            axes[colidx1][colidx2].set_xticklabels([20, 0, -20, -40, -60, -80])
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
#    for ax in axes[:,0]:
#        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical',  fontsize = 15)
    for ax in axes[-1]:
#        ax.set_xlabel("Relative difference in breakthrough time (%)")
        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10))
    plt.figtext(0.5, 0.08, 'Relative difference in breakthrough time (%)', ha='center', va='center', fontsize = 20)
    plt.figtext(0.08, 0.5, 'Relative difference (%)', ha='center', va='center', rotation='vertical', fontsize = 20)   
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return dfall2, fig

def calcoxiccells(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    oxiccells= np.zeros([len(Trial), 51,1])
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        c = []
        for k in range(51):
            c.append(np.count_nonzero(df[vars[1]-3,np.shape(df)[1]-1,k,:]>limit))
#        print (Trial[j], c)
#        print ("Calculating number of oxic cells: ", c)
        oxiccells [j,:,0] = c
    
    return oxiccells

def plotoxiccellssdo(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars):
    axissize = 14
    ticksize = 12
    titlesize = 14
    oxiccells = calcoxiccells(limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    nrows = 2
    ncols = 4
    figsize = [14,8]
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True, sharey = True)
    count=0
    for j in Trial:
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (j, Het[Trial.index(j)], Anis[Trial.index(j)], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        axes.flat[count].plot(conctime[-1,:,1], 'r-')
        axes.flat[count].set_ylim(0,260)
        axes.flat[count].tick_params(axis='y', colors = 'r', labelsize = ticksize)
        axes.flat[count].set_title("Variance: "+str(Het[Trial.index(j)])+"& Anisotropy: "+str(Anis[Trial.index(j)]), fontsize = titlesize)
        ax2 = axes.flat[count].twinx()
        ax2.plot((oxiccells[Trial.index(j),:,0]/31)*100,'b-')
        ax2.set_ylim(0,105)
        ax2.tick_params(axis = 'y', colors = 'b', labelsize = ticksize)
        if (count+1) % ncols == 0:
            ax2.set_ylabel("% of oxic cells (-)", color='b', fontsize = axissize)
        else:
            ax2.set_yticklabels([]) #turning off secondary y axis yticklabels
        count = count+1
    for ax in axes[:,0]:
        ax.set_ylabel("DO (uM)", color='r', fontsize = axissize)
    for ax in axes[-1]:
        ax.set_xlabel("Y (cm)", fontsize = axissize)
    return fig

def boxV_Afluxtime(dataset1, dataset2, dataset3, chemseries, imgsize):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    Regimes = ["Slow", "Medium"]
    dfall["del2massflux_Time%"] = dfall["del2massflux_Time"]*100
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = imgsize)
    plt.suptitle("Change in removal of carbon and nitrogen in transient conditions", fontsize = titlesize)
    for i in Regimes:
        dfall0 = dfall[dfall['Time'] != 0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "del2massflux%", hue = "Time", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes.flat[Chems.index(k)].set_ylabel(k, fontsize = axissize)
            axes.flat[Chems.index(k)].legend_.remove()
            axes.flat[Chems.index(k)].tick_params(labelsize = ticksize)
            if (k != "Nitrate reducers"):
                axes.flat[Chems.index(k)].set_xlabel('')
                axes.flat[Chems.index(k)].set_xticklabels([])
            else:
                axes.flat[Chems.index(k)].set_xlabel("Variance : Anisotropy", fontsize = axissize)
            col = col+1
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col, fontsize = axissize)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)", fontsize = axissize)
    for ax, row in zip(axes[:,1], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 450,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[4]:
        ax.set_xlabel("Variance : Anisotropy", fontsize = axissize)
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def boxV_Abiotime(dataset1, dataset2, dataset3, imgsize):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True)
#    dfall = pd.concat([equalbc, slowbc, fastbc], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance_b'][i]) + ":" + str(dfall['Anisotropy_b'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance_b','Anisotropy_b'])
    dfall["delbiomass%"] = dfall["Change_umoles"]*100
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = imgsize)
    plt.suptitle("Change in total biomass at steady state")
    for i in Regimes:
        dfall0 = dfall[dfall['Time_b']!=0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "delbiomass%", hue = "Time_b", palette = colseries[Regimes.index(i)], data = dfc, ax = axes[colidx1][colidx2])
            axes[colidx1][colidx2].set_xlabel('')
            axes[colidx1][colidx2].set_ylabel('')
            col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:,0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,1], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 450,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center')
    for ax in axes[2]:
        ax.set_xlabel("Variance : Anisotropy")
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def scatterrestime_biomass_temp(dataset1, dataset2, dataset3, chemseries, plotYaxis, plotYaxislabel):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True, sort = False)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
#    dfall["delbiomass%"] = dfall["Change_umoles"]*100
#    dfall["delbiomass_Time%"] = dfall["Change_umoles_Time"]*100
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    
    bth1, bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'fraction']], on = ['Trial', 'Regime'])
    xlabels = np.arange(round(min(np.unique(dfall2['fraction'])),1), round(max(np.unique(dfall2['fraction'])),1)+0.2, 0.2)
#    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
#    markerstyles = ["o", "d","^", "s"]
    ncols = len(Regimes)
    nrows = len(Chems)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobilized", "Mobile"]
    activity = ["Active", "Inactive"]
    ncols = len(species)
    nrows = len(position)+len(activity)
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = [len(Chems)*1.5,len(Chems)*1.5])
    #plt.suptitle("Change in total biomass with respect to homogeneous scenario at steady state", fontsize = 20)
    for k in Chems:
        dfc = dfall2[dfall2['Chem']==k]
        colidx1 = Chems.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp['Regime']==i]
            print (i)
            axes.flat[colidx1].scatter("fraction", plotYaxis,  c = "Time", cmap = colseries[Regimes.index(i)], data = dfcr, label = "Time")
            #        axes.flat[colidx1].plot(pred, slope + intercept*pred, '--')
            #        axes.flat[colidx1].set_ylabel(k, fontsize = 15)
            axes.flat[colidx1].tick_params(axis = 'y', labelsize = 15)
            #        axes.flat[colidx1].set_xscale('log')
            if(Chems.index(k) < len(Chems)-3):
                axes.flat[colidx1].set_xlabel('')
                axes.flat[colidx1].set_xticklabels([])
            else:
                axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
#            axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case', fontsize = 15)
    plt.legend(loc = "best", fontsize = 12)
    for ax,typsp in zip(axes[0,:], species):
        ax.set_title(typsp, fontsize = 15)
    axes[0,-1].annotate(position[0], xy = (0, 0.5), xytext = (300,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    axes[1,-1].annotate(position[0], xy = (0, 0.5), xytext = (300,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    axes[2,-1].annotate(position[1], xy = (0, 0.5), xytext = (300,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    axes[3,-1].annotate(position[1], xy = (0, 0.5), xytext = (300,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    axes[0,0].set_ylabel(activity[0], fontsize = 15)
    axes[1,0].set_ylabel(activity[1], fontsize = 15)
    axes[2,0].set_ylabel(activity[0], fontsize = 15)
    axes[3,0].set_ylabel(activity[1], fontsize = 15)
    plt.annotate(plotYaxislabel, xy = (0, 2.2), xytext = (-800,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    plt.annotate("Fraction of breakthrough time in base case", xy = (-0.4, -0.3), xytext = (-100,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'center', va = 'center', fontsize = 15)

#    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10]), sharex='col')
#    plt.suptitle("Change in total biomass in transient conditions", fontsize = titlesize)
#    for i in Regimes:
#        dfall0 = dfall2[dfall2['Time_b']!=0]
#        df = dfall2[dfall2['Regime']==i]
#        col = 0
##        for t in range(4):
#            dft = df[df['Time']==t]
#            m = markerstyles[t]
#            for k in Chems:
#                dfc = df[df['Chemb']==k]
#                colidx1 = Chems.index(k)
##                colidx2 = Regimes.index(i)
#                axes[colidx1][colidx2].scatter("fraction", plotYaxis, c = "Time", cmap = colseries[Regimes.index(i)], data = dfc)
#                    facecolors = 'none', edgecolors = colseries[Regimes.index(i)])
#                axes[colidx1][colidx2].set_xlabel('')
#                axes[colidx1][colidx2].set_ylabel('')
#                axes[colidx1][colidx2].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
#                axes[colidx1][colidx2].set_xticklabels([])
#                axes[colidx1][colidx2].set_xticklabels([20, 0, -20, -40, -60, -80])
##                col = col+1
#    for ax, col in zip(axes[0], Regimes):
#        ax.set_title(col, fontsize = axissize)
##    for ax in axes[:,0]:
##        ax.set_ylabel("Relative difference (%)")
#    for ax, row in zip(axes[:,2], Chems):
#        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical',  fontsize = 15)
#    for ax in axes[-1]:
##        ax.set_xlabel("Relative difference in breakthrough time (%)")
#        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10), size = ticksize)
#    plt.figtext(0.5, 0.08, 'Relative difference in breakthrough time (%)', ha='center', va='center', fontsize = axissize)
#    plt.figtext(0.08, 0.5, plotYaxislabel, ha='center', va='center', rotation='vertical', fontsize = axissize)   
#    fig.subplots_adjust(left=0.15, top=0.9)
#    
    return fig

def scatterrestime_flux_temp_singleaxis(dataset1, dataset2, dataset3, chemseries, plotYaxis, plotYaxislabel):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    bth1, bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'fraction']], on = ['Trial', 'Regime'])
    
    Regimes = ["Slow", "Medium", "Fast"]
    
#    dfall2["del2massflux_Time%"] = dfall2["del2massflux_Time"]*100
    Chems = chemseries
#    colseries = ["indianred", "g", "steelblue"]
    colseries = ["Reds", "Greens", "Blues"]
    markerstyles = ["o", "d", "^", "s"]
    ncols = len(Regimes)
    nrows = len(Chems)
    xlabels = np.arange(round(min(np.unique(dfall2['fraction'])),1), round(max(np.unique(dfall2['fraction'])),1)+0.2, 0.2)
    print (xlabels)
    fig, axes = plt.subplots(ncols = 1, nrows = nrows, figsize = [8,5])#, sharex = True)
    plt.suptitle("Change in removal of carbon and nitrogen in transient conditions", fontsize = titlesize)
    for i in Regimes:
#        dfall0 = dfall2[dfall2['Time_b']!=0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for t in range(4):
            dft = df[df['Time']==t]
            m = markerstyles[t]
            for k in Chems:
                dfc = df[df['Chem']==k]
                colidx1 = Chems.index(k)
                colidx2 = Regimes.index(i)
                axes[colidx1].scatter("fraction", plotYaxis, c = "Time", cmap =colseries[Regimes.index(i)], data = dfc)
                #    facecolors = 'none', edgecolors = colseries[Regimes.index(i)])
                axes[colidx1].set_xlabel('')
                axes[colidx1].set_ylabel('')
                axes[colidx1].set_xticks(xlabels)
                if (k != Chems[-1]):
                    axes[colidx1].set_xticklabels([])
                else:
                    axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
                    axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case', fontsize = 15)
#    plt.legend(loc = "best", fontsize = 12)
    plt.annotate(plotYaxislabel, xy = (0, 2.2), xytext = (- 80,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)

    return fig

def scatterrestime_flux_temp(dataset1, dataset2, dataset3, chemseries, plotYaxis, plotYaxislabel):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15 
    dfall = pd.concat([dataset2, dataset3], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    bth1, bth = tracerstudies()
    
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'fraction']], on = ['Trial', 'Regime'])
    
    Regimes = ["Slow", "Medium", "Fast"]
    
#    dfall2["del2massflux_Time%"] = dfall2["del2massflux_Time"]*100
    Chems = chemseries
    toplot = str(plotYaxis)
    toplotlabel = str(plotYaxislabel)
#    colseries = ["indianred", "g", "steelblue"]
    colseries = ["Reds", "Greens", "Blues"]
    markerstyles = ["o", "d", "^", "s"]
    ncols = len(Regimes)
    nrows = len(Chems)
    xlabels = np.arange(round(min(np.unique(dfall2['fraction'])),1), round(max(np.unique(dfall2['fraction'])),1)+0.2, 0.2)
    print (xlabels)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10])#, sharex = True)
    plt.suptitle("Change in removal of carbon and nitrogen in transient conditions", fontsize = titlesize)
    for i in Regimes:
#        dfall0 = dfall2[dfall2['Time_b']!=0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for t in range(4):
            dft = df[df['Time']==t]
            m = markerstyles[t]
            for k in Chems:
                dfc = df[df['Chem']==k]
                colidx1 = Chems.index(k)
                colidx2 = Regimes.index(i)
                axes[colidx1][colidx2].scatter("fraction", toplot, c = "Time", cmap =colseries[Regimes.index(i)], data = dfc)
                #    facecolors = 'none', edgecolors = colseries[Regimes.index(i)])
                axes[colidx1][colidx2].set_xlabel('')
                axes[colidx1][colidx2].set_ylabel('')
                axes[colidx1][colidx2].set_xticks(xlabels)
                if (k != Chems[-1]):
#                    axes[colidx1][colidx2].set_xticklabels(xlabels, size = ticksize)    
#                else:
                    axes[colidx1][colidx2].set_xticklabels([])
                col = col+1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col, fontsize = axissize)
#    for ax in axes[:,0]:
#        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:,2], Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 250,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical',  fontsize = 15)
#    for ax in axes[-1]:
#        ax.set_xticklabels(xlabels, size = ticksize)
    plt.figtext(0.5, 0.08, 'Ratio of breakthrough time with the base case', ha='center', va='center', fontsize = axissize)
    plt.figtext(0.08, 0.5, toplotlabel, ha='center', va='center', rotation='vertical', fontsize = axissize)   
    fig.subplots_adjust(left=0.15, top=0.9)
    
    return fig

def calcsum_temp(Trial, Het, Anis, gw, newd, fpre, fsuf, yin, yout, xleft, xright, biomassvars, biomassgvarnames):
    yout = 50
    yin = 0
    xleft = 0
    xright = 30
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    df = np.load(newd+fpre+str(Trial)+fsuf+fpre+str(Trial)+"_df.npy")
    sumalltime = np.zeros([np.shape(df)[1]-1,len(biomassvars)])
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
    for b in range(len(biomassgvarnames)):
        sumalltime[:, b] = ((df[biomassvars[b]-3,1:,yin,xleft]*satiledg + df[biomassvars[b]-3,1:,yout,xleft]*satoledg + df[biomassvars[b]-3,1:,yin,xright]*satiredg + df[biomassvars[b]-3,1:,yout,xright]*satoredg)*vedge*vedge + (np.sum(df[biomassvars[b]-3,1:,yin,xleft+1:xright]*satielem, axis = -1) + np.sum(df[biomassvars[b]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1) + np.sum(df[biomassvars[b]-3,1:,yin+1:yout,xleft]*satlelem, axis = -1) + np.sum(df[biomassvars[b]-3,1:,yin+1:yout,xright]*satrelem, axis = -1))*vedge*velem) + np.sum(np.sum(df[biomassvars[b]-3,1:,yin+1:yout,xleft+1:xright]*satelem, axis = -1), axis = -1)*velem*velem

    return sumalltime

def concprofilewithtime_temp (Trial, intsce, d, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars):
    titlesize = 20
    legendsize = 12
    ticksize = 15
    axissize = 15
    xindex = list(range(0,1095*5,5))
    if ("DOC" in gvarnames):
        figbig, axes = plt.subplots(nrows=5, ncols=len(intsce), figsize=[30, 12], sharex = True)
        plt.suptitle("Concentration profile with time", fontsize = titlesize)
        for i in intsce:
            newd = d + j
            Titles = "Variance:"+str(Het[Trial.index(i)])+" & Anisotropy:"+str(Anis[Trial.index(i)])
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
            axes[0][intsce.index(i)].plot(xindex, head, label = "Velocity (m/d)", color = "black")
            axes[0][intsce.index(i)].set_title(Titles, fontsize = titlesize)
            axes[1][intsce.index(i)].plot(xindex, conctime[1:,yout,0], label = gvarnames[0], color = "black")
            axes[2][intsce.index(i)].plot(xindex, conctime[1:,yout,1], label = gvarnames[1], color = "red")
            axes[3][intsce.index(i)].plot(xindex, conctime[1:,yout,2], label = gvarnames[2], color = "blue")
            axes[4][intsce.index(i)].plot(xindex, conctime[1:,yout,3], label = gvarnames[3], color = "green")
            
            ax = axes[:,intsce.index(i)]
            if (intsce.index(i)==0):
                off = 75
                ax[0].annotate("Velocity", xy = (0,0.5), xytext = (-ax[0].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[1].annotate(str(gvarnames[0]), xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[2].annotate(str(gvarnames[1]), xy = (0,0.5), xytext = (-ax[2].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[3].annotate(str(gvarnames[2]), xy = (0,0.5), xytext = (-ax[3].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[4].annotate(str(gvarnames[3]), xy = (0,0.5), xytext = (-ax[4].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
            else:
                ax[0].set_yticklabels([])
                ax[1].set_yticklabels([])
                ax[2].set_yticklabels([])
                ax[3].set_yticklabels([])
                ax[4].set_yticklabels([])

        for ax in axes[:,0]:
            ax.tick_params(labelsize = ticksize)
        for ax in axes[4]:
            ax.set_xlabel("Time (days)", fontsize = axissize)      
            ax.tick_params(labelsize = ticksize)
    elif ("Aerobes" in gvarnames):
        figbig, axes = plt.subplots(nrows=4, ncols=len(intsce), figsize=[30, 12], sharex = True)
        plt.suptitle("Concentration profile with time", fontsize = titlesize)
        for i in intsce:
            newd = d + j
            Titles = "Variance:"+str(Het[Trial.index(i)])+" & Anisotropy:"+str(Anis[Trial.index(i)])
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
            conctime = calcsum_temp(Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
            axes[0][intsce.index(i)].plot(xindex, head, label = "Velocity (m/d)", color = "black")
            axes[0][intsce.index(i)].set_title(Titles, fontsize = titlesize)
            axes[1][intsce.index(i)].plot(xindex, conctime[:,0], label = gvarnames[0], color = "black")
            axes[2][intsce.index(i)].plot(xindex, conctime[:,1], label = gvarnames[1], color = "red")
            axes[3][intsce.index(i)].plot(xindex, conctime[:,2], label = gvarnames[2], color = "blue")
            
            ax = axes[:,intsce.index(i)]
            if (intsce.index(i)==0):
                off = 75
                ax[0].annotate("Velocity", xy = (0,0.5), xytext = (-ax[0].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[1].annotate(str(gvarnames[0]), xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[2].annotate(str(gvarnames[1]), xy = (0,0.5), xytext = (-ax[2].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[3].annotate(str(gvarnames[2]), xy = (0,0.5), xytext = (-ax[3].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
            else:
                ax[0].set_yticklabels([])
                ax[1].set_yticklabels([])
                ax[2].set_yticklabels([])
                ax[3].set_yticklabels([])

        for ax in axes[:,0]:
            ax.tick_params(labelsize = ticksize)
        for ax in axes[3]:
            ax.set_xlabel("Time (days)", fontsize = axissize)      
            ax.tick_params(labelsize = ticksize)
    return figbig

def shiftoxic_anoxic_temp (chem, Trial, intsce, d, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars):
    titlesize = 20
    legendsize = 12
    ticksize = 15
    axissize = 15
    yindex = list(range(51))
    nrows = 2
    ncols = 3
    figsize = [14,8]
    colors = ['black','red', 'blue', 'green', 'orange', 'grey']
    times = [2*5, 180*5, 360*5, 500*5, 1000*5]
    #Plot all the scenarios excluding steady state scenarios
#    intsce = [76, 73, 80, 84, 44, 63]
    col=0
    figbig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    plt.suptitle(chem + " profile at different time points", fontsize = titlesize)
    for i in intsce:
        newd = d + j 
        df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
        for c in times:
            Titles = "Variance:"+str(Het[Trial.index(i)])+" & Anisotropy:"+str(Anis[Trial.index(i)])
            axes[int(math.floor(col/ncols))][col%ncols].plot(conctime[int(np.shape(conctime)[0]-c/5),:,gvarnames.index(chem)], yindex, label = int(np.shape(conctime)[0]-c/5)*5, color = colors[times.index(c)])
            axes[int(math.floor(col/ncols))][col%ncols].set_title(Titles, fontsize = axissize)
            axes[int(math.floor(col/ncols))][col%ncols].tick_params(labelsize = ticksize)
        axes[int(math.floor(col/ncols))][col%ncols].invert_yaxis()
        col = col+1  
    axes.flat[5].legend((int(5*(np.shape(conctime)[0]-times[0]/5)), int(5*(np.shape(conctime)[0]-times[1]/5)),
                int(5*(np.shape(conctime)[0]-times[2]/5)), int(5*(np.shape(conctime)[0]-times[3]/5)),
                int(5*(np.shape(conctime)[0]-times[4]/5))), fontsize = legendsize)
    for ax in axes[1]:
        ax.set_xlabel("Concentration of "+ chem + " (uM)", size = axissize)
    for ax in axes[:,0]:
        ax.set_ylabel("Y (cm)", size = axissize) 
    
    return figbig

def localmaxmin(y1, y2, init):
    ymax = np.max(np.max(np.append(y1,y2),axis = 0))
    ymin = np.min(np.min(np.append(y1,y2), axis = 0))
    xposmax = (np.argmax(y1) + init)*5
    xposmin = (np.argmin(y1) + init)*5
    print(ymax, ymin, xposmax, xposmin)
    return ymax, ymin, xposmax, xposmin

def generate_timeseries(Regimes, initseries, lastseries, Trial, Het, Anis, gw, d, fpre, vars, gvarnames, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames):
    for Reg in Regimes:
        titlesize = 35
        legendsize = 25
        ticksize = 30
        axissize = 25
        line = '--'
        xindex = list(range(0,1095*5,5))
        Tforfpre = ['0/', '1/', '2/', '5/']
        df0,massendtime0, masstime0, conctime0, Velocity0, head0 = calcconcmasstime (Trial, Het, Anis, gw, d+Reg + "AR_"+Tforfpre[0], fpre, fsuf, yin, yout, xleft, xright, vars)
        sumalltime0 = calcsum_temp (Trial, Het, Anis, gw, d+Reg + "AR_"+Tforfpre[0], fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
        h = []
        c = np.zeros([1095,np.shape(conctime0)[2]])
        s = np.zeros([1095,3])
        for idx in range(1095):
            h.append(head0[-1])
            s[idx,:] = sumalltime0[-1,:]
            c[idx,:] = conctime0[-1,yout,:]
        figbig, axes = plt.subplots(nrows=8, ncols=len(Tforfpre)-1, figsize=[30, 20])
        plt.suptitle("Variance: "+ str(Het)+" & Anisotropy: "+str(Anis), fontsize = titlesize)
        for j, init, last in zip(Tforfpre[1:], initseries, lastseries):
            newd = d + Reg + "AR_" + j
            Titles = "Time series: "+str(Tforfpre.index(j))
            df,massendtime, masstime, conctime, Velocity, head = calcconcmasstime (Trial, Het, Anis, gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
            sumalltime = calcsum_temp(Trial, Het, Anis, gw, newd, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
            axes[0][Tforfpre.index(j)-1].plot(xindex[init:last], head[init:last], label = "Groundwater head (m)", color = "black")
            axes[0][Tforfpre.index(j)-1].plot(xindex[init:last], h[init:last], label = "Groundwater head (m)", color = "black", linestyle = line)        
            ymax, ymin, xposmax, xposmin = localmaxmin(head[init:last],h[init:last], init)
            axes[0][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[0][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[0][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax - (ymax-ymin)/2), size = legendsize)
            axes[0][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin + (ymax-ymin)/2), size = legendsize)
            axes[0][Tforfpre.index(j)-1].set_title(Titles, fontsize = titlesize)
            axes[1][Tforfpre.index(j)-1].plot(xindex[init:last], conctime[init+1:last+1,yout,0], label = gvarnames[0], color = "black")
            axes[1][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,1], label = gvarnames[0], color = "black", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,0], c[init:last,1], init)
            axes[1][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[1][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[1][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[1][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[3][Tforfpre.index(j)-1].plot(xindex[init:last], conctime[init+1:last+1,yout,1], label = gvarnames[1], color = "red")
            axes[3][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,2], label = gvarnames[1], color = "red", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,1], c[init:last,2], init)
            axes[3][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[3][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[3][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[3][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[5][Tforfpre.index(j)-1].plot(xindex[init:last], conctime[init+1:last+1,yout,2], label = gvarnames[2], color = "blue")
            axes[5][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,3], label = gvarnames[2], color = "blue", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,2], c[init:last,3], init)
            axes[5][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[5][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[5][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[5][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)        
            axes[7][Tforfpre.index(j)-1].plot(xindex[init:last], conctime[init+1:last+1:,yout,3], label = gvarnames[3], color = "green")
            axes[7][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,4], label = gvarnames[3], color = "green", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,3], c[init:last,4], init)
            axes[7][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[7][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[7][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[7][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[2][Tforfpre.index(j)-1].plot(xindex[init:last], sumalltime[init:last,0], label = AFbiomassgvarnames[0], color = "red")
            axes[2][Tforfpre.index(j)-1].plot(xindex[init:last], s[init:last,0], label = AFbiomassgvarnames[0], color = "red", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(sumalltime[init:last,0], s[init:last,0], init)
            axes[2][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[2][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[2][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[2][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[4][Tforfpre.index(j)-1].plot(xindex[init:last], sumalltime[init:last,1], label = AFbiomassgvarnames[1], color = "blue")
            axes[4][Tforfpre.index(j)-1].plot(xindex[init:last], s[init:last,1], label = AFbiomassgvarnames[1], color = "blue", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(sumalltime[init:last,1], s[init:last,1], init)
            axes[4][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[4][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[4][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[4][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[6][Tforfpre.index(j)-1].plot(xindex[init:last], sumalltime[init:last,2], label = AFbiomassgvarnames[2], color = "green")
            axes[6][Tforfpre.index(j)-1].plot(xindex[init:last], s[init:last,2], label = AFbiomassgvarnames[2], color = "green", linestyle = line)
            ymax, ymin, xposmax, xposmin = localmaxmin(sumalltime[init:last,2], s[init:last,2], init)
            axes[6][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            axes[6][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            axes[6][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            axes[6][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[0][Tforfpre.index(j)-1].set_xticklabels([])
            axes[1][Tforfpre.index(j)-1].set_xticklabels([])
            axes[2][Tforfpre.index(j)-1].set_xticklabels([])
            axes[3][Tforfpre.index(j)-1].set_xticklabels([])
            axes[4][Tforfpre.index(j)-1].set_xticklabels([])
            axes[5][Tforfpre.index(j)-1].set_xticklabels([])
            axes[6][Tforfpre.index(j)-1].set_xticklabels([])
            ax = axes[:,Tforfpre.index(j)-1]
            if (Tforfpre.index(j)-1==0):
                off = 150
                ax[0].annotate("Velocity", xy = (0,0.5), xytext = (-ax[0].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[1].annotate(str(gvarnames[0]), xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[3].annotate(str(gvarnames[1]), xy = (0,0.5), xytext = (-ax[3].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[5].annotate(str(gvarnames[2]), xy = (0,0.5), xytext = (-ax[5].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[7].annotate(str(gvarnames[3]), xy = (0,0.5), xytext = (-ax[7].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[2].annotate(str(AFbiomassgvarnames[0]), xy = (0,0.5), xytext = (-ax[2].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[4].annotate('Ammonia\noxidizers', xy = (0,0.5), xytext = (-ax[4].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
                ax[6].annotate('Nitrate\nreducers', xy = (0,0.5), xytext = (-ax[6].yaxis.labelpad-off,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = axissize)
            else:
                ax[0].set_yticklabels([])
                ax[1].set_yticklabels([])
                ax[2].set_yticklabels([])
                ax[3].set_yticklabels([])
                ax[4].set_yticklabels([])
                ax[5].set_yticklabels([])
                ax[6].set_yticklabels([])
                ax[7].set_yticklabels([])

        for ax in axes[:,0]:
            ax.tick_params(labelsize = ticksize)
        for ax in axes[7]:
            ax.set_xlabel("Time (days)", fontsize = axissize)      
            ax.tick_params(labelsize = ticksize)
        picname = "Z:/Saturated_flow/diffusion_transient/chemsamndbiomasswithvel_temp_"+str(Tforfpre.index(j))+"_"+Reg+"_"+str(Trial)+"_ZOOMED.png"
        plt.savefig(picname, dpi = 300, bbox_inches='tight', pad_inches = 0)
    
    return None

def biomassvscarryingcapacity(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars):
    maxcap = 500
    yindex = list(range(51))
    cc = np.zeros([51,1])
    df,massendtime, ma, bioconctime = biomasstimefunc (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars)    
    c = ["firebrick", "darkblue", "darkgreen", "lightcoral", "cornflowerblue", "mediumseagreen"]
    l = ["Active aerobes", "Active ammonia oxidizers", "Active nitrate reducers",
         "Inactive aerobes", "Inactive ammonia oxidizers", "Inactive nitrate reducers"]
    bioconctimet = bioconctime[-1, :, :].transpose()
    cc[:,0] = maxcap
    plt.figure()       
    plt.stackplot(yindex, bioconctimet, colors = c, labels = l)
    plt.xlabel("Y (cm)")
    plt.ylabel("Average concentration (uM C)")
#    plt.plot(yindex, cc, label = "Carrying capacity")
    plt.xlim(0, 50)
#    plt.xlim(300,510)
    plt.title ("Variance: "+str(Het)+" & Anisotropy: "+str(Anis))
#    plt.gca().invert_yaxis()
    plt.legend(loc = "lower right")        
    
    return plt
    
def removalperbiomass():
    removal = pd.read_csv("X:/massflux_withbreakthrough_revisedwithmobilebiomass.csv", delimiter="\t")
    Biomass = pd.read_csv("X:/biomass_withbreakthrough.csv", delimiter = "\t")
    newseries = removal[removal["Chem"]=="DO"]
    newseries['Rem'] = 100*(newseries.loc[:,"Inlet_massflux"] - newseries.loc[:,"Outlet_massflux"])/newseries.loc[:,"Inlet_massflux"]
    newb1series = Biomass[Biomass["Chem"]=="Active fixed Aerobes"]
    newb2series = Biomass[Biomass["Chem"]=="Active mobile Aerobes"]
    newbseries = newb1series.loc[:,'Total_biomass'].reset_index().add(newb2series.loc[:,'Total_biomass'].reset_index())
    num = newseries.loc[:,'Rem'].reset_index().to_numpy
    b = newbseries.loc[:,'Total_biomass'].reset_index().to_numpy
    for i in range(len(num)):
        num[:,'Rem'][i]/b[i]
    newseries.loc[:,'Rem'].reset_index().div(newbseries.loc[:,'Total_biomass'].reset_index())
    
    num.div(newbseries.loc[:,'Total_biomass'].reset_index())
    for i in range(len(newseries)):
        newseries.loc[:,'Rem'].reset_index().div(newbseries.loc[:,'Total_biomass'].reset_index())
    return mf