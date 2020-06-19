# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:56:16 2020

@author: khurana
"""

import numpy  as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import signal

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
    mft = np.zeros([(len(Trial)*len(Tforfpre))*len(gvarnames),12])
    #dire = d + fpre + str(Trialhead)
    sumall = np.zeros([(len(Trial)*len(Tforfpre))*len(biomassgvarnames),14])
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

def calcconcmasstime (Trial, Het, Anis, gw, directory, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
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
    di =  directory+fpre+str(Trial)+fsuf
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
#                masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
#                masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
#                masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
                conctime[1:,yin,i] = ((df[vars[i]-3,1:,yin,xleft]*satiledg*veliledg + df[vars[i]-3,1:,yin,xright]*satiredg*veliredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem*velielem, axis=-1))*velem)/(vedge*(veliredg+veliledg) + np.sum(velem*velielem, axis = -1))
                conctime[1:,yout,i] = ((df[vars[i]-3,1:,yout,xleft]*satoledg*veloledg + df[vars[i]-3,1:,yout,xright]*satoredg*veloredg)*(vedge) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem*veloelem, axis = -1))*velem)/(vedge*(veloredg+veloledg) + np.sum(velem*veloelem, axis = -1))
                conctime[1:,yin+1:yout, i] = (np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*velem*velelem,axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem*vellelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem*velrelem)*vedge)/(vedge*(vellelem+velrelem) + np.sum(velem*velelem, axis = -1))
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
#           massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[int(np.shape(df)[1])-2])*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[int(np.shape(df)[1])-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
#           massendtimey[yin,i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
#           massendtimey[yout,i] = (df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
#           massendtimey[yin+1:yout, i] = sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
#           masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
#           masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
#           masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
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

def autocorrelation (directory, timevari, triali, Trial, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, AFbiomassvars, AFbiomassgvarnames):
    df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(triali, Het[Trial.index(triali)], Anis[Trial.index(triali)], gw, directory+timevari, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
    sumall0 = calcsum_temp(triali, Het[Trial.index(triali)], Anis[Trial.index(triali)], gw, directory+timevari, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
    normavgconcout = np.zeros([np.shape(df0)[1],len(gvarnames)])
    normsum = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
    spratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
    normspratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
    corrchem = np.zeros([2*np.shape(Headinlettime0)[0],(len(gvarnames))])
    corrbiomass = np.zeros([2*np.shape(Headinlettime0)[0]-1,(len(AFbiomassgvarnames))])
    for k in range(len(gvarnames)):
        normavgconcout[:,k] = (conctime0[:,yout,k] - np.mean(conctime0[:,yout,k]))
        normheadin = Headinlettime0 - np.mean(Headinlettime0)
        for k in range(len(gvarnames)):
            corrchem[:,k] = np.correlate(normavgconcout[:,k],normheadin,'full')/((np.std(conctime0[:,yout,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
        for k in range(len(AFbiomassgvarnames)):
            normsum[:,k] = (sumall0[:,k] - np.mean(sumall0[:,k]))
            spratio[:,k] = sumall0[:,k]/np.sum(sumall0[:,:], axis = -1)
        for k in range(len(AFbiomassgvarnames)):
            normspratio[:,k] = spratio[:,k] - np.mean(spratio[:,k])
        for k in range(len(AFbiomassgvarnames)):
            corrbiomass[:,k] = np.correlate(normspratio[:,k],normheadin,'full')/((np.std(spratio[:,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
    
    return corrchem, corrbiomass, Headinlettime0


def correlationanalysis (directory, numberofTempscenarios, Trial, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, AFbiomassvars, AFbiomassgvarnames):
    Regimes = ["Slow","Equal", "Fast"]
    amplitudec = np.zeros([len(Trial)*(numberofTempscenarios)*len(Regimes)*(len(gvarnames)),8])
    amplitudeb = np.zeros([len(Trial)*(numberofTempscenarios)*len(Regimes)*(len(AFbiomassgvarnames)),8])
    conccorr = np.zeros([len(Trial)*(numberofTempscenarios)*len(Regimes)*(len(gvarnames)),8])
    conccorrb = np.zeros([len(Trial)*(numberofTempscenarios)*len(Regimes)*(len(AFbiomassgvarnames)),8])
    plotg = ["DOC","DO","TOC","Nitrogen"]
    for Reg in Regimes:
        print (Reg)
        Tforfpre = [Reg+'AR_0', Reg+'AR_1', Reg+'AR_2',Reg+'AR_5']
        for p in Tforfpre[1:]:
            print (p)
            for t in Trial:
                print (t)
#                df = np.load(directory+p+fpre+str(t)+fpre+str(t)+"_df.npy")
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
                sumall0 = calcsum_temp(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
                normavgconcout = np.zeros([np.shape(df0)[1],len(gvarnames)])
                normsum = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                spratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                normspratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                for k in range(len(gvarnames)):
                    normavgconcout[:,k] = (conctime0[:,yout,k] - np.mean(conctime0[:,yout,k]))
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [11,8], sharex = True)
                plt.suptitle(Reg + p + str(t))
                for k in range(len(gvarnames)):
                    idxc = Regimes.index(Reg)*(len(Tforfpre)-1)*len(Trial)*(len(gvarnames)) + (Tforfpre.index(p)-1)*len(Trial)*(len(gvarnames)) + Trial.index(t)*(len(gvarnames))+k
                    if (t=='H'):
                        conccorr[idxc,0] = Trial.index(t)
                    else:
                        conccorr[idxc,0] = t
                    conccorr[idxc,1] = Het[Trial.index(t)]
                    conccorr[idxc,2] = Anis[Trial.index(t)]
                    conccorr[idxc,3] = Tforfpre.index(p)
                    conccorr[idxc,4] = Regimes.index(Reg)
                    conccorr[idxc,5] = k
                    conccorrarray = np.correlate(normavgconcout[:,k],normheadin,'full')/((np.std(conctime0[:,yout,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
                    conccorr[idxc,6] = conccorrarray[np.shape(Headinlettime0)[0]-1+np.argmax(np.abs(conccorrarray[np.shape(Headinlettime0)[0]-1:]))]
                    conccorr[idxc,7] = (np.argmax(np.abs(conccorrarray[np.shape(Headinlettime0)[0]-1:])))*5
                    amplitudec[idxc,:6] = conccorr[idxc, :6]
                    amplitudec[idxc,6] = normavgconcout[np.argmin(normavgconcout[1:,k]),k]/np.mean(conctime0[:,yout,k])
                    amplitudec[idxc,7] = normavgconcout[np.argmax(normavgconcout[1:,k]),k]/np.mean(conctime0[:,yout,k])
                    if gvarnames[k] in plotg:
                        axes.flat[plotg.index(gvarnames[k])].stem(conccorrarray[np.shape(Headinlettime0)[0]-1:])
                for k in range(len(AFbiomassgvarnames)):
                    normsum[:,k] = (sumall0[:,k] - np.mean(sumall0[:,k]))
                    spratio[:,k] = sumall0[:,k]/np.sum(sumall0[:,:], axis = -1)
                for k in range(len(AFbiomassgvarnames)):
                    normspratio[:,k] = spratio[:,k] - np.mean(spratio[:,k])
                for k in range(len(AFbiomassgvarnames)):
                    idxb = Regimes.index(Reg)*(len(Tforfpre)-1)*len(Trial)*(len(AFbiomassgvarnames)) + (Tforfpre.index(p)-1)*len(Trial)*(len(AFbiomassgvarnames)) + Trial.index(t)*(len(AFbiomassgvarnames))+k
                    if (t=='H'):
                        conccorrb[idxb,0] = Trial.index(t)
                    else:
                        conccorrb[idxb,0] = t
                    conccorrb[idxb,1] = Het[Trial.index(t)]
                    conccorrb[idxb,2] = Anis[Trial.index(t)]
                    conccorrb[idxb,3] = Tforfpre.index(p)
                    conccorrb[idxb,4] = Regimes.index(Reg)
                    conccorrb[idxb,5] = k
                    conccorrarrayb = np.correlate(normspratio[:,k],normheadin,'full')/((np.std(spratio[:,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
                    conccorrb[idxb,6] = conccorrarrayb[np.shape(Headinlettime0)[0]-1+np.argmax(np.abs(conccorrarrayb[np.shape(Headinlettime0)[0]-1:]))]
                    conccorrb[idxb,7] = (np.argmax(np.abs(conccorrarrayb[np.shape(Headinlettime0)[0]-1:])))*5
                    amplitudeb[idxb,:6] = conccorrb[idxb, :6]
                    amplitudeb[idxb,6] = normsum[np.argmin(normsum[1:,k]),k]/np.mean(sumall0[:,k])
                    amplitudeb[idxb,7] = normsum[np.argmax(normsum[1:,k]),k]/np.mean(sumall0[:,k])
    
    return amplitudec, conccorr, amplitudeb, conccorrb

def norm_amplitude(directory, numberofTempscenarios, Trial, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, AFbiomassvars, AFbiomassgvarnames):
    Regimes = ["Slow","Equal", "Fast"]
    Ampchem = np.zeros([len(Regimes)*len(Trial)*numberofTempscenarios*len(gvarnames),8])
    Ampbio = np.zeros([len(Regimes)*len(Trial)*numberofTempscenarios*len(AFbiomassgvarnames),8])
    for Reg in Regimes:
        print (Reg)
        Tforfpre = [Reg+'AR_0', Reg+'AR_1', Reg+'AR_2',Reg+'AR_5']
        for p in Tforfpre[1:]:
            print (p)
            for t in Trial:
                print (t)
                path = directory+p+fpre+str(t)+fsuf+fpre+str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
                sumall0 = calcsum_temp(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
                normsum = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                spratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                normspratio = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
                normavgconcout = np.zeros([np.shape(df0)[1],len(gvarnames)])
                normheadin = Headinlettime0 - np.mean(Headinlettime0)        
                for k in range(len(gvarnames)):
                    normavgconcout[:,k] = (conctime0[:,yout,k]/np.mean(conctime0[:,yout,k]))
                for k in range(len(gvarnames)):   
                    f, Pxx_spec = signal.periodogram(normavgconcout[:,k], scaling='spectrum')
                    idxc = Regimes.index(Reg)*(len(Tforfpre)-1)*len(Trial)*(len(gvarnames)) + (Tforfpre.index(p)-1)*len(Trial)*(len(gvarnames)) + Trial.index(t)*(len(gvarnames))+k
                    if (t=='H'):
                        Ampchem[idxc,0] = Trial.index(t)
                    else:
                        Ampchem[idxc,0] = t
                    Ampchem[idxc,1] = Het[Trial.index(t)]
                    Ampchem[idxc,2] = Anis[Trial.index(t)]
                    Ampchem[idxc,3] = Tforfpre.index(p)
                    Ampchem[idxc,4] = Regimes.index(Reg)
                    Ampchem[idxc,5] = k
                    Ampchem[idxc,6] = np.sqrt(Pxx_spec.max())
                    Ampchem[idxc,7] = np.sqrt(Pxx_spec.mean())
                for k in range(len(AFbiomassgvarnames)):
                    normsum[:,k] = (sumall0[:,k]/np.mean(sumall0[:,k]))
                    spratio[:,k] = sumall0[:,k]/np.sum(sumall0[:,:], axis = -1)
                for k in range(len(AFbiomassgvarnames)):
                    normspratio[:,k] = spratio[:,k]/np.mean(spratio[:,k])
                for k in range(len(AFbiomassgvarnames)):
                    f, Pxx_spec = signal.periodogram(normspratio[:,k], scaling='spectrum')
                    idxb = Regimes.index(Reg)*(len(Tforfpre)-1)*len(Trial)*(len(AFbiomassgvarnames)) + (Tforfpre.index(p)-1)*len(Trial)*(len(AFbiomassgvarnames)) + Trial.index(t)*(len(AFbiomassgvarnames))+k
                    if (t=='H'):
                        Ampbio[idxb,0] = Trial.index(t)
                    else:
                        Ampbio[idxb,0] = t
                    Ampbio[idxb,1] = Het[Trial.index(t)]
                    Ampbio[idxb,2] = Anis[Trial.index(t)]
                    Ampbio[idxb,3] = Tforfpre.index(p)
                    Ampbio[idxb,4] = Regimes.index(Reg)
                    Ampbio[idxb,5] = k
                    Ampbio[idxb,6] = np.sqrt(Pxx_spec.max())
                    Ampbio[idxb,7] = np.sqrt(Pxx_spec.mean())
                    
    return Ampchem, Ampbio
                    
                    
def head_amplitude (directory, numberofTempscenarios, Trial, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames):
    Regimes = ["Slow","Equal", "Fast"]
    amplitude = np.zeros([len(Trial)*(numberofTempscenarios)*len(Regimes),8])
    idxc = 0
    for Reg in Regimes:
        print (Reg)
        Tforfpre = [Reg+'AR_0', Reg+'AR_1', Reg+'AR_2',Reg+'AR_5']
        for p in Tforfpre[1:]:
            print (p)
            for t in Trial:
                print (t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                if (t=='H'):
                    amplitude[idxc,0] = Trial.index(t)
                else:
                    amplitude[idxc,0] = t
                amplitude[idxc,1] = Het[Trial.index(t)]
                amplitude[idxc,2] = Anis[Trial.index(t)]
                amplitude[idxc,3] = Tforfpre.index(p)
                amplitude[idxc,4] = Regimes.index(Reg)
                amplitude[idxc,5] = 0
                amplitude[idxc,6] = normheadin[np.argmin(normheadin[1:])]/np.mean(Headinlettime0)
                amplitude[idxc,7] = normheadin[np.argmax(normheadin[1:])]/np.mean(Headinlettime0)
                idxc = idxc + 1
   
    return amplitude