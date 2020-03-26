# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""
import numpy  as np
import csv
import Pythonfunctions_SK as sk
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt

#Saturated flow regime
Reg = "Slow"
d = r"Z:/Saturated_flow/diffusion_transient/"
fpre = '/NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]
Tforfpre = [Reg+'AR_0', Reg+'AR_1', Reg+'AR_2',Reg+'AR_5']
fsuf  = r"/"
gw = 1
filename = 'model_domain_quad.tec'

#Scan for files
for j in Tforfpre:
    print (j)
    for k in Trial:
            file_name = fpre[1:]+str(k)+"_df.npy"
            #np.load(d+j+fpre+str(k)+fsuf+fpre+str(k)+"_df.npy")
            cur_dir = d+Reg+"AR_"+j+fpre+str(k)+fsuf # Dir from where search starts can be replaced with any path
            file_list = os.listdir(cur_dir)
            if file_name in file_list:
#                print (file_name, " File Exists in: ", cur_dir)
                continue
            else:
                print (file_name, "File not found: ", cur_dir)
                
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

#Variations:
Trial = masterTrial[:masterTrial.index(60)]
Trial = []
notlist = [57, 75]
for i in masterTrial:
    if i not in notlist:
        Trial.append(i)
indices = list(masterTrial.index(i) for i in Trial)
Het = list(masterHet[i] for i in indices)
Anis = list(masterAnis[i] for i in indices)
#Constants
yout = 50
yin = 0
xleft = 0
xright = 30
#Assign index to Variable
doc1=10-gw
dox1=11-gw
Amm1=12-gw
nitra1=17-gw
sulpha1=22-gw
tr1 = 29-gw
Bfo1 = 8 - gw
Bfn1 = 15 - gw
Bfs1 = 20 - gw
Bfa1 = 25 - gw
Bmo1 = 9 - gw
Bmn1 = 16 - gw
Bms1 = 21 - gw
Bma1 = 26 - gw
Bifo1 = 13 - gw
Bifn1 = 18 - gw
Bifs1 = 23 - gw
Bifa1 = 27 - gw
Bimo1 = 14 - gw
Bimn1 = 19 - gw
Bims1 = 24 - gw
Bima1 = 28 - gw
vely=5
velx=4
vars = [doc1, dox1, Amm1, nitra1, sulpha1, Bmo1, Bma1, Bmn1, Bimo1, Bima1, Bimn1]
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate", "Sulphate", "Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers",
                      "Inactive mobile Aerobes", "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers","Nitrogen","TOC"]
AFbiomassvars = [Bfo1, Bfa1, Bfn1, Bmo1, Bma1, Bmn1, Bifo1, Bifa1, Bifn1, Bimo1, Bima1, Bimn1]
AFbiomassgvarnames = ["Active fixed Aerobes", "Active fixed Ammonia oxidizers", "Active fixed Nitrate reducers",
                      "Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers",
                      "Inactive fixed Aerobes", "Inactive fixed Ammonia oxidizers", "Inactive fixed Nitrate reducers",
                      "Inactive mobile Aerobes", "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers"]

for Reg in ["Slow","Equal"]:
    Tforfpre = [Reg+'AR_0',Reg+'AR_1', Reg+'AR_2', Reg+'AR_5']
    mft, calcsum = sk.calcmft_temp(Tforfpre, Trial, vars, gvarnames, AFbiomassvars, AFbiomassgvarnames, d, fpre, fsuf, Het, Anis, gw)
#    mftbuffer, calcsummob = sk.calcmft_temp(Tforfpre, Trial, vars, gvarnames, AMbiomassvars, AMbiomassgvarnames, d, fpre, fsuf, Het, Anis, gw)
    #Writing the results in a csv file     
    fnamemf = "X:/diffusion_"+Reg+"_Transient_"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"MF_"+str(gw) + ".csv"
    csvfile= open(fnamemf, "w")
    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Inlet_total_mass_flux", "Outlet_mass_flux", "Removal", "Removal_bc_Ratio", "RemovalRatio_homo_Time","RemovalRatio_hetero_Time","Time", "Regime"])
    idx = 0
    for l in range(len(Tforfpre)):
        for j in range(len(Trial)):
            for i in range(len(gvarnames)):
                idx = (l*len(Trial)+j)*len(gvarnames)+i
                if (Reg == "Equal"):
                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,10],mft[idx,8],mft[idx,9],"Medium"])
                else:
                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,10],mft[idx,8],mft[idx,9],Reg])
    csvfile.close()
        
    fnamebiosum = "X:/diffusion_"+Reg+"_Transient"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"biomass_"+str(gw) +".csv"
    csvfile= open(fnamebiosum, "w")
    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_biomass_umoles", "TotalRatio_bc", "TotalRatio_hetero_Time",  "TotalRatio_homo_Time", "Species_fraction","SpeciesRatio_bc", "SpeciesRatio_hetero_Time", "SpeciesRatio_homo_Time","Time", "Regime"])
    idxb = 0
    for l in range(len(Tforfpre)):
        for j in range(len(Trial)):
            for i in AFbiomassgvarnames:
                idx = (l*len(Trial) + j)*len(AFbiomassgvarnames)+AFbiomassgvarnames.index(i)
                if (Reg == "Equal"):
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7],calcsum[idx,9],calcsum[idx,10],calcsum[idx,11],calcsum[idx,12],calcsum[idx,8], "Medium"])
                else:
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7],calcsum[idx,9],calcsum[idx,10],calcsum[idx,11],calcsum[idx,12],calcsum[idx,8], Reg])
    csvfile.close()
    
    print ("Files written, processing as dataframes ...")
    if Reg == "Fast":
#        Reg = "Fast"
#        fnamemf = "Fast_Transient_H_59_MF_1_time.csv"
#        fnamebiosum = "Fast_Transient_AnaerobicH_59_biomass_1_time.csv"
        fastbc = sk.processchembiomassfiles(fnamebiosum, Reg)
#        fastmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        fastc = sk.processchemfiles(fnamemf, Reg)
    elif Reg == "Slow":
#        Reg = "Slow"
#        fnamemf = "X:/Saturated_flow/Anaerobic_Slow_Transient_H_59_MF_1.csv"
#        fnamebiosum = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59_biomass_1.csv"
#        fnamebiosummob = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59__mobile_biomass_1.csv"
        slowbc = sk.processchembiomassfiles(fnamebiosum, Reg)
#        slowmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        slowc = sk.processchemfiles(fnamemf, Reg)
    elif Reg == "Equal":
#        Reg = "Equal"
#        fnamemf = "Equal_Transient_H_59_MF_1_time.csv"
#        fnamebiosum = "Equal_Transient_AnaerobicH_59_biomass_1_time.csv"
        equalbc = sk.processchembiomassfiles(fnamebiosum, "Medium")
#        equalmbc = sk.processchembiomassfiles(fnamebiosummob, "Medium")
        equalc = sk.processchemfiles(fnamemf, "Medium")
    print ("Dataframes processed, move to next series")


#Aggregated results
#Importing data

#Mass flux and biomass with average residence time
dummy = sk.scatterrestime_flux_temp(0, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "Removal_bc_Ratio", "Ratio of removal with base case in uniform flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_uf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_uf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_flux_temp(0, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_homo_Time", "Ratio of removal with base case in varying flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_flux_temp(0, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_hetero_Time", "Ratio of removal with heterogeneous case in varying flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_hc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_hc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_flux_temp_singleaxis(0, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_homo_Time", "Ratio of removal with base case in varying flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bv_vf_1col.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf_1col.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_biomass_temp(0, slowbc, equalbc, AFbiomassgvarnames, "SpeciesRatio_hetero_Time", "Ratio of species fraction with heterogeneous case in varying flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_hc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_hc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_biomass_temp(0, slowbc, equalbc, AFbiomassgvarnames, "SpeciesRatio_homo_Time", "Ratio of species fraction with homogeneous case in varying flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy = sk.scatterrestime_biomass_temp(0, slowbc, equalbc, AFbiomassgvarnames, "SpeciesRatio_bc", "Ratio of species fraction with homogeneous case in uniform flow rate")
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_uf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_uf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)

#Normalizing time series for correlation
from scipy.signal import correlate
Regimes = ["Slow","Equal"]
corr = np.zeros([len(Trial)*(len(Tforfpre)-1)*len(Regimes)*(len(gvarnames)+1),8])
conccorr = np.zeros([len(Trial)*(len(Tforfpre)-1)*len(Regimes)*(len(gvarnames)),8])
conccorrb = np.zeros([len(Trial)*(len(Tforfpre)-1)*len(Regimes)*(len(AFbiomassgvarnames)),8])
for Reg in Regimes:
    print (Reg)
    Tforfpre = [Reg+'AR_0', Reg+'AR_1', Reg+'AR_2',Reg+'AR_5']
    for p in Tforfpre[1:]:
        print (p)
        for t in Trial:
#            print (t)
#            numpyfile = d+p+fpre+str(t)+fpre+str(t)+"_df.npy"
#            df = np.load(numpyfile)
##            massfluxout = np.zeros([len(vars)+1,np.shape(df)[1]])
#            normmassfluxout = np.zeros([len(vars)+1,np.shape(df)[1]])
#            massfluxout[-1,:] = np.sum(df[vely-3,:,yin,xleft+1:xright], axis = -1)*0.005*0.01/por + (df[vely-3,:,yin,xleft]+df[vely-3,:,yin,xright])*0.005*0.005/por
#            avgconcout = np.zeros([len(vars),np.shape(df)[1]-1])
            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sk.calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, d+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
            sumall0 = sk.calcsum_temp(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, d+p, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
            normavgconcout = np.zeros([np.shape(df0)[1],len(gvarnames)])
            normsum = np.zeros([np.shape(df0)[1]-1,len(AFbiomassgvarnames)])
            for k in range(len(gvarnames)):
               normavgconcout[:,k] = (conctime0[:,yout,k] - np.mean(conctime0[:,yout,k]))
            normheadin = Headinlettime0 - np.mean(Headinlettime0)
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
                conccorrarray = correlate(normheadin,normavgconcout[:,k],'full')/((np.std(conctime0[:,yout,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
#                print(conccorrarray)
#                concabsmaxcorr = np.max(np.abs((conccorrarray)))
#                print(concabsmaxcorr)
#                concdelay = np.where(np.abs(conccorrarray) == concabsmaxcorr)[0]
#                print (concdelay)
#                concmaxcorr = conccorrarray[concdelay]
                conccorr[idxc,6] = np.max(conccorrarray)
                conccorr[idxc,7] = (np.argmax(conccorrarray) - (np.shape(df0)[1]-1))*5
            for k in range(len(AFbiomassgvarnames)):
                normsum[:,k] = (sumall0[:,k] - np.mean(sumall0[:,k]))
#            normheadin = Headinlettime0 - np.mean(Headinlettime0)
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
                conccorrarrayb = correlate(normheadin,normsum[:,k],'full')/((np.std(sumall0[:,k]))*(np.std(Headinlettime0)*np.shape(Headinlettime0)[0]))
#                print(conccorrarray)
#                concabsmaxcorrb = np.max(np.abs((conccorrarrayb)))
#                print(concabsmaxcorr)
#                concdelayb = np.where(np.abs(conccorrarrayb) == concabsmaxcorrb)[0]
#                print (concdelay)
 #               concmaxcorrb = conccorrarrayb[concdelayb]
                conccorrb[idxb,6] = np.max(conccorrarrayb)
                conccorrb[idxb,7] = (np.argmax(conccorrarray) - (np.shape(df0)[1]-1))*5
                

corrdf = pd.DataFrame(conccorr, columns = ["Trial", "Variance", "Anisotropy", "Time", "Regime", "Chem","Correlation","Delay"])
                                       
#                                       "DOC", "DO", "Ammonium", 
#                                       "Nitrate", "Sulphate", "Active mobile Aerobes", "Active mobile Ammonia oxidizers",
#                                       "Active mobile Nitrate reducers","Inactive mobile Aerobes",
#                                       "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers", "Velocity"])
bth1, bth = sk.tracerstudies()
corrdf["Regime"] = corrdf["Regime"].astype(int)
corrdf["Regime"] = corrdf["Regime"].replace([0],"Slow")
corrdf["Regime"] = corrdf["Regime"].replace([1],"Medium")
corrdf["Regime"] = corrdf["Regime"].astype(str)
corrdf["Trial"] = corrdf["Trial"].astype(int)
corrdf["Trial"] = corrdf["Trial"].replace([0], "H")
corrdf["Trial"] = corrdf["Trial"].astype(str)
corrdf["Chem"]= corrdf["Chem"].astype(int)
for k in range(len(AFbiomassgvarnames)):
    corrdf["Chem"] = corrdf["Chem"].replace([k], AFbiomassgvarnames[k])
corrdf["Chem"] = corrdf["Chem"].replace([11], "Velocity")
corrdf["Chem"]= corrdf["Chem"].astype(str)
    
dfall2 = pd.merge(corrdf, bth[['Trial', 'Regime', 'fraction']], on = ['Trial', 'Regime'])
colseries = ["Reds", "Greens", "Blues"]
corrdffilter=dfall2[dfall2["Trial"]!='75']
corrdffilter2=corrdffilter[corrdffilter["Time"]>0]
corrdf2=corrdffilter2[corrdffilter2["Trial"]!='57']
ncols = 1
Regimes = ["Slow", "Medium"]

c = ["Correlation", "Delay"]

for ctype in c:
    Chems = ["DOC", "DO", "Nitrogen", "TOC"]
    Chems = ["Active fixed Aerobes", "Active fixed Ammonia oxidizers", "Active fixed Nitrate reducers"]
    nrows = len(Chems)
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = [4,10])
    for i in Regimes:
        dfr = corrdf2[corrdf2["Regime"]==i]
        for c in Chems:
            dfrc = dfr[corrdf2["Chem"]==c]
            ax.flat[Chems.index(c)].scatter(x = "fraction", y = ctype, c = "Time", cmap = colseries[Regimes.index(i)], data =dfrc)
    ax.flat[Chems.index(c)].set_xlabel("Ratio of breakthrough time with the base case", fontsize = 15)
    for a,typsp in zip(ax, Chems):
        a.set_ylabel(typsp, fontsize = 15)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_AIbiomass_head_normconc.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_AIbiomass_head_normconc.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

    Chems = ["Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers"]
    Chems = ["Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers"]
    nrows = len(Chems)
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = [4,15])
    for i in Regimes:
        dfr = corrdf2[corrdf2["Regime"]==i]
        for c in Chems:
            dfrc = dfr[corrdf2["Chem"]==c]
            ax.flat[Chems.index(c)].scatter(x = "fraction", y = ctype, c = "Time", cmap = colseries[Regimes.index(i)], data = dfrc)
    ax.flat[Chems.index(c)].set_xlabel("Ratio of breakthrough time with the base case", fontsize = 15)
    for a,typsp in zip(ax, Chems):
        a.set_ylabel(typsp, fontsize = 15)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_AMbiomass_head_normconc.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_AMbiomass_head_normconc.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

    Chems = ["Inactive mobile Aerobes","Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers"]
    nrows = len(Chems)
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = [4,15])
    for i in Regimes:
        dfr = corrdf2[corrdf2["Regime"]==i]
        for c in Chems:
            dfrc = dfr[corrdf2["Chem"]==c]
            ax.flat[Chems.index(c)].scatter(x = "fraction", y = ctype, c = "Time", cmap = colseries[Regimes.index(i)], data = dfrc)
    ax.flat[Chems.index(c)].set_xlabel("Ratio of breakthrough time with the base case", fontsize = 15)
    for a,typsp in zip(ax, Chems):
        a.set_ylabel(typsp, fontsize = 15)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_inactivemobilebiomass_head_normconc.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
    plt.savefig("Z:/Saturated_flow/diffusion_transient/"+ctype+"_inactivemobilebiomass_head_normconc.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

fnamemf = "X:/Saturated_flow/time_v2_sumalltime_correlation"+Reg+"_Transient_"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"MF_"+str(gw) + ".csv"
csvfile= open(fnamemf, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Time", "Regime", "Chem", "Max_correlation", "Delay", "Breakthrough_time_fraction"])
for j in range(len(corrdf2)):
    writer.writerow([j, dfall2["Trial"][j], dfall2["Variance"][j], dfall2["Anisotropy"][j],
                                dfall2["Time"][j], dfall2["Regime"][j], dfall2["Chem"][j],
                                dfall2["Correlation"][j], dfall2["Delay"][j],dfall2["fraction"][j]])
csvfile.close()

Regimes = ["Slow", "Equal"]
intsce = ['H', 44, 76, 73, 80, 84, 63]
for i in intsce:
    initseries = [500, 430, 600] 
    lastseries = [700, 630, 800]
    sk.generate_timeseries(Regimes, initseries, lastseries, Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, vars, gvarnames, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)

fig,ax = plt.subplots(nrows = len(Regimes), ncols = 1, figsize = [8,10])
for a, Reg in zip(ax, Regimes):
    f = "chemsamndbiomasswithvel_temp_3_"+Reg+"_H_ZOOMED.png"
    im = plt.imread(d+f)
    a.imshow(im)
    if Reg=="Equal":
        a.annotate("Medium flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
    else:
        a.annotate(Reg+" flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
#    a.set_title (Reg+" flow", fontsize = 15)
    a.axis('off')
plt.savefig("//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_chem_time.pdf", dpi = 300, bbox_inches = 'tight')
plt.savefig("//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_chem_time.png", dpi = 300, bbox_inches = 'tight')

intsce = [80, 84, 73, 63]
fig,ax = plt.subplots(nrows = 2, ncols = 2, figsize = [15,10])
axlinear = ax.flatten()
for a, i in zip(axlinear,intsce):
    f = "chemsamndbiomasswithvel_temp_3_Equal_"+str(i)+"_ZOOMED.png"
    im = plt.imread(d+f)
    a.imshow(im)
    a.axis('off')
plt.savefig("//tsclient/D/Saturated_flow/diffusion/Medium_heterogeneous_profiles_chem_time.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig("//tsclient/D/Saturated_flow/diffusion/Medium_Lheterogeneous_profiles_chem_time.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)

for Reg in Regimes:
    newd = d + Reg + "AR_"
    for j in Tforfpre[1:]:
        #Concentration profiles
        intsce = ["H", 76, 73, 80, 84, 63]
        for chem in ["DO", "DOC"]:
            dum = sk.shiftoxic_anoxic_temp (chem, Trial, intsce, newd, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars)
            picname = "Z:/Saturated_flow/diffusion_transient/chemprofile_temp_"+chem+str(Tforfpre.index(j))+"_"+Reg+"_.png"
            dum.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0)
                
        dummy = sk.concprofilewithtime_temp (Trial, intsce, newd, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars)
        picname = "Z:/Saturated_flow/diffusion_transient/chemswithvel_temp_"+str(Tforfpre.index(j))+"_"+Reg+"_.png"
        dummy.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0)

fig,ax = plt.subplots(ncols = len(Regimes), nrows = 1, figsize = [14,8])
for a, Reg in zip(ax, Regimes):
    f = "chemprofile_temp_DO3_"+Reg+"_.png"
    im = plt.imread(d+f)
    a.imshow(im)
    if Reg=="Equal":
#        a.annotate("Medium flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
        a.set_title ("Medium flow", fontsize = 10)
    else:
#        a.annotate(Reg+" flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
        a.set_title (Reg+" flow", fontsize = 10)
#    a.set_title (Reg+" flow", fontsize = 15)
    a.axis('off')
plt.savefig("//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_DO_time.pdf", dpi = 300, bbox_inches = 'tight')
plt.savefig("//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_DO_time.png", dpi = 300, bbox_inches = 'tight')
    
#calculating average velocities
for Reg in ["Slow","Equal"]:
    print (Reg)
    Tforfpre = [Reg+'AR_0',Reg+'AR_1', Reg+'AR_2', Reg+'AR_5']
    for j in Tforfpre:
        print(j)
        if (Tforfpre.index(j)==0):
            pass
        else:
            intsce = [43,44,45,52,53,54,61,62,63,73,74]
            for i in intsce:
                newd = d + j
                df,massendtime, masstime, conctime, Velocity, head = sk.calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
                Velocity = df[2,1:,:,:]  
                print(i, ": ",np.mean(Velocity))