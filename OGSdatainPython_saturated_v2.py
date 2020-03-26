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

#Saturated flow regime
Reg = "Fast"
d = r"X:/Saturated_flow/diffusion/"+ Reg+ "AR_changedkindox/"
fpre = 'NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]

fsuf  = r"/"
gw = 1

filename = 'model_domain_quad.tec'

#setup what we really want to investigate
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

#Variations:
#Trial = [45]
#indices = list(masterTrial.index(i) for i in Trial)
#Het = list(masterHet[i] for i in indices)
#Anis = list(masterAnis[i] for i in indices)
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
vars = [doc1, dox1, Amm1, nitra1, sulpha1]
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate", "Sulphate", "Nitrogen"]
AFbiomassvars = [Bfo1, Bfa1, Bfn1]
AFbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AMbiomassvars = [Bmo1, Bma1, Bmn1]
AMbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]

#Reading and storing in numpy array
for j in range(len(Trial)):
    print (str(Trial[j]))
    di =  d+fpre+str(Trial[j])+fsuf
    fwithd = di+filename
    print("Reading tech file....")
    size, steps, Headers, D = sk.readTecfile(fwithd)
    print("Converting to array....")
    df = sk.Converttomarr(D)
    print("Saving numpy array....")
#    np.save(di+fpre+str(Trial[j])+'_D', D)
    np.save(di+fpre+str(Trial[j])+'_df', df)
    #Test for correct orientation of the data
    for i in range(np.shape(df)[0]):
        print(Headers[i+3], np.mean(df[i,steps-1,0,:]))

mastermf = np.zeros([1,8])
masterbiomasssum = np.zeros([1,6])
masterbiomasssummob = np.zeros([1,6])
Regimes = ["Slow", "Equal", "Fast"]
for Reg in Regimes:
    d = r"X:/Saturated_flow/Anaerobic/"+Reg+ "AR_changedkindox/"
    fpre = "NS-A"
    mf = sk.calcmassflux (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
    biomasssum = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
    biomasssummob = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AMbiomassvars, AMbiomassgvarnames)
    mastermf = np.append(mastermf, mf, axis = 0)
    masterbiomasssum = np.append(masterbiomasssum, biomasssum, axis =0)
    masterbiomasssummob = np.append(masterbiomasssummob, biomasssummob, axis =0)
mastermf = np.delete(mastermf, 0, 0)
masterbiomasssum = np.delete(masterbiomasssum, 0, 0)
masterbiomasssummob = np.delete(masterbiomasssummob, 0, 0)

del4massflux = np.zeros([882,1])
del2biomass = np.zeros([441,1])
del2biomassmob = np.zeros([441,1])
regindex = np.zeros([882,1])
regindexb = np.zeros([441,1])
for Reg in Regimes:
    for k in Trial:
        for i in gvarnames:
            idx = Regimes.index(Reg)*len(Trial)*len(gvarnames) + Trial.index(k)*len(gvarnames) + gvarnames.index(i)
            del4massflux[idx,0] = (mastermf[idx,6]-mastermf[gvarnames.index(i),6])/mastermf[gvarnames.index(i),6]
            mastermf[idx, 0] = int(Trial.index(k)+36)
        for b in AFbiomassgvarnames:
            idxb = Regimes.index(Reg)*len(Trial)*len(AFbiomassgvarnames) + Trial.index(k)*len(AFbiomassgvarnames) + AFbiomassgvarnames.index(b)
            del2biomass[idxb,0] = (masterbiomasssum[idxb,4]-masterbiomasssum[AFbiomassgvarnames.index(b),4])/masterbiomasssum[AFbiomassgvarnames.index(b),4]
            del2biomassmob[idxb,0] = (masterbiomasssummob[idxb,4]-masterbiomasssummob[AFbiomassgvarnames.index(b),4])/masterbiomasssummob[AFbiomassgvarnames.index(b),4]
            masterbiomasssummob[idxb, 0] = int(Trial.index(k)+36)
            masterbiomasssum[idxb, 0] = int(Trial.index(k)+36)
print (idx, idxb)

for Reg in Regimes:
    for k in Trial:
        for i in gvarnames:
            idx = Regimes.index(Reg)*len(Trial)*len(gvarnames) + Trial.index(k)*len(gvarnames) + gvarnames.index(i)
            regindex[idx,0] = Regimes.index(Reg)
        for b in AFbiomassgvarnames:
            idxb = Regimes.index(Reg)*len(Trial)*len(AFbiomassgvarnames) + Trial.index(k)*len(AFbiomassgvarnames) + AFbiomassgvarnames.index(b)
            regindexb[idxb,0] = Regimes.index(Reg)


mastermf = np.append(mastermf, del4massflux, 1)
masterbiomasssum = np.append(masterbiomasssum, del2biomass, 1)
masterbiomasssummob = np.append(masterbiomasssummob, del2biomassmob, 1)
mastermf = np.append(mastermf, regindex, 1)
masterbiomasssum = np.append(masterbiomasssum, regindexb, 1)
masterbiomasssummob = np.append(masterbiomasssummob, regindexb, 1)

dfallmf = pd.DataFrame(data = mastermf, columns = ["Trial", "Variance", "Anisotropy", "Chem", "Inlet_total_mass_flux", "Outlet_mass_flux", "delmassflux", "del2massflux", "del4massflux", "Regime"])
dfallbiomass = pd.DataFrame(data = masterbiomasssum, columns = ["Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles", "del2biomass", "Regime"])
dfallbiomassmob = pd.DataFrame(data = masterbiomasssummob, columns = ["Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles", "del2biomass", "Regime"])

dfallmf["Regime"] = dfallmf["Regime"].replace([0, 1, 2], ["Slow", "Medium", "Fast"])
dfallmf["Regime"] = dfallmf["Regime"].astype(str)
dfallmf["Trial"] = dfallmf["Trial"].astype(int)
dfallmf["Trial"] = dfallmf["Trial"].replace([36], "H")
dfallmf["Trial"] = dfallmf["Trial"].astype(str)
dfallbiomass["Regime"] = dfallbiomass["Regime"].replace([0, 1, 2], ["Slow", "Medium", "Fast"])
dfallbiomass["Regime"] = dfallbiomass["Regime"].astype(str)
dfallbiomass["Trial"] = dfallbiomass["Trial"].astype(int)
dfallbiomass["Trial"] = dfallbiomass["Trial"].replace([36], "H")
dfallbiomass["Trial"] = dfallbiomass["Trial"].astype(str)
dfallbiomassmob["Regime"] = dfallbiomassmob["Regime"].replace([0, 1, 2], ["Slow", "Medium", "Fast"])
dfallbiomassmob["Regime"] = dfallbiomassmob["Regime"].astype(str)
dfallbiomassmob["Trial"] = dfallbiomassmob["Trial"].astype(int)
dfallbiomassmob["Trial"] = dfallbiomassmob["Trial"].replace([36], "H")
dfallbiomassmob["Trial"] = dfallbiomassmob["Trial"].astype(str)

filename = "X:/Saturated_flow/Steady_state/Tracer_studies/Combined_tracer.csv"
breakthrough = pd.read_csv(filename, delimiter = ",")
breakthrough ['Heterogeneity'] = breakthrough ['Variance'] 
breakthrough['VA']=breakthrough['Heterogeneity']*breakthrough['Anisotropy']
breakthrough['%ofhomogeneous']=breakthrough['del']*100
breakthrough['%sfraction']=breakthrough['fraction_withslow']*100
breakthrough['%efraction']=breakthrough['fraction_withequal']*100
breakthrough['%ffraction']=breakthrough['fraction_withfast']*100

dfall2 = pd.merge(dfallmf, breakthrough[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous', 'fraction_withslow', '%sfraction', '%efraction', '%ffraction']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
dfall2["%del4massflux"] = dfall2["del4massflux"]*100
Regimes = ["Slow", "Medium", "Fast"]
Chems = ["DOC", "DO", "Nitrogen"]
colseries = ["indianred", "g", "steelblue"]
ncols = len(Regimes)
nrows = len(Chems)

import matplotlib.pyplot as plt
fig, axes = plt.subplots(nrows = nrows, figsize = [10,10])
plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state", fontsize = 20)
for k in Chems:
    dfc = dfall2[dfall2['Chem'].astype(int)==Chems.index(k)]
    colidx1 = Chems.index(k)
    pred = dfc["%sfraction"].to_numpy()
    resp = dfc["%del4massflux"].to_numpy()
    line= np.polyfit(pred, resp,2)
    print (line)
    for i in Regimes:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp['Regime']==i]
        print (i)
        axes.flat[colidx1].scatter("%sfraction", "%del4massflux", color = colseries[Regimes.index(i)], data = dfcr)
#        axes.flat[colidx1].plot(pred, slope + intercept*pred, '--')
        axes.flat[colidx1].set_ylabel(k, fontsize = 15)
        axes.flat[colidx1].tick_params(axis = 'y', labelsize = 15)
        axes.flat[colidx1].set_xscale('log')
        if(Chems.index(k) != len(Chems)-1):
            axes.flat[colidx1].set_xlabel('')
            axes.flat[colidx1].set_xticklabels([])
        else:
            axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
            axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case in slow flow regime (%)', fontsize = 15)
fig.savefig("X:/Saturated_flow/Anaerobic/steadystate_mf_together_nitrogen.png")

dfall2 = pd.merge(dfallbiomass, breakthrough[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous', '%sfraction', '%efraction', '%ffraction']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
dfall2["%del2biomass"] = dfall2["del2biomass"]*100
Regimes = ["Slow", "Medium", "Fast"]
Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
colseries = ["indianred", "g", "steelblue"]
ncols = len(Regimes)
nrows = len(Chems)
fig, axes = plt.subplots(nrows = nrows, figsize = [10,10])
plt.suptitle("Change in total biomass with respect to homogeneous scenario at steady state", fontsize = 20)
for k in Chems:
    dfc = dfall2[dfall2['Chem'].astype(int)==Chems.index(k)]
    colidx1 = Chems.index(k)
    pred = dfc["%sfraction"].to_numpy()
    resp = dfc["%del2biomass"].to_numpy()
    line= np.polyfit(pred, resp,2)
    print (line)
    for i in Regimes:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp['Regime']==i]
        print (i)
        axes.flat[colidx1].scatter("%sfraction", "%del2biomass", color = colseries[Regimes.index(i)], data = dfcr)
#        axes.flat[colidx1].plot(pred, slope + intercept*pred, '--')
        axes.flat[colidx1].set_ylabel(k, fontsize = 15)
        axes.flat[colidx1].tick_params(axis = 'y', labelsize = 15)
        axes.flat[colidx1].set_xscale('log')
        if(Chems.index(k) != len(Chems)-1):
            axes.flat[colidx1].set_xlabel('')
            axes.flat[colidx1].set_xticklabels([])
        else:
            axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
            axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case in slow flow regime (%)', fontsize = 15)
fig.savefig("X:/Saturated_flow/Anaerobic/steadystate_biomass_together.png")


dfall2 = pd.merge(dfallbiomassmob, breakthrough[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous', '%sfraction', '%efraction', '%ffraction']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
dfall2["%del2biomass"] = dfall2["del2biomass"]*100
Regimes = ["Slow", "Medium", "Fast"]
Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
colseries = ["indianred", "g", "steelblue"]
ncols = len(Regimes)
nrows = len(Chems)
fig, axes = plt.subplots(nrows = nrows, figsize = [10,10])
plt.suptitle("Change in total biomass with respect to homogeneous scenario at steady state", fontsize = 20)
for k in Chems:
    dfc = dfall2[dfall2['Chem'].astype(int)==Chems.index(k)]
    colidx1 = Chems.index(k)
    pred = dfc["%sfraction"].to_numpy()
    resp = dfc["%del2biomass"].to_numpy()
    line= np.polyfit(pred, resp,2)
    print (line)
    for i in Regimes:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp['Regime']==i]
        print (i)
        axes.flat[colidx1].scatter("%sfraction", "%del2biomass", color = colseries[Regimes.index(i)], data = dfcr)
#        axes.flat[colidx1].plot(pred, slope + intercept*pred, '--')
        axes.flat[colidx1].set_ylabel(k, fontsize = 15)
        axes.flat[colidx1].tick_params(axis = 'y', labelsize = 15)
        axes.flat[colidx1].set_xscale('log')
        if(Chems.index(k) != len(Chems)-1):
            axes.flat[colidx1].set_xlabel('')
            axes.flat[colidx1].set_xticklabels([])
        else:
            axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
            axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case in slow flow regime (%)', fontsize = 15)
fig.savefig("X:/Saturated_flow/Anaerobic/steadystate_mobile_biomass_together.png")


mf = sk.calcmassflux (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
biomasssum = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
biomasssummob = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AMbiomassvars, AMbiomassgvarnames)
meanresultsarr, stdresultsarr, covresultsarr = sk.calcaggres (mf)

#Writing the results into csv files
fnamemf = d+fpre+"Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"MF_"+str(gw) + datetime.now().strftime("%d%m%Y%H%M")+".csv"
csvfile= open(fnamemf, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Inlet_total_mass_flux", "Outlet_mass_flux", "delmassflux", "del2massflux"])
for j in Trial:
    for i in gvarnames:
        idx = Trial.index(j)*len(gvarnames)+gvarnames.index(i)
        writer.writerow([idx+1, j, Het[Trial.index(j)], Anis [Trial.index(j)], i, mf[idx,4], mf[idx,5], mf[idx,6], mf[idx,7]])
csvfile.close()

fnamebiosum = d+fpre+"Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"biomass_"+str(gw) + datetime.now().strftime("%d%m%Y%H%M")+".csv"
csvfile= open(fnamebiosum, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles"])
for j in Trial:
    for i in AFbiomassgvarnames:
        idx = Trial.index(j)*len(AFbiomassgvarnames)+AFbiomassgvarnames.index(i)
        writer.writerow([idx+1, j, Het[Trial.index(j)], Anis [Trial.index(j)], i, biomasssum[idx,4], biomasssum[idx,5]])
csvfile.close()

fnamebiosummob = d+fpre+"Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"biomass_"+str(gw) + datetime.now().strftime("%d%m%Y%H%M")+".csv"
csvfile= open(fnamebiosummob, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles"])
for j in Trial:
    for i in AMbiomassgvarnames:
        idx = Trial.index(j)*len(AMbiomassgvarnames)+AMbiomassgvarnames.index(i)
        writer.writerow([idx+1, j, Het[Trial.index(j)], Anis [Trial.index(j)], i, biomasssum[idx,4], biomasssum[idx,5]])
csvfile.close()

fnamemf_summary = d+fpre+"Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"new_MF_"+str(gw) +datetime.now().strftime("%d%m%Y%H%M")+ "_summary.csv"
csvfile= open(fnamemf_summary, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno", "Variance", "Anisotropy", "Chem", "Removal", "Change_removalrate", "Std_var_RR"])
for j in range(len(meanresultsarr)):
        writer.writerow([j+1, meanresultsarr[j,2], meanresultsarr[j,3], meanresultsarr[j,4], meanresultsarr[j,0]*100, meanresultsarr[j,1]*100, stdresultsarr[j,0]])
csvfile.close()

#1D analysis:
#Importing data
if Reg == "Fast":
#    Reg = "Fast"
#    fnamemf = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_84_MF_1030120201028.csv"
#    fnamebiosum = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_84_biomass_1030120201028.csv"
    fastb = sk.processchembiomassfiles(fnamebiosum, Reg)
    fastc = sk.processchemfiles(fnamemf, Reg)
elif Reg == "Slow":
    Reg = "Slow"
#    fnamemf = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_83_MF_1070120201302.csv"
#    fnamebiosum = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_83_biomass_1070120201707.csv"
    slowb = sk.processchembiomassfiles(fnamebiosum, Reg)
    slowc = sk.processchemfiles(fnamemf, Reg)
elif Reg == "Equal":
    Reg = "Equal"
#    fnamemf = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_84_MF_1030120201027.csv"
#    fnamebiosum = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_84_biomass_1030120201027.csv"
    equalb = sk.processchembiomassfiles(fnamebiosum, "Medium")
    equalc = sk.processchemfiles(fnamemf, "Medium")
    
#Scatter plot Flux with residence time
sgvarnames = [elem for elem in gvarnames if elem != "Sulphate"] 
sngvarnames = [elem for elem in sgvarnames if elem != "Nitrogen"] 
comb,dummy = sk.scatterrestime_flux(fastc, slowc, equalc, sngvarnames)
picname = "X:/Saturated_flow/Anaerobic/flux_restime_SCATTER_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname, bbox_inches = 'tight')
#Scatter plot biomass wtih residence time
comb,dummy = sk.scatterrestime_biomass(fastb, slowb, equalb)
picname = "X:/Saturated_flow/Anaerobic/mobile_biomass_restime_SCATTER_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
    
#Concentration profiles
nrows = 8
ncols = 6
figsize = [21,28]
colors = ['black','red', 'blue', 'green', 'orange', 'grey']
#Plot all the scenarios excluding homogeneous
intsce = Trial
sk.plotconcallatss(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors, nrows, ncols, figsize)

#changing oxic-anoxic interface with heterogeneity
intsce = ["H", 50, 76, 73, 80, 84, 44, 63]
colorseries = ["Reds", "Greens", "Blues"]
chem = "DO"
Regimes = ["Slow", "Equal", "Fast"]
#def concprofileatss (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors):
dummy = sk.plotconcwithhet (Regimes, intsce, chem, colorseries, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
picname = "X:/Saturated_flow/diffusion/"+str(intsce[0])+ "_" + str(intsce[-1])+"_"+datetime.now().strftime("%d%m%Y%H%M")+"_"+chem+"with_het.png"
dummy.savefig(picname)
        
#Zoom into select heterogeneous scenarios
intsce = ['H']
#intsce = Trial
#intsce = ["H"]
sk.concprofileatss(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors)
sk.concprofileatss(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames, intsce, colors)

#Zoom into select heterogeneous scenarios
intsce = [50, 76, 73, 80, 84, 44, 63]
for i in intsce:
    chemspecies = sk.concprofilewithH(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors,i)
    picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_chem_concwithH.png"
    chemspecies.savefig(picname, bbox_inches = 'tight')
    biodist = sk.concprofilewithH(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames, intsce, colors, i)
    picname = d+fpre+str(Trial[Trial.index(i)])+ "_" + datetime.now().strftime("%d%m%Y%H%M")+"_biomass_conc_withH.png"
    biodist.savefig(picname, bbox_inches = 'tight')
        
intsce = ['H', 63, 84, 76]

sk.concprofileatssX(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, colors)
sk.concprofileatssX(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, biomassvars, biomassgvarnames, intsce, colors)

#Plot tracer with time
intsce = [37,38,39,46,47,48,55,56,57,64,65,66,40,41,42,49,50,51,58,59,60,67,68,69,70,71,72,76,77,78,79,80,81,82,83,84,43,44,45,52,53,54,61,62,63,73,74,75]
inttime = [2,5,8,15,18,20]
chem = "DO"
lylim = 0
uylim = 5
#colors = "Red"
nrows = 2
ncols = 2
sk.plottimeseries(Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames, intsce, inttime, colors, nrows, ncols, figsize, chem, lylim, uylim)

#Classify oxic cells and then plot along Y axis
fig = sk.plotoxiccellssdo(20, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
picname = d+"DOandoxiccells_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
fig.savefig(picname)

for k in intsce:
    df = np.load(d+fpre+str(k)+fsuf+fpre+str(k)+'_df.npy')
    title = "Variance "+str(Het[Trial.index(k)])+" : Anisotropy "+str(Anis[Trial.index(k)])
    sk.heatmapconcdist (df, vars, k, gvarnames, d, fpre, title)
    sk.heatmapconcdist (df, biomassvars, k, biomassgvarnames, d, fpre, title)

start = 0
end = 15
chem = "DOC"
#Calculating average velocities
intsce = ["H", 37,38,39,46,47,48,55,56,57,64,65,66,40,41,42,49,50,51,58,59,60,67,68,69,70,71,72,76,77,78,79,80,81,82,83,84,43,44,45,52,53,54,61,62,63,73,74,75]
for j in range(len(Trial)):
    df,massendtime, masstime, conctime, Velocity = sk.calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    print ("Scenario:", Trial[j], ", Velocity: ", Velocity)

 
tr1 = 8-gw
vely=5
velx=4
vars = [tr1]

for j in range(len(Trial)):
    df,massendtime, masstime, conctime, Velocity = sk.calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    Time = np.where(np.round(conctime[:,yout, 1],3)>10)
    Time2 = np.where(np.round(df[-1, :, 50,:],3)>10)
    print (0.01*(1+Time[0][0]), 0.01*(1+Time2[0][0]))


#Boxplots:
#Mass flux with variance
dummy = sk.boxV_Aflux (fastc, slowc, equalc, gvarnames, imgsize = [24,18])
picname = "X:/Saturated_flow/Anaerobic/flux_onlyvariance_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Biomass with variance
dummy = sk.boxV_Abio (fastbc, slowbc, equalbc, imgsize = [22,15])
picname = "X:/Saturated_flow/Anaerobic/biomass_onlyvariance_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Biomass with mass flux
dummy = sk.scatterbioflux (fastbc, slowbc, equalbc)
picname = "X:/Saturated_flow/Anaerobic/biomass_flux_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Mass flux wtih varxanis
dummy = sk.boxVAflux(fastc, slowc, equalc, gvarnames)
picname = "X:/Saturated_flow/Anaerobic/flux_VA_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Biomass wtih varxanis
dummy = sk.boxVAbio (fastbc, slowbc, equalbc)
picname = "X:/Saturated_flow/Anaerobic/biomass_VA_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Flux wtih residence time
dummy = sk.boxrestime_flux(fastc, slowc, equalc, gvarnames)
picname = "X:/Saturated_flow/Anaerobic/flux_restime_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
#Flux wtih residence time
dummy = sk.boxrestime_biomass(fastbc, slowbc, equalbc)
picname = "X:/Saturated_flow/Anaerobic/biomass_restime_"+datetime.now().strftime("%d%m%Y%H%M")+".png"
dummy.savefig(picname)
    
# plot scatter plot for mean values
#Var0.1smeans = {'DOC': meanresultsarr[np.where(meanresultsarr[:,4]==0 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0], 'Ammonium': meanresultsarr[np.where(meanresultsarr[:,2]==0 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0], 'Nitrate': meanresultsarr[np.where(meanresultsarr[:,4]==3 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0]}

#Var0.1smeans = {'DOC': meanresultsarr[np.where(meanresultsarr[:,4]==0 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0], 'Ammonium': meanresultsarr[np.where(meanresultsarr[:,2]==0 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0], 'Nitrate': meanresultsarr[np.where(meanresultsarr[:,4]==3 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0]}
#Var0.1smeans = {'DOC': meanresultsarr[np.where(meanresultsarr[:,4]==0),:][np.where(meanresultsarr[:,2]==0.1),:][np.where(meanresultsarr[:,3]==2),0]}, 'Ammonium': meanresultsarr[np.where(meanresultsarr[:,2]==0 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0], 'Nitrate': meanresultsarr[np.where(meanresultsarr[:,4]==3 & (meanresultsarr[:,2]==0.1) & (meanresultsarr[:,3]==2)),0]}

#    varmeans = {'mhm_001': harmean_D_out_mhm_001, 'whitenoise_001': harmean_D_out_wn_001, 'mhm_00001': harmean_D_out_mhm_00001, 'whitenoise_00001': harmean_D_out_wn_00001}
#arimeans = {'mhm_001': arimean_D_out_mhm_001, 'whitenoise_001': arimean_D_out_wn_001, 'mhm_00001': arimean_D_out_mhm_00001,  'whitenoise_00001': arimean_D_out_wn_00001}
#names_species = list(speciesmeans.keys())
#values_species = list(speciesmeans.values())
#names_har = list(harmeans.keys())
#values_har = list(harmeans.values())
#names_ari = list(arimeans.keys())
#values_ari = list(arimeans.values())
#plt.figure(figsize=(9,9))
#plt.scatter(names_species, values_species, label = "speciesmeans", color="blue")
#plt.scatter(names_har, values_har, label = "harmean", color="orange")
#plt.scatter(names_ari, values_ari, label = "arimean", color="red")
#plt.yscale('log')
#plt.legend()
#plt.savefig("/Users/houben/Desktop/baseflow_sa/in_vs_out_means.png", dpi=300)

#Heatmap
#name = sk.heatmap(df, vars)    
#picname = di+fpre+str(Trial[j])+"_chematss.png"
#name.get_figure().savefig(picname)

#Scatter plot
#import plotly.express as px
#from plotly.offline import plot
#df = pd.DataFrame(meanresultsarr, columns = ['Removal_Rate', 'Removal_Rate_change', 'Variance', 'Anisotropy', 'Chemical'])
#mult5 = [0,5,10,15,20,25,30,35,40,45]
#for i in mult5:
#    df['Chemical'][i] = "DOC"
#    df['Chemical'][i+1] = "DO"
#    df['Chemical'][i+2] = "Amm"
#    df['Chemical'][i+3] = "Nitrate"
#    df['Chemical'][i+4] = "Sulphate"

#df['Removal_Rate'] = abs(df['Removal_Rate']*100)
#df['Removal_Rate_change'] = df['Removal_Rate_change']*100 

#df = df[df['Chemical']!='Sulphate']
#df = df[df['Chemical']!='DO']

#marker_size = (5, 10, 20, 30)
#marker_shape = ['.', 'o', '^', 'D']
#marker_color = "Blues"
#fig = px.scatter(df, x="Chemical", y="Removal_Rate_change", size="Removal_Rate", color = "Variance", symbol = "Anisotropy")
#plot(fig)