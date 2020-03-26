# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""

import numpy  as np
import csv
import Pythonfunctions_SK as sk
from datetime import datetime
import matplotlib.pyplot as plt
import os

#Saturated flow regime
d = ""
Reg = "Equal"
Tforfpre = ['Z:/Saturated_flow/Anaerobic/'+Reg+'AR_changedkindox','Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_1', 'Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_2','Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_5']

fpre = '/NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]
#Tforfpre = ["H", 45, 54, 63, 75]
#Trial = ['','-T7']
Cov = [1, 2, 5]
fsuf  = r"/"
gw = 1

filename = 'model_domain_quad.tec'

#Scan for files
for j in Tforfpre:
    print (j)
    for k in Trial:
            file_name = fpre[1:]+str(k)+"_df.npy"
            #np.load(d+j+fpre+str(k)+fsuf+fpre+str(k)+"_df.npy")
            cur_dir = d+j+fpre+str(k)+fsuf # Dir from where search starts can be replaced with any path
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
notlist = [55, 56, 57, 58, 59, 60, 61, 62, 63, 64]
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
vars = [doc1, dox1, Amm1, nitra1, sulpha1]
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate", "Sulphate", "Nitrogen"]
AFbiomassvars = [Bfo1, Bfa1, Bfn1]
AFbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AMbiomassvars = [Bmo1, Bma1, Bmn1]
AMbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]

for Reg in ["Slow","Equal", "Fast"]:
    Tforfpre = ['Z:/Saturated_flow/Anaerobic/'+Reg+'AR_changedkindox','Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_1', 'Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_2','Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_5']
    mft, calcsum = sk.calcmft_temp(Tforfpre, Trial, vars, gvarnames, AFbiomassvars, AFbiomassgvarnames, d, fpre, fsuf, Het, Anis, gw)
    mftbuffer, calcsummob = sk.calcmft_temp(Tforfpre, Trial, vars, gvarnames, AMbiomassvars, AMbiomassgvarnames, d, fpre, fsuf, Het, Anis, gw)
    #Writing the results in a csv file     
    fnamemf = "X:/Saturated_flow/Anaerobic_"+Reg+"_Transient_"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"MF_"+str(gw) + ".csv"
    csvfile= open(fnamemf, "w")
    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Inlet_total_mass_flux", "Outlet_mass_flux", "delmassflux", "del2massflux", "del2massflux_Time","Time", "Regime"])
    idx = 0
    for l in range(len(Tforfpre)):
        for j in range(len(Trial)):
            for i in range(len(gvarnames)):
                idx = (l*len(Trial)+j)*len(gvarnames)+i
                if (Reg == "Equal"):
                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,8],mft[idx,9],"Medium"])
                else:
                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,8],mft[idx,9],Reg])
    csvfile.close()
        
    fnamebiosum = "X:/Saturated_flow/Anaerobic_"+Reg+"_Transient_Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"biomass_"+str(gw) +".csv"
    csvfile= open(fnamebiosum, "w")
    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles", "Change_umoles_Time", "Time", "Regime"])
    idxb = 0
    for l in range(len(Tforfpre)):
        for j in range(len(Trial)):
            for i in AFbiomassgvarnames:
                idx = (l*len(Trial) + j)*len(AFbiomassgvarnames)+AFbiomassgvarnames.index(i)
                if (Reg == "Equal"):
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7], "Medium"])
                else:
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7], Reg])
    csvfile.close()
    
    fnamebiosummob =  "X:/Saturated_flow/Anaerobic_"+Reg+"_Transient_Anaerobic"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"_mobile_biomass_"+str(gw) +".csv"
    csvfile= open(fnamebiosummob, "w")
    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_mobile_biomass", "Change_umoles", "Change_umoles_Time", "Time", "Regime"])
    idxb = 0
    for l in range(len(Tforfpre)):
        for j in range(len(Trial)):
            for i in AMbiomassgvarnames:
                idx = (l*len(Trial) + j)*len(AMbiomassgvarnames)+AMbiomassgvarnames.index(i)
                if (Reg == "Equal"):
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsummob[idx,4], calcsummob[idx,5], calcsummob[idx,6],calcsummob[idx,7], "Medium"])
                else:
                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsummob[idx,4], calcsummob[idx,5], calcsummob[idx,6],calcsummob[idx,7], Reg])
    csvfile.close()
    
    print ("Files written, processing as dataframes ...")
    if Reg == "Fast":
#        Reg = "Fast"
#        fnamemf = "Fast_Transient_H_59_MF_1_time.csv"
#        fnamebiosum = "Fast_Transient_AnaerobicH_59_biomass_1_time.csv"
        fastbc = sk.processchembiomassfiles(fnamebiosum, Reg)
        fastmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        fastc = sk.processchemfiles(fnamemf, Reg)
    elif Reg == "Slow":
#        Reg = "Slow"
#        fnamemf = "X:/Saturated_flow/Anaerobic_Slow_Transient_H_59_MF_1.csv"
#        fnamebiosum = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59_biomass_1.csv"
#        fnamebiosummob = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59__mobile_biomass_1.csv"
        slowbc = sk.processchembiomassfiles(fnamebiosum, Reg)
        slowmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        slowc = sk.processchemfiles(fnamemf, Reg)
    elif Reg == "Equal":
#        Reg = "Equal"
#        fnamemf = "Equal_Transient_H_59_MF_1_time.csv"
#        fnamebiosum = "Equal_Transient_AnaerobicH_59_biomass_1_time.csv"
        equalbc = sk.processchembiomassfiles(fnamebiosum, "Medium")
        equalmbc = sk.processchembiomassfiles(fnamebiosummob, "Medium")
        equalc = sk.processchemfiles(fnamemf, "Medium")
    print ("Dataframes processed, move to next series")


#Aggregated results
#Importing data

#Mass flux and biomass with average residence time
dummy = sk.scatterrestime_biomass_temp(fastbc, slowbc, equalbc)
dummy.savefig("X:/Saturated_flow/Anaerobic/scatterrestime_biomass_temp.png")
dummy = sk.scatterrestime_biomass_temp(fastmbc, slowmbc, equalmbc)
dummy.savefig("X:/Saturated_flow/Anaerobic/scatterrestime_mobile_biomass_temp.png")
dummy = sk.scatterrestime_flux_temp(fastc, equalc, slowc, ["DOC", "DO", "Nitrogen"])
dummy.savefig("X:/Saturated_flow/Anaerobic/scatterrestime_flux_temp.png")

Regimes = ["Slow", "Equal", "Fast"]
intsce = ["H", 76, 73, 80, 84, 63]
for i in intsce:
    initseries = [500, 430, 600] 
    lastseries = [700, 630, 800]
    sk.generate_timeseries(Regimes, initseries, lastseries, Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, d, fpre, vars, gvarnames, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)

for j in Tforfpre[1:]:
    #Concentration profiles
    intsce = ["H", 76, 73, 80, 84, 63]
#    for chem in ["DO", "DOC"]:
#        dum = sk.shiftoxic_anoxic_temp (chem, Trial, intsce, d, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars)
#        picname = "X:/Saturated_flow/Anaerobic/chemprofile_temp_"+chem+str(Tforfpre.index(j))+"_"+Reg+"_.png"
#        dum.savefig(picname)

    dummy = sk.concprofilewithtime_temp (Trial, intsce, d, j, gvarnames, Het, Anis, gw, fpre, fsuf, yin, yout, xleft, xright, vars)
    picname = "X:/Saturated_flow/Anaerobic/chemswithvel_temp_"+str(Tforfpre.index(j))+"_"+Reg+"_.png"
    dummy.savefig(picname)
    
    dfall = slowc
    l = []
    for i in range(len(dfall)):
        if (dfall['Variance'][i] == 0.1):
            l.append(str(dfall['Variance'][i]) + ":" + str(dfall['Anisotropy'][i]))
        else:
            l.append(str(int(dfall['Variance'][i])) + ":" + str(dfall['Anisotropy'][i]))
        
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance','Anisotropy'])
    bth = sk.tracerstudies()
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    Regimes = ["Slow"]
    dfall2["del2massflux_Time%"] = dfall2["del2massflux_Time"]*100
    Chems = ["DOC", "DO", "Ammonium", "Nitrate", "Nitrogen"]
    colseries = ["Greens"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [12,12])
    plt.suptitle("Change in removal of carbon and nitrogen in transient conditions", fontsize = 25)
    for i in Regimes:
        dfall0 = dfall2[dfall2['Time'] != 0]
        df = dfall2[dfall2['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chem']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "del2massflux_Time%", hue = "Time", palette = colseries[Regimes.index(i)], data = dfc, ax = axes.flat[Chems.index(k)])
            axes.flat[Chems.index(k)].set_ylabel(k, fontsize = 15)
            axes.flat[Chems.index(k)].legend_.remove()
            axes.flat[Chems.index(k)].tick_params(labelsize = 15)
            if (k != "Nitrogen"):
                axes.flat[Chems.index(k)].set_xlabel('')
                axes.flat[Chems.index(k)].set_xticklabels([])
            else:
                axes.flat[Chems.index(k)].set_xlabel("Variance : Anisotropy", fontsize = 15)
            col = col+1
    fig.savefig("X:/Saturated_flow/Anaerobic/flux_varianceandtime_"+datetime.now().strftime("%d%m%Y%H%M")+"_time.png")
    
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
            axes.flat[Chems.index(k)].scatter(x = "%ofhomogeneous", y = "del2massflux_Time%", c = "Time", cmap = "Greens", data = dfc)
            axes.flat[Chems.index(k)].set_xlabel('')
            axes.flat[Chems.index(k)].set_ylabel('')
            axes.flat[Chems.index(k)].set_xticks(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10.0))
            axes.flat[Chems.index(k)].set_xticklabels([])
            col = col+1
    for ax, col in zip(axes, Regimes):
        ax.set_title(col+" flow", ha = 'center', fontsize = 15)
    for ax, row in zip(axes, Chems):
        ax.annotate(row, xy = (0, 0.5), xytext = (-ax.yaxis.labelpad + 255,0), xycoords = 'axes fraction', textcoords='offset points', size = 'large', ha = 'left', va = 'center', rotation = 'vertical', fontsize = 15)
    for ax in axes:
        ax.set_xticklabels(np.arange(round(min(np.unique(dfall2['%ofhomogeneous'])),-1), round(max(np.unique(dfall2['%ofhomogeneous'])),-1), 10), size = 15)
    plt.figtext(0.5, 0.08, 'Relative difference in breakthrough time (%)', ha='center', va='center', fontsize = 20)
    plt.figtext(0.08, 0.5, 'Relative difference (%)', ha='center', va='center', rotation='vertical', fontsize = 20)   
   
    
    dfall = slowbc
    l = []
    for i in range(len(dfall)):
        if (dfall['Variance_b'][i] == 0.1):
            l.append(str(dfall['Variance_b'][i]) + ":" + str(dfall['Anisotropy_b'][i]))
        else:
            l.append(str(int(dfall['Variance_b'][i])) + ":" + str(dfall['Anisotropy_b'][i]))
    dfall['Xlabels'] = l
    dfall = dfall.sort_values(by=['Variance_b','Anisotropy_b'])
    dfall["delbiomass_Time%"]=dfall["Change_umoles_Time"]*100
    bth = sk.tracerstudies()
    dfall2 = pd.merge(dfall, bth[['Trial', 'Regime', 'Firsthit', '%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Firsthit':'Residencetime'})
    
    Regimes = ["Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Greens"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [12,12])
    plt.suptitle("Change in biomass in transient conditions", fontsize = 25)
    for i in Regimes:
        dfall0 = dfall[dfall['Time_b']!=0]
        df = dfall[dfall['Regime']==i]
        col = 0
        for k in Chems:
            dfc = df[df['Chemb']==k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(x = "Xlabels", y = "delbiomass%", hue = "Time_b", palette = colseries[Regimes.index(i)], data = dfc, ax = axes.flat[Chems.index(k)])
            axes.flat[Chems.index(k)].set_ylabel(k, fontsize = 15)
            axes.flat[Chems.index(k)].legend_.remove()
            axes.flat[Chems.index(k)].tick_params(labelsize = 15)
            if (k != "Nitrate reducers"):
                axes.flat[Chems.index(k)].set_xlabel('')
                axes.flat[Chems.index(k)].set_xticklabels([])
            else:
                axes.flat[Chems.index(k)].set_xlabel("Variance : Anisotropy", fontsize = 15)
            col = col+1
    fig.savefig("X:/Saturated_flow/Anaerobic/biomass_varianceandtime_"+datetime.now().strftime("%d%m%Y%H%M")+"_time2.png")
    
#Calculating average velocities
for Reg in ["Slow","Equal", "Fast"]:
    print (Reg)
    Tforfpre = ['X:/Saturated_flow/Anaerobic/'+Reg+'AR_changedkindox',
                'Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_1',
                'Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_2',
                'Z:/Saturated_flow/changedkindox_Transient/'+Reg+'AR_5']
    for j in Tforfpre:
        print(j)
        if (Tforfpre.index(j)==0):
            pass
        else:
            intsce = Trial
            for i in intsce:
                newd = d + j
                df,massendtime, masstime, conctime, Velocity, head = sk.calcconcmasstime (Trial[Trial.index(i)], Het[Trial.index(i)], Anis[Trial.index(i)], gw, newd, fpre, fsuf, yin, yout, xleft, xright, vars)
                Velocity = df[2,:,:,:]  
                print(i, ": ",np.mean(Velocity))