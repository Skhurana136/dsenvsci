# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:07:39 2019

@author: khurana
"""

import numpy  as np
import csv
import pandas as pd
import data_reader.data_processing as proc
#import seaborn as sns
#import matplotlib.pyplot as plt
#from datetime import datetime

d = r"Z:/Saturated_flow/"
Tforfpre = ['Steady_state/EqualAR/','Transient/EqualAR_0.2/', 'Transient/EqualAR_0.4/']
fpre = 'NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]
fsuf  = r"/"
gw = 1
filename = 'ratesAtFinish.dat'

#setup what we really want to investigate
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

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
POM1 = 30 - gw
vely=5
velx=4
AFbiomassvars = [Bfo1, Bfa1, Bfn1, Bmo1, Bma1, Bmn1, Bifo1, Bifa1, Bifn1, Bimo1, Bima1, Bimn1]
AFbiomassgvarnames = ["Active fixed Aerobes", "Active fixed Ammonia oxidizers", "Active fixed Nitrate reducers",
                      "Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers",
                      "Inactive fixed Aerobes", "Inactive fixed Ammonia oxidizers", "Inactive fixed Nitrate reducers",
                      "Inactive mobile Aerobes", "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers"]

listofcolumns = [6,12]
listofcolumns = []
for i in range(67):
    listofcolumns.append(i+18)

ratenames = ["x_m", "Y","Z",
              "Fixedaeroresp", "Mobaeroresp","Fixedaerogwth", "Mobaerogwth",
              "Fixedactaerodett", "Fixedinaerodett", "Mobactaeroattach", "Mobinaeroattach",
              "FixeddeactlowDOX", "MobdeactlowDOX", "Fixedaeroreact", "Mobaeroreact",
              "Mortfixedactaero","Mortmobactaero", "Mortinfixedaero", "Mortinmobaero",
              "Fixednitraresp", "Mobnitraresp","Fixednitragwth", "Mobnitragwth",
              "Fixedactnitradett","Fixedinnitradett", "Mobactnitraattach", "Mobinnitraattach",
              "FixeddeactlowN", "MobdeactlowN", "Fixednitrareact", "Mobnitrareact",
              "Mortfixedactnitra", "Mortmobactnitra", "Mortinfixednitra", "Mortinmobnitra",
              "Fixedsulpharesp","Mobsulpharesp", "Fixedsulphagwth","Mobsulphagwth",
              "Fixedactsulphadett", "Fixedinsulphadett", "Mobactsulphaattach", "Mobinsulphaattach",
              "FixedDeactlowS", "MobDeactlowS","Fixedsulphareact", "Mobsulphareact",
              "Mortfixedactsulpha", "Mortmobactsulpha", "Mortinfixedsulpha", "Mortinmobsulpha",
              "Fixedammresp", "Mobammresp", "Fixedammgwth", "Mobammgwth",
              "Fixedactammdett", "Fixedinammdett", "Mobactammattach", "Mobinammattach",
              "FixedammdeactlowA", "MobammdeactlowA", "Fixedammreact", "Mobammreact",
              "Mortfixedactamm", "Mortmobactamm", "Mortinfixedamm", "Mortinmobamm",
              "Hydrolysis"]

respindx = np.array([ratenames.index("Fixedaerogwth"),ratenames.index("Mobaerogwth"),ratenames.index("Fixednitragwth"),ratenames.index("Mobnitragwth"),ratenames.index("Fixedsulphagwth"),ratenames.index("Mobsulphagwth"),ratenames.index("Fixedammgwth"),ratenames.index("Mobammgwth")])
respindx = np.array([ratenames.index("Fixedaerogwth"),ratenames.index("Fixedammgwth"),ratenames.index("Fixednitragwth")])
respindx = np.array([ratenames.index("Fixedaeroresp"),ratenames.index("Mobaeroresp"),ratenames.index("Fixednitraresp"),ratenames.index("Mobnitraresp"),ratenames.index("Fixedsulpharesp"),ratenames.index("Mobsulpharesp"),ratenames.index("Fixedammresp"),ratenames.index("Mobammresp")])
gratenames = ["Immobile aerobic respiration", "Mobile aerobic respiration", "Immobile nitrate respiration", "Mobile nitrate respiration","Immobile sulphate respiration", "Mobile sulphate respiration","Immobile ammonia respiration", "Mobile ammonia respiration"]

Regimes = ["Slow", "Equal", "Fast"]
summ = np.zeros([len(Trial)*len(Tforfpre)*len(Regimes)*len(respindx), 7])
act = np.zeros([len(Trial)*len(Tforfpre)*len(Regimes),1581, len(respindx)+3])
normsumm = np.zeros([len(Trial)*len(Tforfpre)*len(Regimes)*len(respindx),7])
#Plot Damkohler numbers in the domain

yin = 0
yout = 50
xleft = 0
xright = 30
vedge = 0.005
velem = 0.01
por = 0.2
gw = 1
rate = np.zeros([len(Trial),len(respindx)])

for Reg in Regimes:
    d = r"Z:/Saturated_flow/diffusion_transient/"+Reg+"AR_"
    Tforfpre = ['0/','1/','2/', '5/']
    print (Reg)
    for k in range(len(Tforfpre)):
        print(Tforfpre[k])
        for j in range(len(Trial)):    
            di =  d+str(Tforfpre[k])+fpre+str(Trial[j])+fsuf
            print (str(Trial[j]))
            fwithd = di+filename
            M = np.loadtxt(fwithd, dtype = float, delimiter = ' ', usecols = 16+respindx)
            df = np.load(di+fpre+str(Trial[j])+"_df.npy")
            for i in range(len(respindx)):
                idx = Regimes.index(Reg)*len(Tforfpre)*len(Trial)*len(respindx)+k*len(Trial)*len(respindx) + j*len(respindx) + i
                summ[idx,6] = sum(M[:,i])
#                act[idx,:,i] = M[:,i]
#                ratearr = sk.RateConverttomarr(M)
#                normsumm[idx,i] = np.mean(M[i,:,:]/df[AFbiomassvars[i]-3, -1, :,:])
                sumbiomass = np.sum(df[AFbiomassvars[i]-3, -1, :,:])
                normsumm[idx,6]= summ[idx,6]/sumbiomass
                if j==0:
                    summ[idx,0] = j
                else:
                    summ[idx,0] = Trial[j]
                summ[idx,3] = Regimes.index(Reg)
                summ[idx,4] = k
                summ[idx,5] = i
                summ[idx,1] = masterHet[masterTrial.index(Trial[j])]
                summ[idx,2] = masterAnis[masterTrial.index(Trial[j])]
                normsumm[idx, :6] = summ[idx,:6]

dfsum = pd.DataFrame(data = summ, columns = ["Trial", "Variance", "Anisotropy", "Regime", "Time_series", "Chem", "Sum_rate"])
dfall2 = proc.processdataframe(dfsum, gratenames)
data["Regime"] = data["Regime"].astype(int)
data["Regime"] = data["Regime"].replace([0,1,2], ["Slow","Medium", "Fast"])
#    data["Regime"] = data["Regime"].replace([0,1], ["Medium", "Fast"])
    data["Regime"] = data["Regime"].astype(str)
    data["Trial"] = data["Trial"].astype(int)
    data["Trial"] = data["Trial"].replace([0], "H")
    data["Trial"] = data["Trial"].astype(str)
    data["Chem"]= data["Chem"].astype(int)
    for k in range(len(variablenames)):
        data["Chem"] = data["Chem"].replace([k], variablenames[k])
    data["Chem"]= data["Chem"].astype(str)

filename = "X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv"
breakthrough = pd.read_csv(filename, delimiter = "\t")
breakthrough ['Heterogeneity'] = breakthrough ['Variance'] 
breakthrough['VA']=breakthrough['Heterogeneity']*breakthrough['Anisotropy']
breakthrough['%ofhomogeneous']=breakthrough['fraction']*100
breakthrough.loc[breakthrough.Regime=="Equal", 'Regime'] = "Medium"

dfall2 = pd.merge(data, breakthrough[['Trial', 'Regime', 'Time','fraction','%ofhomogeneous']], on = ['Trial', 'Regime']).rename(columns={'Time':'Breakthroughtime'})

dfall2.to_csv("Z:/Saturated_flow/diffusion_transient/total_rate_withbreakthrough.csv",sep='\t')
#calculating node wise difference from homoegeneous scenario for the activities/respiration rates
delresprate = np.zeros([len(Trial*len(Tforfpre)),1581, len(respindx)])
delsummresprate = np.zeros([len(Trial)*len(Tforfpre),len(respindx)])
delsummhetresprate = np.zeros([len(Trial)*len(Tforfpre),len(respindx)])
for k in range(len(Tforfpre)):
    for j in range(len(Trial)):
        idx = k*len(Trial) + j
        delresprate[idx,:,:] = (act[idx,:,:] - act[0,:,:])/act[0,:,:]
        delsummresprate[idx,:] = (summ[idx,:] - summ[0,:])/summ[0,:]
        delsummhetresprate[idx,:] = (summ[idx,:] - summ[j,:])/summ[j,:]
#    sk.heatmapratedist(delresprate[j,:,:], respindx)
delrespratearr = sk.RateConverttomarr(delresprate)

for j in range(len(Trial)):
    print(str(Trial[j]))
    name = sk.heatmapratedist(df[:,j,:,:], respindx, ratenames)
    picname = d+fpre+str(Trial[j])+datetime.now().strftime("%d%m%Y%H%M")+"_ratesatss.png"
    name.get_figure().savefig(picname)
    
#Writing the results in a csv file     
fname = d+fpre+str(Trial[0])+"_"+str(Trial[-1])+"_"+"Resp_"+ datetime.now().strftime("%d%m%Y%H%M")+"_Equal.csv"
csvfile= open(fname, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Het", "Anis", "Rate_type", "Total_rate", "Change_in_rate_H", "Change_in_rate_Time", "Time"])
for i in range(len(Tforfpre)):
    for j in range(len(Trial)):
        for k in range(len(respindx)):
            idx = i*len(Trial) + j
            writer.writerow([idx+k, Trial[j], Het[j], Anis[j], gratenames[k], summ[idx,k] ,delsummresprate[idx,k], delsummhetresprate[idx,k],i])
csvfile.close()
print ("File written")