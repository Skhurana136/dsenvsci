# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:07:39 2019

@author: khurana
"""

import numpy  as np
import csv
import Pythonfunctions_SK as sk
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime

d = r"Z:/Saturated_flow/"
Tforfpre = ['Steady_state/EqualAR/','Transient/EqualAR_0.2/', 'Transient/EqualAR_0.4/']
fpre = 'NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]
fsuf  = r"/"

filename = 'ratesAtFinish.dat'

#setup what we really want to investigate
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis
Trial = ['H']
Het = [0]
Anis = [1]
#Trial = []
#Het = []
#Anis = []
for k in range(len(masterHet)):
    if (masterHet[k]>1):
        Trial.append(masterTrial[k])
        Het.append(masterHet[k])
        Anis.append(masterAnis[k])
listofcolumns = [6,12]
for i in range(67):
    listofcolumns.append(i+18)

ratenames = ["x_m", "Y","Z",
              "Fixedaeroresp", "Mobaeroresp","Fixedaerogwth", "Mobileaerogwth",
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

gratenames = ["Aerobic", "Nitrate", "Sulphate", "Ammonia"]

respindx = np.array([3, 19, 35, 51])
summ = np.zeros([len(Trial)*len(Tforfpre),len(respindx)])
act = np.zeros([len(Trial)*len(Tforfpre),1581, len(respindx)])
#Plot Damkohler numbers in the domain
for k in range(len(Tforfpre)):
    for j in range(len(Trial)):    
        di =  d+str(Tforfpre[k])+fpre+str(Trial[j])+fsuf
        print (str(Trial[j]))
        fwithd = di+filename
        M = np.loadtxt(fwithd, dtype = float, delimiter = ' ', usecols = 16+respindx)
        for i in range(len(respindx)):
            idx = k*len(Trial) + j
            summ[idx,i] = sum(M[:,i])
            act[idx,:,i] = M[:,i]

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