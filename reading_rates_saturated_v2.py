# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:05:45 2019

@author: khurana
"""

import numpy  as np
import csv
import Pythonfunctions_SK as sk
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.colors import LogNorm
import math
from matplotlib.ticker import FuncFormatter

Reg = "Equal"
d = r"Z:/Saturated_flow/diffusion_transient/"+Reg+"AR_0/"
fpre = 'NS-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]
fsuf  = r"/"

filename = 'ratesAtFinish.dat'

#setup what we really want to investigate
#Default:
#setup what we really want to investigate
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

#Variations:
#Trial = masterTrial[masterTrial.index(37):]
#Het = masterHet[masterTrial.index(37):]
#Anis = masterAnis[masterTrial.index(37):]

#Trial = ['H']
#Het = [0]
#Anis = [1]
#for k in range(len(masterHet)):
#    if (masterHet[k]>4):
#        Trial.append(masterTrial[k])
#        Het.append(masterHet[k])
#        Anis.append(masterAnis[k])

listofcolumns = [6,12]
for i in range(67):
    listofcolumns.append(i+18)

ratenames = ["x_m", "Y","Z",
              "Fixederoresp", "Mobaeroresp","Fixedaerogwth", "Mobileaerogwth",
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
              "Hydrolysis","POMgen"]

gratenames = ["Aerobic", "Nitrate", "Ammonia"]
respindx = np.array([ratenames.index("Fixederoresp"), ratenames.index("Mobaeroresp"),ratenames.index("Fixednitraresp"),ratenames.index("Mobnitraresp"), ratenames.index("Fixedsulpharesp"),ratenames.index("Mobsulpharesp"),ratenames.index("Fixedammresp"),ratenames.index("Mobammresp"),ratenames.index("Hydrolysis"),
                     ratenames.index("Fixedaerogwth"),ratenames.index("Mobileaerogwth"),ratenames.index("Fixednitragwth"),ratenames.index("Mobnitragwth"),ratenames.index("Fixedsulphagwth"),ratenames.index("Mobsulphagwth"),ratenames.index("Fixedammgwth"),ratenames.index("Mobammgwth")])
mobindx = np.array([ratenames.index("Fixedactaerodett"), ratenames.index("Fixedinaerodett"),ratenames.index("Fixedactammdett"),ratenames.index("Fixedinammdett"), ratenames.index("Fixedactnitradett"), ratenames.index("Fixedinnitradett")])
summ = np.zeros([len(Trial),len(respindx)])
act = np.zeros([len(Trial),1581, len(respindx)])
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
for j in range(len(Trial)):    
    di =  d+fpre+str(Trial[j])+fsuf
    print (str(Trial[j]))
    fwithd = di+filename
    print(fwithd)
    M = np.loadtxt(fwithd, dtype = float, delimiter = ' ', usecols = 16 + respindx)
    for i in range(len(respindx)):
        summ[j,i] = sum(M[:,i])
        act[j,:,i] = M[:,i]
    Numberofrates = np.shape(M)[1]
    ratedf = np.ndarray([Numberofrates, 51,31])
    ratedf2 = np.ndarray([Numberofrates, 51,31])
    Da = np.ndarray([Numberofrates, 51,31])
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    for i in range(Numberofrates):
            for k in range(51):
                for l in range(31):             
                    ratedf2[i,k,l] = M[(k)*31+l, i]
    for k in range(51):
        ratedf[:,50-k,:] = ratedf2[:,k,:]
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
    #rate[j,:] = (((ratedf[:,yout,xleft]*satoledg + ratedf[:,yout,xright]*satoredg)*vedge*vedge  + (np.sum(ratedf[:,yout,xleft+1:xright]*satoelem, axis = -1) + np.sum(ratedf[:,yin+1:yout,xleft]*satlelem, axis = -1) + np.sum(ratedf[:,yin+1:yout,xright]*satrelem, axis = -1))*vedge*velem) + np.sum(np.sum(ratedf[:,yin+1:yout,xleft+1:xright]*satelem, axis = -1), axis = -1)*velem*velem)*por
    rate[j,:] = (((np.sum(ratedf[:,yin+1:yout,xleft]*satlelem, axis = -1) + np.sum(ratedf[:,yin+1:yout,xright]*satrelem, axis = -1))*vedge*velem) + np.sum(np.sum(ratedf[:,yin+1:yout,xleft+1:xright]*satelem, axis = -1), axis = -1)*velem*velem)*por

f = "X:/massbalance_file2_"+Reg+".csv"
csvfile= open(f, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Fixedaeroresp_umolesperday","Mobaeroresp_umolesperday","Fixednitraresp_umolesperday","Mobnitraresp_umolesperday",
                 "Fixedsulpharesp_umolesperday","Mobsulpharesp_umolesperday","Fixedammresp_umolesperday","Mobammresp_umolesperday", "Hydrolysis_umolesperday",
                 "Fixedaerogrowth_umolesperday","Mobaerogrowth_umolesperday","Fixednitragrowth_umolesperday","Mobnitragrowth_umolesperday",
                 "Fixedsulphagrowth_umolesperday","Mobsulphagrowth_umolesperday","Fixedammgrowth_umolesperday","Mobammgrowth_umolesperday"])
for j in range(len(Trial)):
        writer.writerow([j, Trial[j], rate[j,0],rate[j,1],rate[j,2],rate[j,3],rate[j,4],rate[j,5],rate[j,6], rate[j,7], rate[j,8],rate[j,9],rate[j,10],rate[j,11],rate[j,12],rate[j,13],rate[j,14],rate[j,15],rate[j,16]])
csvfile.close()
    
    for resp in range(Numberofrates):
        Da[resp,:,:] = ratedf[resp,:,:]/abs(df[2,np.shape(df)[1]-1,:,:])
    df2 = np.load (di+fpre+str(Trial[j])+"_df.npy")
    vel = df2[2,np.shape(df)[1]-1,:,:]*-1
    fig, axes = plt.subplots(nrows = 1, ncols = len(gratenames) + 1, figsize=(26,10), sharex = True, sharey = True)
    plt.suptitle ("Variance "+str(Het[j])+" : Anisotropy "+str(Anis[j]), ha='center', va='center', fontsize = 30)
    for i in range(len(gratenames)):
        log_norm = LogNorm(vmin=Da[i,:,:].min().min(), vmax=Da[i,:,:].max().max())
        sns.set(font_scale = 2)
        cbar_ticks = [math.pow(10, x) for x in range(int(math.floor(math.log10(Da[i,:,:].min().min()))), 1+int(math.ceil(math.log10(Da[i,:,:].max().max()))))]
        name = sns.heatmap(Da[i,:,:], square = False,
                           norm = log_norm, cmap = "YlGnBu", 
                           cbar_kws = {"ticks":cbar_ticks}, vmin = Da[i,:,:].min().min(), vmax = Da[i,:,:].max().max(),
                           xticklabels=False, yticklabels=False, ax=axes.flat[i])
        axes.flat[i].set_title(gratenames[i], fontsize = 30)
    v = sns.heatmap(vel, square = False, norm = LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()), cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=axes.flat[3],
                            cbar_kws={'format': FuncFormatter(fmtvel), 'ticks':[math.pow(10, x) for x in range(int(math.floor(math.log10(np.min(abs(vel))))), 1+int(math.ceil(math.log10(np.max(abs(vel))))))]},
                           vmin = abs(vel).min().min(), vmax = abs(vel).max().max())
    axes.flat[3].set_title("Velocity", fontsize = 0)
    picname = d+fpre+str(Trial[j])+ "_" + "heatmap" + datetime.now().strftime("%d%m%Y%H%M")+"_LOG_Da.png"
    pictosave = name.get_figure()
    pictosave.savefig(picname)

#calculating node wise difference from homoegeneous scenario for the activities/respiration rates
actarr = sk.RateConverttomarr(act)
actarrscaled = np.zeros(np.shape(actarr))
for j in range(len(Trial)):
    for k in range(np.shape(actarrscaled)[0]):
        minrate = np.min(actarr[k,j,:,:])
        maxrate = np.max(actarr[k,j,:,:])
        actarrscaled[k,j,:,:] = (actarr[k,j,:,:] - minrate)/(maxrate- minrate)

delresprate = np.zeros([len(Trial),1581, len(respindx)])
delsummresprate = np.zeros([len(Trial),len(respindx)])
for j in range(len(Trial)):
    delresprate[j,:,:] = (act[0,:,:] - act[j,:,:])/act[0,:,:]
    delsummresprate[j,:] = (summ[0,:] - summ[j,:])/summ[0,:]
#    sk.heatmapratedist(delresprate[j,:,:], respindx)
delrespratearr = sk.RateConverttomarr(delresprate)
delrespratearrlog = np.zeros(np.shape(delrespratearr))
delrespratearrlog = np.zeros([4,49,51,31])
for j in range(len(Trial)):
    for k in range(np.shape(delrespratearr)[0]):
        for l in range(np.shape(delrespratearr)[2]):
            for m in range(np.shape(delrespratearr)[3]):
                if delrespratearr[k,j,l,m]==0:
                    pass
                else:
                    delrespratearrlog[k,j,l,m] = np.log(delrespratearr[k,j,l,m])

delrespratearrscaled = np.zeros([4,49,51,31])
for j in range(len(Trial)-1):
    for k in range(np.shape(delrespratearr)[0]):
        minrate = np.min(delrespratearr[k,j+1,:,:])
        maxrate = np.max(delrespratearr[k,j+1,:,:])
        delrespratearrscaled[k,j,:,:] = (delrespratearr[k,j,:,:] - minrate)/(maxrate- minrate)

intrespindx = np.array([3, 19, 51])
nrows = 1
ncols = 3
figsize = [10,5]
for j in range(len(Trial)):
    print(str(Trial[j]))
    name = sk.heatmapratedistLOG(delrespratearr[:,j,:,:], intrespindx, gratenames, nrows, ncols, figsize)
    picname = d+fpre+str(Trial[j])+datetime.now().strftime("%d%m%Y%H%M")+"_delrespratesatss_log.png"
    name.get_figure().savefig(picname)
    
#Writing the results in a csv file     
fname = d+fpre+str(Trial[0])+"_"+str(Trial[-1])+"_"+"Resp_"+ datetime.now().strftime("%d%m%Y%H%M")+".csv"
csvfile= open(fname, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Het", "Anis", "Rate_type", "Total_rate", "Change_in_rate"])
for j in range(len(Trial)):
    for k in range(len(respindx)):
        writer.writerow([j, Trial[j], Het[j], Anis[j], gratenames[k], summ[j,k] ,delsummresprate[j,k]])
csvfile.close()
print ("File written")

#Investigating attachment/detachment
respindx = np.array([ratenames.index("Fixedactaerodett"), ratenames.index("Fixedactammdett"),ratenames.index("Fixedactnitradett"),ratenames.index("Fixedinaerodett"),ratenames.index("Fixedinammdett"),  ratenames.index("Fixedinnitradett")])
gratenames = ["AA", "AM", "AN", "IA", "IM", "IN"]
summ = np.zeros([len(Trial),len(respindx)])
act = np.zeros([len(Trial),1581, len(respindx)])

respindx = np.array([ratenames.index("Mobactaeroattach"), ratenames.index("Mobactammattach"),ratenames.index("Mobactnitraattach"),ratenames.index("Mobinaeroattach"),ratenames.index("Mobinammattach"),  ratenames.index("Mobinnitraattach")])

intsce = ['H', 84, 63]
for j in range(len(intsce)):    
    di =  d+fpre+str(intsce[j])+fsuf
    print (str(intsce[j]))
    fwithd = di+filename
    M = np.loadtxt(fwithd, dtype = float, delimiter = ' ', usecols = 16 + respindx)
    for i in range(len(respindx)):
        summ[j,i] = sum(M[:,i])
        act[j,:,i] = M[:,i]
    Numberofrates = np.shape(M)[1]
    ratedf = np.ndarray([Numberofrates, 51,31])
    ratedf2 = np.ndarray([Numberofrates, 51,31])
    Da = np.ndarray([Numberofrates, 51,31])
    fmtvel = lambda x, pos: '{:1.1e}'.format(x)
    for i in range(Numberofrates):
            for k in range(51):
                for l in range(31):
                    ratedf2[i,k,l] = M[(k)*31+l, i]
    for k in range(51):
        ratedf[:,50-k,:] = ratedf2[:,k,:]
#    df2 = np.load (di+fpre+str(Trial[j])+"_df.npy")
#    vel = df2[2,-1,:,:]*-1
    fig, axes = plt.subplots(nrows = 2, ncols = 3, figsize=(26,10), sharex = True, sharey = True)
    plt.suptitle ("Variance "+str(Het[Trial.index(intsce[j])])+" : Anisotropy "+str(Anis[Trial.index(intsce[j])]), ha='center', va='center', fontsize = 30)
    for i in range(len(gratenames)):
#        log_norm = LogNorm(vmin=Da[i,:,:].min().min(), vmax=Da[i,:,:].max().max())
#        sns.set(font_scale = 2)
#        cbar_ticks = [math.pow(10, x) for x in range(int(math.floor(math.log10(Da[i,:,:].min().min()))), 1+int(math.ceil(math.log10(Da[i,:,:].max().max()))))]
        name = sns.heatmap(ratedf2[i,:,:], square = False,
                           #norm = log_norm, 
                           cmap = "YlGnBu", 
                          # cbar_kws = {"ticks":cbar_ticks}, vmin = Da[i,:,:].min().min(), vmax = Da[i,:,:].max().max(),
                           xticklabels=False, yticklabels=False, ax=axes.flat[i])
        axes.flat[i].set_title(gratenames[i], fontsize = 30)
#    v = sns.heatmap(vel, square = False, norm = LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()), cmap = "YlGnBu", xticklabels=False, yticklabels=False, ax=axes.flat[3],
#                            cbar_kws={'format': FuncFormatter(fmtvel), 'ticks':[math.pow(10, x) for x in range(int(math.floor(math.log10(np.min(abs(vel))))), 1+int(math.ceil(math.log10(np.max(abs(vel))))))]},
#                           vmin = abs(vel).min().min(), vmax = abs(vel).max().max())
#    axes.flat[3].set_title("Velocity", fontsize = 0)
    picname = d+fpre+str(intsce[j])+ "_" + "heatmap_" + datetime.now().strftime("%d%m%Y%H%M")+"_attachment.png"
    pictosave = name.get_figure()
    pictosave.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0)
