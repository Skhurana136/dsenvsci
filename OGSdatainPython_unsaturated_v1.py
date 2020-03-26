# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy  as np
import csv
import Pythonfunctions_SK as sk

#Saturated flow regime
Reg = "Slow"
d = r"X:/Richards_flow/"+ Reg+ "AR_0/"
fpre = 'RF-A'
masterTrial = ['H',37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84]
masterHet = [0,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,10,10,10,0.1,0.1,0.1,1,1,1,5,5,5,10,10,10,5,5,5,5,5,5,5,5,5]
masterAnis = [1,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10,2,5,10]

fsuf  = r"/"
gw = 0

filename = 'model_domain_quad.tec'

#setup what we really want to investigate
#Default:
Trial = masterTrial
Het = masterHet
Anis = masterAnis

#Variations:
Trial = []
notlist = [45, 49, 51, 56, 58]
for i in masterTrial:
    if i not in notlist:
        Trial.append(i)
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
vars = [doc1, dox1, Amm1, nitra1, sulpha1, Bmo1, Bma1, Bmn1, Bimo1, Bima1, Bimn1]
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate", "Sulphate", "Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers",
                      "Inactive mobile Aerobes", "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers","Nitrogen", "TOC", ]
AFbiomassvars = [Bfo1, Bfa1, Bfn1]
AFbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AFbiomassvars = [Bfo1, Bfa1, Bfn1, Bmo1, Bma1, Bmn1, Bifo1, Bifa1, Bifn1, Bimo1, Bima1, Bimn1]
AFbiomassgvarnames = ["Active fixed Aerobes", "Active fixed Ammonia oxidizers", "Active fixed Nitrate reducers",
                      "Active mobile Aerobes", "Active mobile Ammonia oxidizers", "Active mobile Nitrate reducers",
                      "Inactive fixed Aerobes", "Inactive fixed Ammonia oxidizers", "Inactive fixed Nitrate reducers",
                      "Inactive mobile Aerobes", "Inactive mobile Ammonia oxidizers", "Inactive mobile Nitrate reducers"]
AMbiomassvars = [Bmo1, Bma1, Bmn1]
AMbiomassgvarnames = ["Mobile Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AFbiomassvars = [Bifo1, Bifa1, Bifn1]
AFbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AMbiomassvars = [Bimo1, Bima1, Bimn1]
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
    np.save(di+fpre+str(Trial[j])+'_df_2', df)
    #Test for correct orientation of the data
    for i in range(np.shape(df)[0]):
        print(Headers[i+3], np.mean(df[i,steps-1,0,:]))

mf = sk.calcmassflux (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
biomasssum = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
biomasssummob = sk.calcsum (Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, AMbiomassvars, AMbiomassgvarnames)
#Tracer studies       
tr1 = 8-gw
vely=5
velx=4
vars = [tr1]
steps = [1, 0.1, 0.02]
Regimes = ["Slow", "Equal", "Fast"]
f = "X:/Richards_flow/Steady_state/Tracer_studies/tracer_combined.csv"
csvfile= open(f, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Time", "fraction","Regime"])
idx = 1
for Reg,step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/"+Reg+ "AR/"
    fpre = "NS-A"
    df,massendtime, masstime, conctime, Velocity, head = sk.calcconcmasstime (Trial[0], Het[0], Anis[0], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    Time = np.where(np.round(conctime[:,yout, 1],3)>10)
    initial = step*Time[0][0]
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity, head = sk.calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        Time = np.where(np.round(conctime[:,yout, 1],3)>10)
        Time2 = np.where(np.round(df[-1, :, 50,:],3)>10)
        print (step*Time[0][0], step*Time2[0][0], initial, (step*Time[0][0])/initial)
        writer.writerow([idx,Trial[j], Het[j], Anis[j], "Tracer", step*Time[0][0], (step*Time[0][0])/initial, Reg])
        idx = idx + 1
csvfile.close()


for j in Trial:
    df = np.load(d+fpre+str(j)+fsuf+fpre+str(j)+"_df.npy")
    print (d+str(j), ": last time step: ", (np.shape(df)[1]-1)*5, 5475 - (np.shape(df)[1]-1)*5)
    for h in Headers[3:]:
        fname = "X:/Saturated_flow/changedkindox_transient/"+Reg+"AR_1/NS-A"+str(j)+fsuf+"model_MASS_TRANSPORT_"+h[2:-1]+"_primary_value.asc"
        csvfile= open(fname, "w")
        writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
        for y in range(51):
            for x in range(31):
                writer.writerow([df[Headers.index(h)-3, -1, yout - y, x], df[Headers.index(h)-3, -1, yout - y, x]])
        csvfile.close()
