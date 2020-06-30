# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 16:00:17 2020

@author: khurana
"""
import numpy  as np
import csv
import matplotlib.pyplot as plt
import plots.saturated_steady_state as sssp
import analyses.saturated_transient as sta

#Saturated flow regime
Reg = "Fast"
directory = r"Z:/Saturated_flow/diffusion_transient/"+Reg+"AR_0/"
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

#Constants
yin = 0
yout = 50
xleft = 0
xright = 30
tr1 = 8-gw
vely=5
velx=4
vars = [tr1]
steps = [500*0.005, 2*0.005, 2*0.0005]
Regimes = ["Slow","Equal", "Fast"]
f = "X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv"
csvfile= open(f, "w")
writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Time", "fraction","Regime"])
idx = 1
for Reg,step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/"+Reg+ "AR/"
    fpre = "NS-A"
    df,massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime (Trial[0], Het[0], Anis[0], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
    print (np.mean(df[vely-3,1:,:,:]))
    Time = np.where(np.round(conctime[:,yout, 0],3)>10)
    initial = step*Time[0][0]
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        print (np.mean(df[vely-3,1:,:,:]))
        Time = np.where(np.round(conctime[:,yout, 0],3)>10)
        Time2 = np.where(np.round(df[-1, :, 50,:],3)>10)
        s = step
        if (Reg == "Slow"):
            if (Trial[j]==54):
                s = 100*0.01
        print (s*Time[0][0], s*Time2[0][0], initial, (s*Time[0][0])/initial)
        writer.writerow([idx,Trial[j], Het[j], Anis[j], "Tracer", s*Time[0][0], (s*Time[0][0])/initial, Reg])
        idx = idx + 1
csvfile.close()
#plotting boxplots to see variance of breakthrough from homogeneous scenario
tracerplot = sssp.plot_tracer()
tracerplot.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.png", dpi = 300, pad_inches = 0.01)
tracerplot = sssp.plot_tracer()
tracerplot.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.png", dpi = 300, pad_inches = 0.01)

#plotting time series of tracer breakthrough in all flow regimes and trials:
ncols = 8
nrows = 6
for Reg, step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/"+Reg+ "AR/"
    fpre = "NS-A"
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = [15,10],sharex = True, sharey = True)
    plt.suptitle("Tracer breakthrough curve" + Reg + " flow regime")
    col = 0
    for j in range(len(Trial)):
        df,massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime (Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars)
        xindex = list(x*step for x in list(range(np.shape(conctime)[0])))
        axes.flat[j].plot(xindex, conctime[:,yout,0])
        axes.flat[j].set_xlabel('')
#        axes[colidx1][colidx2].set_xticklabels([])
        axes.flat[j].set_ylabel('')
        axes.flat[j].set_title(Trial[j])
        col = col+1
    for ax in axes[:,0]:
        ax.set_ylabel("Tracer (uM)")
    for ax in axes[5]:
        ax.set_xlabel("Time (days)")
    fig.subplots_adjust(left=0.15, top=0.9)
    fig.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughcurve_"+Reg+".png", dpi = 300)
    fig.savefig("X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughcurve_"+Reg+".pdf", dpi = 300)