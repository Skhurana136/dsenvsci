# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 16:00:17 2020

@author: khurana
"""
import numpy as np
import csv
import matplotlib.pyplot as plt
import plots.saturated_steady_state as sssp
import analyses.saturated_transient as sta
import data_reader.data_processing as proc

# Saturated flow regime
Reg = "Fast"
directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
fpre = "NS-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

yin = 0
xleft = 0
xright = 30
trialist = ['H','37', '38', '39', '40', '41', '42', '43', '44', '45']
Trial = list(t for t,values in scdict.items() if t in trialist)
Het = list(values['Het'] for t,values in scdict.items() if t in Trial)
Anis = list(values['Anis'] for t,values in scdict.items() if t in Trial)

fsuf = r"/"
gw = 1

filename = "model_domain_quad.tec"

# plotting time series of tracer breakthrough in all flow regimes and trials:
ncols = 2
nrows = 5
Regimes = ["Slow"]
steps = [200 * 0.01]
yin = 0
xleft = 0
xright = 30
trialist = ['H','37', '38', '39', '40', '41', '42', '43', '44', '45']
Trial = list(t for t,values in scdict.items() if t in trialist)

yout = 125
tr1 = 6 - gw
vely = 5
velx = 4
vars = [tr1]
gvarnames = ["Tracer"]
for Reg, step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/Big_" + Reg + "AR/"
    fpre = "NS-A"
    fig, axes = plt.subplots(
        ncols=ncols, nrows=nrows, figsize=[15, 10], sharex=True#, sharey=True
    )
    plt.suptitle("Tracer breakthrough curve" + Reg + " flow regime")
    col = 0
    for j in Trial:
        df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
            j, scdict[j]["Het"], scdict[j]["Anis"], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames
        )
#        for t in [0,1,5,-1]:
#            axes.flat[Trial.index(j)].plot(np.mean(df[2, t, :, :], axis = -1), label = t)
        xindex = list(x * step for x in list(range(np.shape(conctime)[0])))
#        for t in [0,1,15,-1]:
#            axes.flat[Trial.index(j)].plot(conctime[t, :, 0], label = t)
        axes.flat[Trial.index(j)].plot(xindex, conctime[:,yout,0])
        axes.flat[Trial.index(j)].set_xlabel("")
        #        axes[colidx1][colidx2].set_xticklabels([])
        axes.flat[Trial.index(j)].set_ylabel("")
        axes.flat[Trial.index(j)].set_title(j)
        col = col + 1
    for ax in axes[:, 0]:
        ax.set_ylabel("Tracer (uM)")
    for ax in axes[-1]:
        ax.set_xlabel("Time (days)")
    fig.subplots_adjust(left=0.15, top=0.9)
    plt.legend()
    fig.savefig(
        "X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughcurve_"
        + Reg
        + ".png",
        dpi=300,
    )
    fig.savefig(
        "X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughcurve_"
        + Reg
        + ".pdf",
        dpi=300,
    )

f = "X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv"
csvfile = open(f, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    ["Sno", "Trial", "Variance", "Anisotropy", "Chem", "Time", "fraction", "Regime"]
)
idx = 1
for Reg, step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/" + Reg + "AR/"
    fpre = "NS-A"
    df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
        Trial[0], Het[0], Anis[0], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
    )
    print(np.mean(df[vely - 3, 1:, :, :]))
    Time = np.where(np.round(conctime[:, yout, 0], 3) > 10)
    initial = step * Time[0][0]
    for j in range(len(Trial)):
        df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
            Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
        )
        print(np.mean(df[vely - 3, 1:, :, :]))
        Time = np.where(np.round(conctime[:, yout, 0], 3) > 10)
        Time2 = np.where(np.round(df[-1, :, 50, :], 3) > 10)
        s = step
        if Reg == "Slow":
            if Trial[j] == 54:
                s = 100 * 0.01
        print(s * Time[0][0], s * Time2[0][0], initial, (s * Time[0][0]) / initial)
        writer.writerow(
            [
                idx,
                Trial[j],
                Het[j],
                Anis[j],
                "Tracer",
                s * Time[0][0],
                (s * Time[0][0]) / initial,
                Reg,
            ]
        )
        idx = idx + 1
csvfile.close()
# plotting boxplots to see variance of breakthrough from homogeneous scenario
tracerplot = sssp.plot_tracer()
tracerplot.savefig(
    "X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.png",
    dpi=300,
    pad_inches=0.01,
)
tracerplot = sssp.plot_tracer()
tracerplot.savefig(
    "X:/Saturated_flow/Steady_state/Tracer_studies/breakthroughfraction.png",
    dpi=300,
    pad_inches=0.01,
)