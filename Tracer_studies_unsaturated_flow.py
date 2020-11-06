# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import analyses.unsaturated_transient as uta
import plots.general_plots as gp

# Unsaturated flow regime
Reg = "Fast"
directory = r"X:/Richards_flow/Tracer_studies/" + Reg + "AR/"
fpre = "RF-A"
fsuf = r"/"
gw = 0

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

# setup what we really want to investigate
# Default:
# Variations:
#notlist = ['43','46', '78', '79']
#Trial = list(t for t,values in scdict.items() if t not in notlist)
#Het = list(values['Het'] for t,values in scdict.items() if t not in notlist)
#Anis = list(values['Anis'] for t,values in scdict.items() if t not in notlist)

# Constants
yout = 50
yin = 0
xleft = 0
xright = 30
tr1 = 8 - gw
vely = 5
velx = 4
vars = [tr1]
gvarnames = ["Tracer"]
vedge = 0.005
velem = 0.01
vbc = 0.3
por = 0.2
Regimes = ["Equal", "Fast", "Slow"]
steps = [20*0.005, 5 * 0.0002, 200*0.005]
#steps = [500 * 0.005, 200*0.005, 500 * 0.0002]
#Regimes = ["Slow", "Equal", "Fast"]

breakthrough2 = []
for Reg, step in zip(Regimes, steps):
    s = step
    if Reg == "Equal":
        r = "Medium"
    else:
        r = Reg
    d = r"X:/Richards_flow/Tracer_studies/" + Reg + "AR_eve/"
    fpre = "RF-A"
    df, conctime, masstime, Velocity, head = uta.calcconcmasstime('H', scdict['H']['Het'], scdict['H']['Het'], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
    print(np.mean(df[vely - 3, 1:, :, :]))
    Time = np.where(np.round(conctime[:, yout, 0], 3) > 10)
    initial = step * Time[0][0]
    for t in Trial:
        h = scdict[t]['Het']
        a = scdict[t]['Anis']
        df, conctime, masstime, Velocity, head = uta.calcconcmasstime(t, h, a, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
        print(np.mean(df[vely - 3, 1:, :, :]))
        Time = np.where(np.round(conctime[:, yout, 0], 3) > 10)
        Time2 = np.where(np.round(df[-1, :, 50, :], 3) > 10)
        s = step
        print(s * Time[0][0], s * Time2[0][0], initial, (s * Time[0][0]) / initial)
        breakthrough2.append([t, h, a, "Tracer", s*Time[0][0], (s*Time[0][0])/initial, r])

data = pd.DataFrame(breakthrough2, columns = ["Trial", "Variance", "Anisotropy", "Chem", "Time", "fraction", "Regime"])
f = r"X:/Richards_flow/Tracer_studies/tracer_equalfast_06112020.csv"
data.to_csv(f, sep = '\t')
        
# plotting boxplots to see variance of breakthrough from homogeneous scenario
tracerplot = gp.plot_tracer(f)
tracerplot.savefig("X:/Richards_flow/Tracer_studies/tracer_breakthrough_impact.png", dpi = 300, pad_inches = 0.1, bbox_inches = 'tight')

#data_exploration of tracer studies
import matplotlib.pyplot as plt
import seaborn as sns

for Reg in ["Slow"]:
    d = r"X:/Richards_flow/Tracer_studies/" + Reg + "AR/"
    fpre = "RF-A"
    count = 0
    for j in range(len(Trial)):
        df = np.load(d + fpre + str(Trial[j]) + fsuf + fpre + str(Trial[j]) + "_df.npy")
        #        v1 = np.sqrt(np.square(df[2,-2,:,:]) + np.square(df[2,-2,:,:]))
        #        v2 = np.sqrt(np.square(df[2,-1,:,:]) + np.square(df[2,-1,:,:]))
        v1 = df[2, -2, :, :]
        v2 = df[2, -1, :, :]
        diff = (np.abs(v2 - v1)) * 100 / 0.0038
        #        plt.figure()
        #        sns.heatmap(diff)
        #        plt.title(Trial[j])
        #        plt.savefig(d+str(Trial[j])+"_diff_velocities_finergrid_100days.png",dpi = 300, pad_inches = 0)
        arr = v1 == v2
        if arr.any():
            print(Trial[j], " steady")
            count = count + 1
        else:
            print(Trial[j], np.max(np.abs(v2 - v1)), np.argmax(np.abs(v1 - v2)))
    print(count)
    for j in range(len(Trial)):
        for t in [-10, -1]:
            plt.figure()
            sns.heatmap(df[2, t, :, :])
    plt.figure()
    for j in range(len(Trial)):
        df = np.load(d + fpre + str(Trial[j]) + fsuf + fpre + str(Trial[j]) + "_df.npy")
        ussp.plotdataindex(df, 0, t, "Pressure", Reg + str(Trial[j]))
    #        plt.legend()
    plt.figure()
    for j in range(len(Trial)):
        df = np.load(d + fpre + str(Trial[j]) + fsuf + fpre + str(Trial[j]) + "_df.npy")
        ussp.plotdataindex(df, -2, t, "Saturation", str(Trial[j]))
        plt.legend()
    plt.figure()
    for j in range(len(Trial)):
        df = np.load(d + fpre + str(Trial[j]) + fsuf + fpre + str(Trial[j]) + "_df.npy")
        veliredg = df[2, 1:, yin, xright]
        veliledg = df[2, 1:, yin, xleft]
        veloredg = df[2, 1:, yout, xright]
        veloledg = df[2, 1:, yout, xleft]
        veloelem = df[2, 1:, yout, xleft + 1 : xright]
        velielem = df[2, 1:, yin, xleft + 1 : xright]
        velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
        vellelem = df[2, 1:, yin + 1 : yout, xleft]
        velrelem = df[2, 1:, yin + 1 : yout, xright]
        satielem = df[4, 1:, yin, xleft + 1 : xright]
        satoelem = df[4, 1:, yout, xleft + 1 : xright]
        satlelem = df[4, 1:, yin + 1 : yout, xleft]
        satrelem = df[4, 1:, yin + 1 : yout, xright]
        satiredg = df[4, 1:, yin, xright]
        satiledg = df[4, 1:, yin, xleft]
        satoledg = df[4, 1:, yout, xleft]
        satoredg = df[4, 1:, yout, xright]
        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
        v = np.zeros([np.shape(df)[1] - 1, np.shape(df)[2]])
        v[:, yin] = (veliredg * satiredg + veliledg * satiledg) * vedge + np.sum(
            (velielem * satielem) * velem, axis=-1
        )  # /((satiledg+satiredg)*vedge + np.sum(satielem, axis = -1)*velem)#0.3
        v[:, yout] = (veloredg * satoredg + veloledg * satoledg) * vedge + np.sum(
            (veloelem * satoelem) * velem, axis=-1
        )  # /((satoledg+satoredg)*vedge + np.sum(satoelem, axis = -1)*velem)#0.3
        v[:, yin + 1 : yout] = (
            velrelem * satrelem + vellelem * satlelem
        ) * vedge + np.sum(
            (velelem * satelem) * velem, axis=-1
        )  # /((satlelem+satrelem)*vedge + np.sum(satelem, axis =-1)*velem)#0.3
        v = np.zeros([np.shape(df)[1] - 1, np.shape(df)[2]])
        v[:, yin] = (veliredg + veliledg) * vedge + np.sum(velielem * velem, axis=-1)
        v[:, yout] = (veloredg + veloledg) * vedge + np.sum(veloelem * velem, axis=-1)
        v[:, yin + 1 : yout] = (velrelem + vellelem) * vedge + np.sum(
            velelem * velem, axis=-1
        )
        plt.plot(v[t, :], label=Reg + str(Trial[j]))
        #        plt.ylim((-0.00113, -0.00115))#-0.004,0)
        plt.ylabel("Volumetric flow rate (m3/d)")
        plt.xlabel("Y (cm)")
        plt.legend()
    plt.figure()
    for j in range(len(Trial)):
        df, conctime, mftime, Velocity, Headinlettime = uta.calcconcmasstime(
            Trial[j],
            Het[j],
            Anis[j],
            gw,
            directory,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
            gvarnames,
        )
        plt.plot(conctime[-1, :], label=Reg + str(Trial[j]))
        plt.ylabel("Concentration of tracer (uM)")
        plt.xlabel("Y (cm)")
        #        plt.ylim((0.0148,0.015))
        plt.title("Tracer concentration in the domain at steady state conditions")
        plt.legend()
    plt.figure()
    for j in range(len(Trial)):
        df, conctime, mftime, Velocity, Headinlettime = uta.calcconcmasstime(
            Trial[j],
            Het[j],
            Anis[j],
            gw,
            directory,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
            gvarnames,
        )
        plt.plot(conctime[1:, yout], label=Reg + str(Trial[j]))
        plt.ylabel("Tracer (uM)")
        plt.xlabel("Time step")
        plt.title("Tracer concentration at outlet at steady state conditions")
        plt.ylim((0, 20.1))
        #        plt.legend()

        plt.plot(np.mean(df[0, -1, :, :], axis=-1), label=Reg + str(Trial[j]))
        plt.ylabel("Pressure")
        plt.xlabel("Y (cm)")
        plt.plot(np.mean(df[-1, -1, :, :], axis=-1), label=Reg + str(Trial[j]))
        plt.ylabel("Tracer (uM)")
        plt.xlabel("Y (cm)")

        plt.plot(np.mean(df[-2, -1, :, :], axis=-1), label=Reg + str(Trial[j]))
        plt.ylabel("Saturation")
        plt.xlabel("Y (cm)")

        plt.plot(v[1, :], label=Reg + str(Trial[j]))
        plt.plot(v[2, :], label=Reg + str(Trial[j]))
        plt.plot(v[15, :], label=Reg + str(Trial[j]))
        plt.legend()

        plt.plot(conctime[-1, :], label=Reg + str(Trial[j]))
        plt.plot(conctime[1, :], label=Reg + str(Trial[j]))
        plt.plot(conctime[2, :], label=Reg + str(Trial[j]))
        plt.plot(conctime[15, :], label=Reg + str(Trial[j]))
        plt.legend()

        plt.plot(v[-1, :], label=Reg + str(Trial[j]))
        plt.legend()

        if ((Reg == "Equal") & (avgvel > -0.0035)) or (
            (Reg == "Fast") & (avgvel < -0.038)
        ):
            print(avgvel)
            print(
                np.max(
                    (
                        (v[1:, yin] + v[1:, yout]) * 0.005
                        + np.sum(v[1:, yin + 1 : yout], axis=-1) * 0.01
                    )
                    / 0.5
                )
            )

        #        print(v[-1,yin])
        #        print(v[-1,yout-1])
        #        print(v[-1,yout])
        #        print(np.mean(df[-2,1,:,:]))
        print(avgvel)
        plt.plot(np.mean(v[1:, :], axis=-1), label=Reg + str(Trial[j]))
        plt.legend()

        plt.plot(conctime[:, yout], label=Reg + str(Trial[j]))
        plt.legend()
        print(np.mean(df[-2, -1, :, :]))  # , axis = -1))
#        plt.plot(np.mean(df[0,-1,:, :], axis = -1),label = Reg + str(Trial[j]))
#        plt.legend()
#       print(np.mean(df[-2,0,yout, :], axis = -1))
#        print(conctime[-1,yout])
