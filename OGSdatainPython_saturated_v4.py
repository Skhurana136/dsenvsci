# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""
import numpy as np
import csv
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import analyses.saturated_steady_state as sssa
import plots.saturated_steady_state as sssp
import data_reader.reader as rdr
import data_reader.data_processing as proc
import analyses.saturated_transient as sta

# Saturated flow regime
Reg = "Fast"
directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
fpre = "NS-A"
fsuf = r"/"
gw = 1

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())


# setup what we really want to investigate
# Default:
# Variations:
#notlist = [43,54]
#Trial = list(t for t,values in scdict.items() if t not in notlist)
#Het = list(values['Het'] for t,values in scdict.items() if t not in notlist)
#Anis = list(values['Anis'] for t,values in scdict.items() if t not in notlist)

# Constants
yout = 50
yin = 0
xleft = 0
xright = 30
# Assign index to Variable
doc1 = 10 - gw
dox1 = 11 - gw
Amm1 = 12 - gw
nitra1 = 17 - gw
sulpha1 = 22 - gw
tr1 = 29 - gw
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
vely = 5
velx = 4
vars = [doc1, dox1, Amm1, nitra1, sulpha1, Bmo1, Bma1, Bmn1, Bimo1, Bima1, Bimn1]
gvarnames = [
    "DOC",
    "DO",
    "Ammonium",
    "Nitrate",
    "Sulphate",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
    "Nitrogen",
    "TOC",
]
AFbiomassvars = [
    Bfo1,
    Bfa1,
    Bfn1,
    Bmo1,
    Bma1,
    Bmn1,
    Bifo1,
    Bifa1,
    Bifn1,
    Bimo1,
    Bima1,
    Bimn1,
]
AFbiomassgvarnames = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive fixed Aerobes",
    "Inactive fixed Ammonia oxidizers",
    "Inactive fixed Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
]

mf = sssa.calcmassflux(
    Trial,
    Het,
    Anis,
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
biomasssum = sssa.calcsum(
    Trial,
    Het,
    Anis,
    gw,
    directory,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    AFbiomassvars,
    AFbiomassgvarnames,
)
meanresultsarr, stdresultsarr, covresultsarr = sssa.calcaggres(mf)

# Writing the results into csv files
fnamemf = (
    directory
    + fpre
    + str(Trial[0])
    + "_"
    + str(Trial[-1])
    + "_mass_flux"
    + Reg
    + str(gw)
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".csv"
)
csvfile = open(fnamemf, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    [
        "Sno",
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Inlet_total_mass_flux",
        "Outlet_mass_flux",
        "delmassflux",
        "del2massflux",
        "fdelmassflux",
    ]
)
for j in Trial:
    for i in gvarnames:
        idx = Trial.index(j) * len(gvarnames) + gvarnames.index(i)
        writer.writerow(
            [
                idx + 1,
                j,
                Het[Trial.index(j)],
                Anis[Trial.index(j)],
                i,
                mf[idx, 4],
                mf[idx, 5],
                mf[idx, 6],
                mf[idx, 7],
                mf[idx, 8],
            ]
        )
csvfile.close()

# Mass balance
yin = 0
yout = 50
vedge = 0.005
velem = 0.01
vbc = 0.3
por = 0.2
diff = 10 ** -8
disp = 0.1
Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
mf = np.zeros([len(Trial) * (len(gvarnames)), 9])
for j in range(len(Trial)):
    di = directory + fpre + str(Trial[j]) + fsuf
    print(str(Trial[j]))
    df = np.load(di + fpre + str(Trial[j]) + "_df.npy")
    veliredg = df[2, -1, yin, xright]
    veliledg = df[2, -1, yin, xleft]
    veloredg = df[2, -1, yout, xright]
    veloledg = df[2, -1, yout, xleft]
    veloelem = df[2, -1, yout, xleft + 1 : xright]
    velielem = df[2, -1, yin, xleft + 1 : xright]
    velelem = df[2, -1, yin + 1 : yout, xleft + 1 : xright]
    vellelem = df[2, -1, yin + 1 : yout, xleft]
    velrelem = df[2, -1, yin + 1 : yout, xright]
    if gw == 1:
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
        satielem = df[4, -1, yin, xleft + 1 : xright]
        satoelem = df[4, -1, yout, xleft + 1 : xright]
        satiredg = df[4, -1, yin, xright]
        satiledg = df[4, -1, yin, xright]
        satoledg = df[4, -1, yout, xleft]
        satoredg = df[4, -1, yout, xright]
        satlelem = df[4, -1, yin + 1 : yout, xleft]
        satrelem = df[4, -1, yin + 1 : yout, xright]
        satelem = df[4, -1, yin + 1 : yout, xleft + 1 : xright]
    for i in range(len(gvarnames)):
        idx = j * (len(gvarnames)) + i
        mf[idx, 0] = j
        mf[idx, 1] = Het[j]
        mf[idx, 2] = Anis[j]
        mf[idx, 3] = i
        if gvarnames[i] == "Nitrogen":
            ninlet = 0
            noutlet = 0
            for n in Nspecies:
                ninlet = (
                    ninlet
                    + df[n - 3, -1, yin, xleft]
                    * satiledg
                    * (diff + disp * abs(veliledg))
                    + df[n - 3, -1, yin, xright]
                    * satiredg
                    * (diff + disp * abs(veliredg))
                    + sum(
                        df[n - 3, -1, yin, xleft + 1 : xright]
                        * satielem
                        * (diff + disp * abs(velielem))
                    )
                    + abs(
                        df[n - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                        + df[n - 3, -1, yin, xright] * satiredg * veliredg * vedge
                        + sum(
                            df[n - 3, -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem
                        )
                    )
                )
                noutlet = (
                    noutlet
                    + df[n - 3, -1, yout, xleft]
                    * satoledg
                    * (diff + disp * abs(veloledg))
                    + df[n - 3, -1, yout, xright]
                    * satoredg
                    * (diff + disp * abs(veloredg))
                    + sum(
                        df[n - 3, -1, yout, xleft + 1 : xright]
                        * satoelem
                        * (diff + disp * abs(veloelem))
                    )
                    + abs(
                        df[n - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                        + df[n - 3, -1, yout, xright] * satoredg * veloredg * vedge
                        + sum(
                            df[n - 3, -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem
                        )
                    )
                )
            ninlet = (
                ninlet / 5
                + df[Amm1 - 3, -1, yin, xleft]
                * satiledg
                * (diff + disp * abs(veliledg))
                + df[Amm1 - 3, -1, yin, xright]
                * satiredg
                * (diff + disp * abs(veliredg))
                + sum(
                    df[Amm1 - 3, -1, yin, xleft + 1 : xright]
                    * satielem
                    * (diff + disp * abs(velielem))
                )
                + df[nitra1 - 3, -1, yin, xleft]
                * satiledg
                * (diff + disp * abs(veliledg))
                + df[nitra1 - 3, -1, yin, xright]
                * satiredg
                * (diff + disp * abs(veliredg))
                + sum(
                    df[nitra1 - 3, -1, yin, xleft + 1 : xright]
                    * satielem
                    * (diff + disp * abs(velielem))
                )
                + (
                    df[Amm1 - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                    + df[Amm1 - 3, -1, yin, xright] * satiredg * veliredg * vedge
                    + sum(
                        df[Amm1 - 3, -1, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem
                    )
                    + df[nitra1 - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                    + df[nitra1 - 3, -1, yin, xright] * satiredg * veliredg * vedge
                    + sum(
                        df[nitra1 - 3, -1, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem
                    )
                )
            )
            noutlet = (
                noutlet / 5
                + df[Amm1 - 3, -1, yout, xleft]
                * satoledg
                * (diff + disp * abs(veloledg))
                + df[Amm1 - 3, -1, yout, xright]
                * satoredg
                * (diff + disp * abs(veloredg))
                + sum(
                    df[Amm1 - 3, -1, yout, xleft + 1 : xright]
                    * satoelem
                    * (diff + disp * abs(veloelem))
                )
                + df[nitra1 - 3, -1, yout, xleft]
                * satoledg
                * (diff + disp * abs(veloledg))
                + df[nitra1 - 3, -1, yout, xright]
                * satoredg
                * (diff + disp * abs(veloredg))
                + sum(
                    df[nitra1 - 3, -1, yout, xleft + 1 : xright]
                    * satoelem
                    * (diff + disp * abs(veloelem))
                )
                + (
                    df[Amm1 - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                    + df[Amm1 - 3, -1, yout, xright] * satoredg * veloredg * vedge
                    + sum(
                        df[Amm1 - 3, -1, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem
                    )
                    + df[nitra1 - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                    + df[nitra1 - 3, -1, yout, xright] * satoredg * veloredg * vedge
                    + sum(
                        df[nitra1 - 3, -1, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem
                    )
                )
            )
            mf[idx, 4] = abs(
                ninlet
            )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
            mf[idx, 5] = abs(
                noutlet
            )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
        elif gvarnames[i] == "TOC":
            cinlet = 0
            coutlet = 0
            for c in Cspecies:
                cinlet = (
                    cinlet
                    + df[c - 3, -1, yin, xleft]
                    * satiledg
                    * (diff + disp * abs(veliledg))
                    + df[c - 3, -1, yin, xright]
                    * satiredg
                    * (diff + disp * abs(veliredg))
                    + sum(
                        df[c - 3, -1, yin, xleft + 1 : xright]
                        * satielem
                        * (diff + disp * abs(velielem))
                    )
                    + abs(
                        df[c - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                        + df[c - 3, -1, yin, xright] * satiredg * veliredg * vedge
                        + sum(
                            df[c - 3, -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem
                        )
                    )
                )
                coutlet = (
                    coutlet
                    + df[c - 3, -1, yin, xleft]
                    * satiledg
                    * (diff + disp * abs(veliledg))
                    + df[c - 3, -1, yin, xright]
                    * satiredg
                    * (diff + disp * abs(veliredg))
                    + sum(
                        df[c - 3, -1, yin, xleft + 1 : xright]
                        * satielem
                        * (diff + disp * abs(velielem))
                    )
                    + abs(
                        df[c - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                        + df[c - 3, -1, yout, xright] * satoredg * veloredg * vedge
                        + sum(
                            df[c - 3, -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem
                        )
                    )
                )
            mf[idx, 4] = abs(
                cinlet
            )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
            mf[idx, 5] = abs(
                coutlet
            )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
        else:
            #                mf[idx,4] = (df[vars[i]-3,-1,yin,xleft]-df[vars[i]-3,-1,yin+1,xleft])*satiledg*(disp*abs(veliledg))+(df[vars[i]-3,-1,yin,xright]-df[vars[i]-3,-1,yin+1,xright])*satiredg*(diff+disp*abs(veliredg))+sum((df[vars[i]-3,-1,yin,xleft+1:xright]-df[vars[i]-3,-1,yin+1,xleft+1:xright])*satielem*(diff+disp*abs(velielem)))+abs(((df[vars[i]-3,-1,yin,xleft]*satiledg*veliledg*vedge + df[vars[i]-3,-1,yin,xright]*satiredg*veliredg*vedge + sum(df[vars[i]-3,-1,yin,xleft+1:xright]*satielem*velielem*velem))))/por#/(sum(velielem)*velem + (veliledg+veliredg)*vedge)
            #                mf[idx,5] = (df[vars[i]-3,-1,yout-1,xleft]-df[vars[i]-3,-1,yout,xleft])*satoledg*(diff+disp*abs(veloledg))+(df[vars[i]-3,-1,yout-1,xright]-df[vars[i]-3,-1,yout,xright])*satoredg*(diff+disp*abs(veloredg))+sum((df[vars[i]-3,-1,yout-1,xleft+1:xright]-df[vars[i]-3,-1,yout,xleft+1:xright])*satoelem*(diff+disp*abs(veloelem)))+abs(((df[vars[i]-3,-1,yout,xleft]*satoledg*veloledg*vedge + df[vars[i]-3,-1,yout,xright]*satoredg*veloredg*vedge + sum(df[vars[i]-3,-1,yout,xleft+1:xright]*satoelem*veloelem*velem))))/por#/(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            mf[idx, 4] = (
                (df[vars[i] - 3, -1, yin, xleft] - df[vars[i] - 3, -1, yin + 1, xleft])
                * satiledg
                * disp
                + (
                    df[vars[i] - 3, -1, yin, xright]
                    - df[vars[i] - 3, -1, yin + 1, xright]
                )
                * satiredg
                * disp
                + sum(
                    (
                        df[vars[i] - 3, -1, yin, xleft + 1 : xright]
                        - df[vars[i] - 3, -1, yin + 1, xleft + 1 : xright]
                    )
                    * satielem
                    * disp
                )
                + abs(
                    (
                        (
                            df[vars[i] - 3, -1, yin, xleft]
                            * satiledg
                            * veliledg
                            * vedge
                            + df[vars[i] - 3, -1, yin, xright]
                            * satiredg
                            * veliredg
                            * vedge
                            + sum(
                                df[vars[i] - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                    )
                )
                / por
            )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
            mf[idx, 5] = (
                (
                    df[vars[i] - 3, -1, yout - 1, xleft]
                    - df[vars[i] - 3, -1, yout, xleft]
                )
                * satoledg
                * disp
                + (
                    df[vars[i] - 3, -1, yout - 1, xright]
                    - df[vars[i] - 3, -1, yout, xright]
                )
                * satoredg
                * disp
                + sum(
                    (
                        df[vars[i] - 3, -1, yout - 1, xleft + 1 : xright]
                        - df[vars[i] - 3, -1, yout, xleft + 1 : xright]
                    )
                    * satoelem
                    * disp
                )
                + abs(
                    (
                        (
                            df[vars[i] - 3, -1, yout, xleft]
                            * satoledg
                            * veloledg
                            * vedge
                            + df[vars[i] - 3, -1, yout, xright]
                            * satoredg
                            * veloredg
                            * vedge
                            + sum(
                                df[vars[i] - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                    )
                )
                / por
            )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
        mf[idx, 6] = mf[idx, 4] - mf[idx, 5]

f = "X:/massbalance_file1_" + Reg + ".csv"
csvfile = open(f, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    [
        "Sno",
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Inlet_massflux_umolesperday",
        "Outlet_massflux_umolesperday",
        "Moles_removed_per_day",
    ]
)
for j in range(len(Trial)):
    for i in range(len(gvarnames)):
        idx = j * len(gvarnames) + i
        writer.writerow(
            [
                idx,
                Trial[j],
                Het[j],
                Anis[j],
                gvarnames[i],
                mf[idx, 4],
                mf[idx, 5],
                mf[idx, 6],
            ]
        )
csvfile.close()

mastermf = np.zeros([1, 9])
masterbiomasssum = np.zeros([1, 7])
Regimes = ["Slow", "Equal", "Fast"]
for Reg in Regimes:
    d = "Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    fpre = "NS-A"
    mf = sssa.calcmassflux(
        Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames
    )
    biomasssum = sssa.calcsum(
        Trial,
        Het,
        Anis,
        gw,
        d,
        fpre,
        fsuf,
        yin,
        yout,
        xleft,
        xright,
        AFbiomassvars,
        AFbiomassgvarnames,
    )
    mastermf = np.append(mastermf, mf, axis=0)
    masterbiomasssum = np.append(masterbiomasssum, biomasssum, axis=0)
mastermf = np.delete(mastermf, 0, 0)
masterbiomasssum = np.delete(masterbiomasssum, 0, 0)

del4massflux = np.zeros([len(Trial) * len(Regimes) * len(gvarnames), 1])
del2biomass = np.zeros([len(Trial) * len(Regimes) * len(AFbiomassgvarnames), 1])
regindex = np.zeros([len(Trial) * len(Regimes) * len(gvarnames), 1])
regindexb = np.zeros([len(Trial) * len(Regimes) * len(AFbiomassgvarnames), 1])
for Reg in Regimes:
    for k in Trial:
        for i in gvarnames:
            idx = (
                Regimes.index(Reg) * len(Trial) * len(gvarnames)
                + Trial.index(k) * len(gvarnames)
                + gvarnames.index(i)
            )
            del4massflux[idx, 0] = (mastermf[idx, 6]) / mastermf[gvarnames.index(i), 6]
            if k == "H":
                mastermf[idx, 0] = 0
            else:
                mastermf[idx, 0] = k
        for b in AFbiomassgvarnames:
            idxb = (
                Regimes.index(Reg) * len(Trial) * len(AFbiomassgvarnames)
                + Trial.index(k) * len(AFbiomassgvarnames)
                + AFbiomassgvarnames.index(b)
            )
            del2biomass[idxb, 0] = (masterbiomasssum[idxb, 4]) / masterbiomasssum[
                AFbiomassgvarnames.index(b), 4
            ]
            if k == "H":
                masterbiomasssum[idxb, 0] = 0
            else:
                masterbiomasssum[idxb, 0] = k
print(idx, idxb)

for Reg in Regimes:
    for k in Trial:
        for i in gvarnames:
            idx = (
                Regimes.index(Reg) * len(Trial) * len(gvarnames)
                + Trial.index(k) * len(gvarnames)
                + gvarnames.index(i)
            )
            regindex[idx, 0] = Regimes.index(Reg)
        for b in AFbiomassgvarnames:
            idxb = (
                Regimes.index(Reg) * len(Trial) * len(AFbiomassgvarnames)
                + Trial.index(k) * len(AFbiomassgvarnames)
                + AFbiomassgvarnames.index(b)
            )
            regindexb[idxb, 0] = Regimes.index(Reg)

mastermf = np.append(mastermf, del4massflux, 1)
masterbiomasssum = np.append(masterbiomasssum, del2biomass, 1)
mastermf = np.append(mastermf, regindex, 1)
masterbiomasssum = np.append(masterbiomasssum, regindexb, 1)

dfallmf = pd.DataFrame(
    data=mastermf,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Inlet_total_mass_flux",
        "Outlet_mass_flux",
        "delmassflux",
        "reldelmassflux",
        "del2massflux",
        "del4massflux",
        "Regime",
    ],
)
dfall2 = proc.processdataframe(dfallmf, gvarnames)
dfallbiomass = pd.DataFrame(
    data=masterbiomasssum,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Total_biomass",
        "Change_umoles",
        "fractionoftotal",
        "del2biomass",
        "Regime",
    ],
)
dfallbiomass2 = proc.processdataframe(dfallbiomass, AFbiomassgvarnames)

filename = "X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv"
breakthrough = pd.read_csv(filename, delimiter="\t")
breakthrough["Heterogeneity"] = breakthrough["Variance"]
breakthrough["VA"] = breakthrough["Heterogeneity"] * breakthrough["Anisotropy"]
breakthrough["%ofhomogeneous"] = breakthrough["fraction"] * 100
breakthrough.loc[breakthrough.Regime == "Equal", "Regime"] = "Medium"

dfall2 = pd.merge(
    dfallmf,
    breakthrough[["Trial", "Regime", "Time", "fraction", "%ofhomogeneous"]],
    on=["Trial", "Regime"],
).rename(columns={"Time": "Breakthroughtime"})
dfall2["%del4massflux"] = dfall2["del4massflux"] * 100
dfall2.to_csv(
    "Z:/Saturated_flow/diffusion_transient/massflux_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)

filename = "Z:/Saturated_flow/diffusion_transient/massflux_withbreakthrough_forMartin_v4_complete.csv"
dfall2 = pd.read_csv(filename, delimiter="\t")
dummy = sssp.norm_mf(dfall2, gvarnames)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_impact_massflux.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_impact_massflux.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)
dummy = sssp.norm_mf_2x2(dfall2, gvarnames)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_impact_massflux_square.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_impact_massflux_square.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)


dfall2biomass = pd.merge(
    dfallbiomass2,
    breakthrough[["Trial", "Regime", "Time", "fraction", "%ofhomogeneous"]],
    on=["Trial", "Regime"],
).rename(columns={"Time": "Breakthroughtime"})
dfall2biomass["%del2biomass"] = dfall2biomass["del2biomass"] * 100
dfall2biomass.to_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)

activebiomassgvarnames = AFbiomassgvarnames[:6]
inactivebiomassgvarnames = AFbiomassgvarnames[6:]
dummy = sssp.steadystate_active_biomass(dfall2biomass, activebiomassgvarnames)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_active_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_active_biomass.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)

dummy = sssp.steadystate_active_biomass(dfall2biomass, inactivebiomassgvarnames)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_inactive_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)
dummy.savefig(
    "Z:/Saturated_flow/diffusion_transient/steadystate_inactive_biomass.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)

fnamebiosum = (
    directory
    + fpre
    + "Anaerobic"
    + str(Trial[0])
    + "_"
    + str(Trial[-1])
    + "_"
    + "biomass_"
    + str(gw)
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".csv"
)
csvfile = open(fnamebiosum, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    ["Sno", "Trial", "Variance", "Anisotropy", "Chem", "Total_biomass", "Change_umoles"]
)
for j in Trial:
    for i in AFbiomassgvarnames:
        idx = Trial.index(j) * len(AFbiomassgvarnames) + AFbiomassgvarnames.index(i)
        writer.writerow(
            [
                idx + 1,
                j,
                Het[Trial.index(j)],
                Anis[Trial.index(j)],
                i,
                biomasssum[idx, 4],
                biomasssum[idx, 5],
            ]
        )
csvfile.close()

fnamemf_summary = (
    directory
    + fpre
    + "Anaerobic"
    + str(Trial[0])
    + "_"
    + str(Trial[-1])
    + "_"
    + "new_MF_"
    + str(gw)
    + datetime.now().strftime("%d%m%Y%H%M")
    + "_summary.csv"
)
csvfile = open(fnamemf_summary, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    [
        "Sno",
        "Variance",
        "Anisotropy",
        "Chem",
        "Removal",
        "Change_removalrate",
        "Std_var_RR",
    ]
)
for j in range(len(meanresultsarr)):
    writer.writerow(
        [
            j + 1,
            meanresultsarr[j, 2],
            meanresultsarr[j, 3],
            meanresultsarr[j, 4],
            meanresultsarr[j, 0] * 100,
            meanresultsarr[j, 1] * 100,
            stdresultsarr[j, 0],
        ]
    )
csvfile.close()

# 1D analysis:
# Importing data
if Reg == "Fast":
    #    Reg = "Fast"
    #    fnamemf = "X:/Saturated_flow/diffusion/"+Reg+"AR_changedkindox/NS-AH_84_mass_fluxFast1030320201111.csv"
    #    fnamebiosum = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_84_biomass_1030120201028.csv"
    fastb = proc.processchembiomassfiles(fnamebiosum, Reg)
    fastc = proc.processchemfiles(fnamemf, Reg)
elif Reg == "Slow":
    Reg = "Slow"
    #    fnamemf = "X:/Saturated_flow/diffusion/"+Reg+"AR_changedkindox/NS-AH_84_mass_fluxSlow1030320201110.csv"
    #    fnamebiosum = "X:/Saturated_flow/Anaerobic/"+Reg+"AR_changedkindox/NS-AAnaerobicH_83_biomass_1070120201707.csv"
    slowb = proc.processchembiomassfiles(fnamebiosum, Reg)
    slowc = proc.processchemfiles(fnamemf, Reg)
elif Reg == "Equal":
    Reg = "Equal"
    #    fnamemf = "X:/Saturated_flow/diffusion/"+Reg+"AR_changedkindox/NS-AH_84_mass_fluxEqual1030320201112.csv"
    equalb = proc.processchembiomassfiles(fnamebiosum, "Medium")
    equalc = proc.processchemfiles(fnamemf, "Medium")

# Scatter plot Flux with residence time
sgvarnames = [elem for elem in gvarnames if elem != "Sulphate"]
sngvarnames = [elem for elem in sgvarnames if elem != "Nitrogen"]
comb, dummy = sssp.scatterrestime_flux(fastc, slowc, equalc, sngvarnames)
picname = (
    "X:/Saturated_flow/Anaerobic/flux_restime_SCATTER_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname, bbox_inches="tight")
# Scatter plot biomass wtih residence time
comb, dummy = sssp.scatterrestime_biomass(fastb, slowb, equalb)
picname = (
    "X:/Saturated_flow/Anaerobic/mobile_biomass_restime_SCATTER_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)

# Concentration profiles
nrows = 8
ncols = 6
figsize = [21, 28]
colors = ["black", "red", "blue", "green", "orange", "grey"]
# Plot all the scenarios excluding homogeneous
intsce = Trial
sssp.plotconcallatss(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
    nrows,
    ncols,
    figsize,
)

# changing oxic-anoxic interface with heterogeneity
intsce = ["H", 50, 76, 73, 80, 84, 44, 63]
colorseries = ["Reds", "Greens", "Blues"]
chem = "DO"
Regimes = ["Slow", "Equal", "Fast"]
dummy = sssp.plotconcwithhet(
    Regimes,
    intsce,
    chem,
    colorseries,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
)
picname1 = (
    "X:/Saturated_flow/diffusion/"
    + str(intsce[0])
    + "_"
    + str(intsce[-1])
    + "_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + "_"
    + chem
    + "with_het.pdf"
)
picname2 = (
    "X:/Saturated_flow/diffusion/"
    + str(intsce[0])
    + "_"
    + str(intsce[-1])
    + "_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + "_"
    + chem
    + "with_het.png"
)
dummy.savefig(picname1, dpi=300)
dummy.savefig(picname2, dpi=300)

# Zoom into select heterogeneous scenarios

Reg = "Equal"
intsce = ["H", 44]
for ctype in ["chem"]:
    fig, ax = plt.subplots(ncols=len(intsce), nrows=1, figsize=[15, 8])
    for a, i in zip(ax, intsce):
        d = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
        f = "NS-A" + str(i) + "__" + ctype + "_conc_withH.png"
        im = plt.imread(d + f)
        a.imshow(im)
        #        a.set_title (Reg+" flow", fontsize = 15)
        a.axis("off")
        plt.savefig(
            "//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_"
            + ctype
            + ".pdf",
            dpi=300,
        )
        plt.savefig(
            "//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_"
            + ctype
            + ".png",
            dpi=300,
        )

Regimes = ["Slow", "Equal", "Fast"]
intsce = ["H"]
i = "H"
for ctype in ["chem"]:
    fig, ax = plt.subplots(ncols=len(Regimes), nrows=1, figsize=[15, 8])
    for a, Reg in zip(ax, Regimes):
        d = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
        f = "NS-A" + str(i) + "__" + ctype + "_conc_withH.png"
        im = plt.imread(d + f)
        a.imshow(im)
        if Reg == "Equal":
            a.set_title("Medium flow", fontsize=15)
        else:
            a.set_title(Reg + " flow", fontsize=15)
        a.axis("off")
        plt.savefig(
            "//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_"
            + ctype
            + ".pdf",
            dpi=300,
        )
        plt.savefig(
            "//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_"
            + ctype
            + ".png",
            dpi=300,
        )

# Zoom into select heterogeneous scenarios
intsce = ["H", 50, 76, 73, 84, 44, 63]
intsce = [63]
Regimes = ["Slow", "Equal", "Fast"]
chemvars = [doc1, dox1, Amm1, nitra1]
gvarnames = ["DOC", "DO", "Ammonium", "Nitrate"]
AMbiomassvars = [Bmo1, Bma1, Bmn1]
AMbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
AFbiomassvars = [Bfo1, Bfa1, Bfn1]
AFbiomassgvarnames = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]

sssp.concprofileatssX(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
)
sssp.concprofileatssX(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    AFbiomassvars,
    AFbiomassgvarnames,
    intsce,
    colors,
)

# Plot tracer with time
intsce = [
    37,
    38,
    39,
    46,
    47,
    48,
    55,
    56,
    57,
    64,
    65,
    66,
    40,
    41,
    42,
    49,
    50,
    51,
    58,
    59,
    60,
    67,
    68,
    69,
    70,
    71,
    72,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    43,
    44,
    45,
    52,
    53,
    54,
    61,
    62,
    63,
    73,
    74,
    75,
]
inttime = [2, 5, 8, 15, 18, 20]
chem = "DO"
lylim = 0
uylim = 5
nrows = 2
ncols = 2
sssp.plottimeseries(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    inttime,
    colors,
    nrows,
    ncols,
    figsize,
    chem,
    lylim,
    uylim,
)

# Plot velocity 2-D profile
intsce = ["H", 37, 38, 39, 40, 41, 51, 82, 80, 84, 53, 54]
indices = list(masterTrial.index(i) for i in intsce)
inthet = list(masterHet[i] for i in indices)
intanis = list(masterAnis[i] for i in indices)
for Reg in ["Slow", "Equal", "Fast"]:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    dummy = sssp.heatmapvelocity(directory, fpre, fsuf, intsce, inthet, intanis)
    picname = directory + "velocity_" + Reg + ".pdf"
    dummy.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)
    picname = directory + "velocity_" + Reg + ".png"
    dummy.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

# Classify oxic cells and then plot along Y axis
intsce = ["H", 50, 76, 73, 80, 84, 44, 63]
for Reg in Regimes:
    d = r"X:/Saturated_flow/diffusion/" + Reg + "AR_changedkindox/"
    fig = sssp.plotoxiccellssdo(
        20, intsce, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
    )
    picname = (
        directory + "DOandoxiccells_" + datetime.now().strftime("%d%m%Y%H%M") + ".png"
    )
    fig.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)
    picname = (
        directory + "DOandoxiccells_" + datetime.now().strftime("%d%m%Y%H%M") + ".pdf"
    )
    fig.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

Regimes = ["Slow", "Equal", "Fast"]
for Reg in Regimes:
    d = r"X:/Saturated_flow/diffusion/" + Reg + "AR_changedkindox/"
    for k in intsce:
        df = np.load(d + fpre + str(k) + fsuf + fpre + str(k) + "_df.npy")
        title = (
            "Variance "
            + str(Het[Trial.index(k)])
            + " : Anisotropy "
            + str(Anis[Trial.index(k)])
        )
        heatmapc = sssp.heatmapconcdist(df, vars, k, gvarnames, d, fpre, title)
        heatmapb = sssp.heatmapconcdistLOG(
            df, AFbiomassvars, k, AFbiomassgvarnames, d, fpre, title
        )
        picname = directory + "mobile_biomass_" + Reg + ".pdf"
        heatmapb.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

intsce = [50, 73, 63]
for ctype in ["chem", "biomass"]:
    fig, ax = plt.subplots(ncols=len(intsce), nrows=len(Regimes), figsize=[15, 15])
    for Reg in Regimes:
        d = r"X:/Saturated_flow/diffusion/" + Reg + "AR_changedkindox/"
        for a, i in zip(ax[Regimes.index(Reg)], intsce):
            f = "NS-A" + str(i) + "__" + ctype + "_conc_withH.png"
            im = plt.imread(d + f)
            a.imshow(im)
            #            a.set_title (Reg+" flow", fontsize = 15)
            a.axis("off")
        if Reg == "Equal":
            ax[Regimes.index(Reg), 0].annotate(
                "Medium flow",
                xy=(-0.04, 0.4),
                xytext=(0, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="left",
                va="center",
                rotation="vertical",
                fontsize=15,
            )
        else:
            ax[Regimes.index(Reg), 0].annotate(
                Reg + " flow",
                xy=(-0.04, 0.4),
                xytext=(0, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="left",
                va="center",
                rotation="vertical",
                fontsize=15,
            )
    plt.subplots_adjust(wspace=0.025, hspace=0.025)
    plt.savefig(
        "X:/Saturated_flow/diffusion/1D_heterogeneous_profiles_" + ctype + ".pdf",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.savefig(
        "X:/Saturated_flow/diffusion/1D_heterogeneous_profiles_" + ctype + ".png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
    )

for ctype in ["chem", "biomass"]:
    fig, ax = plt.subplots(ncols=len(intsce), nrows=len(Regimes), figsize=[15, 15])
    for Reg in Regimes:
        d = r"X:/Saturated_flow/diffusion/" + Reg + "AR_changedkindox/"
        for a, i in zip(ax[Regimes.index(Reg)], intsce):
            f = "NS-A" + str(i) + "_heatmap_" + ctype + ".png"
            im = plt.imread(d + f)
            a.imshow(im)
            #            a.set_title (Reg+" flow", fontsize = 15)
            a.axis("off")
        if Reg == "Equal":
            ax[Regimes.index(Reg), 0].annotate(
                "Medium flow",
                xy=(0, 0.5),
                xytext=(0, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="left",
                va="center",
                rotation="vertical",
                fontsize=15,
            )
        else:
            ax[Regimes.index(Reg), 0].annotate(
                Reg + " flow",
                xy=(0, 0.5),
                xytext=(0, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="left",
                va="center",
                rotation="vertical",
                fontsize=15,
            )
    plt.subplots_adjust(wspace=0.025, hspace=0.025)
    plt.savefig(
        "X:/Saturated_flow/diffusion/2D_heterogeneous_profiles_" + ctype + ".pdf",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.savefig(
        "X:/Saturated_flow/diffusion/2D_heterogeneous_profiles_" + ctype + ".png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
    )

# Checking total fixed microbes vs carrying capacity
Fbiomassvars = [Bfo1, Bfa1, Bfn1, Bifo1, Bifa1, Bifn1]
Fbiomassgvarnames = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Inactive fixed Aerobes",
    "Inactive fixed Ammonia oxidizers",
    "Inactive fixed Nitrate reducers",
]
for Reg in Regimes:
    d = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    for j in intsce:
        picture = sssp.biomassvscarryingcapacity(
            j,
            Het[Trial.index(j)],
            Anis[Trial.index(j)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            Fbiomassvars,
        )
        picture.savefig(
            "//tsclient/D/Saturated_flow/diffusion/"
            + Reg
            + str(j)
            + "_biomassvscarryingcapacity.png",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0,
        )
        picture.savefig(
            "//tsclient/D/Saturated_flow/diffusion/"
            + Reg
            + str(j)
            + "_biomassvscarryingcapacity.pdf",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0,
        )
start = 0
end = 15
chem = "DOC"

# Calculating average velocities
intsce = [
    "H",
    37,
    38,
    39,
    46,
    47,
    48,
    55,
    56,
    57,
    64,
    65,
    66,
    40,
    41,
    42,
    49,
    50,
    51,
    58,
    59,
    60,
    67,
    68,
    69,
    70,
    71,
    72,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    43,
    44,
    45,
    52,
    53,
    54,
    61,
    62,
    63,
    73,
    74,
    75,
]
for j in range(len(Trial)):
    df, massendtime, masstime, conctime, Velocity = sssa.calcconcmasstime(
        Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
    )
    print("Scenario:", Trial[j], ", Velocity: ", Velocity)

# Reading and storing in numpy array
for Reg in ["Slow"]:
    #    d = r"X:/Saturated_flow/diffusion/"+Reg+"AR_changedkindox/"
    d = r"X:/Saturated_flow/diffusion/" + Reg + "AR_changedkindox/"
    for j in range(len(Trial)):
        di = d + fpre + str(Trial[j]) + fsuf
        print(di)
        fwithd = di + filename
        print("Reading tec file....")
        size, steps, Headers, D = rdr.readTecfile(fwithd)
        print("Converting to array....")
        df = rdr.Converttomarr(D)
        print("Saving numpy array....")
        np.save(di + fpre + str(Trial[j]) + "_df", df)
        # Test for correct orientation of the data
        for i in range(np.shape(df)[0]):
            print(Headers[i + 3], np.mean(df[i, steps - 1, 0, :]))

# Final step outputs for Martin:
f = "X:/results_" + Reg + ".csv"
csvfile = open(f, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(["Y_cm", Headers[3:]])
for y in range(51):
    writer.writerow(
        [
            y,
            df[0, -1, y, 15],
            df[1, -1, y, 15],
            df[2, -1, y, 15],
            df[3, -1, y, 15],
            df[4, -1, y, 15],
            df[5, -1, y, 15],
            df[6, -1, y, 15],
            df[7, -1, y, 15],
            df[8, -1, y, 15],
            df[9, -1, y, 15],
            df[10, -1, y, 15],
            df[11, -1, y, 15],
            df[12, -1, y, 15],
            df[13, -1, y, 15],
            df[14, -1, y, 15],
            df[15, -1, y, 15],
            df[16, -1, y, 15],
            df[17, -1, y, 15],
            df[18, -1, y, 15],
            df[19, -1, y, 15],
            df[20, -1, y, 15],
            df[21, -1, y, 15],
            df[22, -1, y, 15],
            df[23, -1, y, 15],
            df[24, -1, y, 15],
            df[25, -1, y, 15],
            df[26, -1, y, 15],
        ]
    )
csvfile.close()

# Boxplots:

# tracer studies
tr1 = 8 - gw
vely = 5
velx = 4
vars = [tr1]
steps = [500 * 0.005, 2 * 0.005, 2 * 0.0005]
Regimes = ["Slow", "Equal", "Fast"]
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
sssp.plot_tracer("X:/Saturated_flow/Steady_state/Tracer_studies/tracer_combined_05032020.csv")
#or sssp.plot_tracer(f)


# plotting time series of tracer breakthrough in all flow regimes and trials:
ncols = 8
nrows = 6
for Reg, step in zip(Regimes, steps):
    d = r"X:/Saturated_flow/Steady_state/Tracer_studies/" + Reg + "AR/"
    fpre = "NS-A"
    fig, axes = plt.subplots(
        ncols=ncols, nrows=nrows, figsize=[15, 10], sharex=True, sharey=True
    )
    plt.suptitle("Tracer breakthrough curve" + Reg + " flow regime")
    col = 0
    for j in range(len(Trial)):
        df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
            Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
        )
        xindex = list(x * step for x in list(range(np.shape(conctime)[0])))
        axes.flat[j].plot(xindex, conctime[:, yout, 0])
        axes.flat[j].set_xlabel("")
        #        axes[colidx1][colidx2].set_xticklabels([])
        axes.flat[j].set_ylabel("")
        axes.flat[j].set_title(Trial[j])
        col = col + 1
    for ax in axes[:, 0]:
        ax.set_ylabel("Tracer (uM)")
    for ax in axes[5]:
        ax.set_xlabel("Time (days)")
    fig.subplots_adjust(left=0.15, top=0.9)
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

# Boxplots:
sgvarnames = ["TOC", "DO", "Nitrogen"]
# Mass flux with variance
dummy = sssp.boxV_Aflux(fastc, slowc, equalc, sgvarnames, imgsize=[24, 18])
picname = (
    "X:/Saturated_flow/diffusion/flux_onlyvariance_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname, dpi=300, bbox_inches="tight")
# Biomass with variance
dummy = sssp.boxV_Abio(fastb, slowb, equalb, imgsize=[22, 15])
picname = (
    "X:/Saturated_flow/Anaerobic/biomass_onlyvariance_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
# Biomass with mass flux
dummy = sssp.scatterbioflux(fastb, slowb, equalb)
picname = (
    "X:/Saturated_flow/Anaerobic/biomass_flux_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
# Mass flux wtih varxanis
dummy = sssp.boxVAflux(fastc, slowc, equalc, gvarnames)
picname = (
    "X:/Saturated_flow/Anaerobic/flux_VA_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
# Biomass wtih varxanis
dummy = sssp.boxVAbio(fastb, slowb, equalb)
picname = (
    "X:/Saturated_flow/Anaerobic/biomass_VA_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
# Flux wtih residence time
dummy = sssp.boxrestime_flux(fastc, slowc, equalc, gvarnames)
picname = (
    "X:/Saturated_flow/Anaerobic/flux_restime_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
# Flux wtih residence time
dummy = sssp.boxrestime_biomass(fastb, slowb, equalb)
picname = (
    "X:/Saturated_flow/Anaerobic/biomass_restime_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".png"
)
dummy.savefig(picname)
