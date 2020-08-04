# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:05:45 2019

@author: khurana
"""
##Don't use this. Use read_rate_data.py
import pandas as pd
import data_reader.data_processing as proc

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.colors import LogNorm
import math
from matplotlib.ticker import FuncFormatter


# Saturated flow regime
Reg = "Fast"
directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
fpre = "NS-A"
fsuf = r"/"
gw = 1
filename = "ratesAtFinish.dat"

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

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
listofcolumns = [6, 12]
for i in range(67):
    listofcolumns.append(i + 18)

ratenames = [
    "x_m",
    "Y",
    "Z",
    "Fixedaeroresp",
    "Mobaeroresp",
    "Fixedaerogwth",
    "Mobaerogwth",
    "Fixedactaerodett",
    "Fixedinaerodett",
    "Mobactaeroattach",
    "Mobinaeroattach",
    "FixeddeactlowDOX",
    "MobdeactlowDOX",
    "Fixedaeroreact",
    "Mobaeroreact",
    "Mortfixedactaero",
    "Mortmobactaero",
    "Mortinfixedaero",
    "Mortinmobaero",
    "Fixednitraresp",
    "Mobnitraresp",
    "Fixednitragwth",
    "Mobnitragwth",
    "Fixedactnitradett",
    "Fixedinnitradett",
    "Mobactnitraattach",
    "Mobinnitraattach",
    "FixeddeactlowN",
    "MobdeactlowN",
    "Fixednitrareact",
    "Mobnitrareact",
    "Mortfixedactnitra",
    "Mortmobactnitra",
    "Mortinfixednitra",
    "Mortinmobnitra",
    "Fixedsulpharesp",
    "Mobsulpharesp",
    "Fixedsulphagwth",
    "Mobsulphagwth",
    "Fixedactsulphadett",
    "Fixedinsulphadett",
    "Mobactsulphaattach",
    "Mobinsulphaattach",
    "FixedDeactlowS",
    "MobDeactlowS",
    "Fixedsulphareact",
    "Mobsulphareact",
    "Mortfixedactsulpha",
    "Mortmobactsulpha",
    "Mortinfixedsulpha",
    "Mortinmobsulpha",
    "Fixedammresp",
    "Mobammresp",
    "Fixedammgwth",
    "Mobammgwth",
    "Fixedactammdett",
    "Fixedinammdett",
    "Mobactammattach",
    "Mobinammattach",
    "FixedammdeactlowA",
    "MobammdeactlowA",
    "Fixedammreact",
    "Mobammreact",
    "Mortfixedactamm",
    "Mortmobactamm",
    "Mortinfixedamm",
    "Mortinmobamm",
    "Hydrolysis",
]

respindx = np.array(
    [
        ratenames.index("Fixedaerogwth"),
        ratenames.index("Mobaerogwth"),
        ratenames.index("Fixednitragwth"),
        ratenames.index("Mobnitragwth"),
        ratenames.index("Fixedsulphagwth"),
        ratenames.index("Mobsulphagwth"),
        ratenames.index("Fixedammgwth"),
        ratenames.index("Mobammgwth"),
    ]
)
respindx = np.array(
    [
        ratenames.index("Fixedaerogwth"),
        ratenames.index("Fixedammgwth"),
        ratenames.index("Fixednitragwth"),
    ]
)
respindx = np.array(
    [
        ratenames.index("Fixedaeroresp"),
        ratenames.index("Mobaeroresp"),
        ratenames.index("Fixednitraresp"),
        ratenames.index("Mobnitraresp"),
        ratenames.index("Fixedsulpharesp"),
        ratenames.index("Mobsulpharesp"),
        ratenames.index("Fixedammresp"),
        ratenames.index("Mobammresp"),
    ]
)
gratenames = [
    "Immobile aerobic respiration",
    "Mobile aerobic respiration",
    "Immobile nitrate respiration",
    "Mobile nitrate respiration",
    "Immobile sulphate respiration",
    "Mobile sulphate respiration",
    "Immobile ammonia respiration",
    "Mobile ammonia respiration",
]

gvarnames = ["DO", "Ammonium", "Nitrate"]

Regimes = ["Slow", "Equal", "Fast"]

AFbiomassvars = [
    Bfo1,
    Bfa1,
    Bfn1]

AFbiomassgvarnames = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
]

# Calculate bulk Damkohler numbers in the domain
respindx = np.array(range(66)[1:])

biomass_path_data = r"Y:\Home\khurana\4. Publications\Restructuring\Paper1\Figurecodes\biomass_withbreakthrough_forMartin_v4_complete.csv"
biomass = pd.read_csv(biomass_path_data, sep = '\t')
biomass.columns

microbes = ["Active fixed Aerobes", "Active fixed Ammonia oxidizers", "Active fixed Nitrate reducers"]

biomass = biomass[biomass['Chem'].isin (microbes)]
biomass.shape

row = []
for Reg in Regimes:
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    if Reg == "Equal":
        Reg = "Medium"
    for t in Trial:
        print(t)
        filepath = directory + fpre + str(t) + fsuf + filename
        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=18 + respindx)
        df = np.load(directory + fpre + t + fsuf + fpre + t + "_df.npy")
        for bioindx, bioname, i,c in zip(AFbiomassvars, microbes, range(len(respindx)), gvarnames):
            biorow = biomass.loc[(biomass.Regime == Reg) & (biomass.Trial == t) & (biomass.Chem == bioname)]
            meanrate = np.mean(M[:, i])
            meanbiomass = np.mean(df[bioindx, - 1, :, :])
            traveltime = biorow.Breakthroughtime.iloc[0]
            da = traveltime*meanrate
            dabio = traveltime*meanrate/meanbiomass
            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate, bioname, meanbiomass, meanrate/meanbiomass, da, dabio, c])

df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Rate_type', 'Meanrate', 'Microbe', 'Meanbiomassconc', 'Rateperbio', 'Da', 'Dabio', 'Chem'])


df.to_csv(r"Z:\Saturated_flow\diffusion_transient\rates_biomassconc_mean_ss.csv", sep = '\t')

respindx = np.array(
    [
        ratenames.index("Fixedaeroresp"),
        ratenames.index("Fixedammresp"),
        ratenames.index("Fixednitraresp")
    ]
)
gratenames = [
    "Immobile aerobic respiration",
    "Immobile ammonia respiration",
    "Immobile nitrate respiration",
]
gvarnames = ["DO", "Ammonium", "Nitrate"]

Regimes = ["Slow", "Equal", "Fast"]
velocities = [0.00038, 0.0038, 0.038]
# Calculate bulk Damkohler numbers in the domain

row = []
for vel, Reg in zip(velocities, Regimes):
    directory = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "AR_0/"
    if Reg == "Equal":
        Reg = "Medium"
    for t in Trial:
        print(t)
        filepath = directory + fpre + str(t) + fsuf + filename
        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=16 + respindx)
        for bio,i,c in zip(microbes, range(len(respindx)), gvarnames):
            biorow = biomass.loc[(biomass.Regime == Reg) & (biomass.Trial == t) & (biomass.Chem == bio)]
            biototal = biorow.Total_biomass.iloc[0]
            ratesum = np.mean(M[:, i])
            da = biorow.Breakthroughtime.iloc[0]*ratesum
            dabio = biorow.Breakthroughtime.iloc[0]*ratesum/(biototal)
            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], ratesum, bio, ratesum/biototal, da, dabio, c])

df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Rate_type', 'Totalrate', 'Microbe', 'Rateperbio', 'Da', 'Dabio', 'Chem'])


df.to_csv(r"Z:\Saturated_flow\diffusion_transient\rates_median_ss.csv", sep = '\t')

# Plot Damkohler numbers in the domain


for j in range(len(Trial)):
    #    j = Trial.index(63)
    di = d + fpre + str(Trial[j]) + fsuf
    print(str(Trial[j]))
    fwithd = di + filename
    M = np.loadtxt(fwithd, dtype=float, delimiter=" ", usecols=16 + respindx)
    for i in range(len(respindx)):
        summ[j, i] = sum(M[:, i])
        act[j, :, i] = M[:, i]
    Numberofrates = np.shape(M)[1]
    ratedf = np.ndarray([Numberofrates, 51, 31])
    ratedf2 = np.ndarray([Numberofrates, 51, 31])
    Da = np.ndarray([Numberofrates, 51, 31])
    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    for i in range(Numberofrates):
        for k in range(51):
            for l in range(31):
                ratedf2[i, k, l] = M[(k) * 31 + l, i]
    for k in range(51):
        ratedf[:, 50 - k, :] = ratedf2[:, k, :]
    df = np.load(di + fpre + str(Trial[j]) + "_df.npy")
    for resp in range(Numberofrates):
        Da[resp, :, :] = ratedf[resp, :, :] / abs(df[2, np.shape(df)[1] - 1, :, :])
    df2 = np.load(di + fpre + str(Trial[j]) + "_df.npy")
    vel = df2[2, np.shape(df)[1] - 1, :, :] * -1
    fig, axes = plt.subplots(
        nrows=1, ncols=len(gratenames) + 1, figsize=(26, 10), sharex=True, sharey=True
    )
    plt.suptitle(
        "Variance " + str(Het[j]) + " : Anisotropy " + str(Anis[j]),
        ha="center",
        va="center",
        fontsize=30,
    )
    for i in range(len(gratenames)):
        log_norm = LogNorm(vmin=Da[i, :, :].min().min(), vmax=Da[i, :, :].max().max())
        sns.set(font_scale=2)
        cbar_ticks = [
            math.pow(10, x)
            for x in range(
                int(math.floor(math.log10(Da[i, :, :].min().min()))),
                1 + int(math.ceil(math.log10(Da[i, :, :].max().max()))),
            )
        ]
        name = sns.heatmap(
            Da[i, :, :],
            square=False,
            norm=log_norm,
            cmap="YlGnBu",
            cbar_kws={"ticks": cbar_ticks},
            vmin=Da[i, :, :].min().min(),
            vmax=Da[i, :, :].max().max(),
            xticklabels=False,
            yticklabels=False,
            ax=axes.flat[i],
        )
        axes.flat[i].set_title(gratenames[i], fontsize=30)
    v = sns.heatmap(
        vel,
        square=False,
        norm=LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()),
        cmap="YlGnBu",
        xticklabels=False,
        yticklabels=False,
        ax=axes.flat[3],
        cbar_kws={
            "format": FuncFormatter(fmtvel),
            "ticks": [
                math.pow(10, x)
                for x in range(
                    int(math.floor(math.log10(np.min(abs(vel))))),
                    1 + int(math.ceil(math.log10(np.max(abs(vel))))),
                )
            ],
        },
        vmin=abs(vel).min().min(),
        vmax=abs(vel).max().max(),
    )
    axes.flat[3].set_title("Velocity", fontsize=0)
    picname = (
        d
        + fpre
        + str(Trial[j])
        + "_"
        + "heatmap"
        + datetime.now().strftime("%d%m%Y%H%M")
        + "_LOG_Da.png"
    )
    pictosave = name.get_figure()
    pictosave.savefig(picname)

# calculating node wise difference from homoegeneous scenario for the activities/respiration rates
actarr = sk.RateConverttomarr(act)
actarrscaled = np.zeros(np.shape(actarr))
for j in range(len(Trial)):
    for k in range(np.shape(actarrscaled)[0]):
        minrate = np.min(actarr[k, j, :, :])
        maxrate = np.max(actarr[k, j, :, :])
        actarrscaled[k, j, :, :] = (actarr[k, j, :, :] - minrate) / (maxrate - minrate)

delresprate = np.zeros([len(Trial), 1581, len(respindx)])
delsummresprate = np.zeros([len(Trial), len(respindx)])
for j in range(len(Trial)):
    delresprate[j, :, :] = (act[0, :, :] - act[j, :, :]) / act[0, :, :]
    delsummresprate[j, :] = (summ[0, :] - summ[j, :]) / summ[0, :]
#    sk.heatmapratedist(delresprate[j,:,:], respindx)
delrespratearr = sk.RateConverttomarr(delresprate)
delrespratearrlog = np.zeros(np.shape(delrespratearr))
delrespratearrlog = np.zeros([4, 49, 51, 31])
for j in range(len(Trial)):
    for k in range(np.shape(delrespratearr)[0]):
        for l in range(np.shape(delrespratearr)[2]):
            for m in range(np.shape(delrespratearr)[3]):
                if delrespratearr[k, j, l, m] == 0:
                    pass
                else:
                    delrespratearrlog[k, j, l, m] = np.log(delrespratearr[k, j, l, m])

delrespratearrscaled = np.zeros([4, 49, 51, 31])
for j in range(len(Trial) - 1):
    for k in range(np.shape(delrespratearr)[0]):
        minrate = np.min(delrespratearr[k, j + 1, :, :])
        maxrate = np.max(delrespratearr[k, j + 1, :, :])
        delrespratearrscaled[k, j, :, :] = (delrespratearr[k, j, :, :] - minrate) / (
            maxrate - minrate
        )

intrespindx = np.array([3, 19, 51])
nrows = 1
ncols = 3
figsize = [10, 5]
for j in range(len(Trial)):
    print(str(Trial[j]))
    name = sk.heatmapratedistLOG(
        delrespratearr[:, j, :, :], intrespindx, gratenames, nrows, ncols, figsize
    )
    picname = (
        d
        + fpre
        + str(Trial[j])
        + datetime.now().strftime("%d%m%Y%H%M")
        + "_delrespratesatss_log.png"
    )
    name.get_figure().savefig(picname)

# Writing the results in a csv file
fname = (
    d
    + fpre
    + str(Trial[0])
    + "_"
    + str(Trial[-1])
    + "_"
    + "Resp_"
    + datetime.now().strftime("%d%m%Y%H%M")
    + ".csv"
)
csvfile = open(fname, "w")
writer = csv.writer(
    csvfile,
    delimiter="\t",
    quotechar="\t",
    quoting=csv.QUOTE_MINIMAL,
    lineterminator="\n",
)
writer.writerow(
    ["Sno", "Trial", "Het", "Anis", "Rate_type", "Total_rate", "Change_in_rate"]
)
for j in range(len(Trial)):
    for k in range(len(respindx)):
        writer.writerow(
            [
                j,
                Trial[j],
                Het[j],
                Anis[j],
                gratenames[k],
                summ[j, k],
                delsummresprate[j, k],
            ]
        )
csvfile.close()
print("File written")

# Checking transfer of velocity function - it looks good
const = 0.0842 * (0.0089 ** 0.58)
for j in range(len(Trial)):
    j = Trial.index(63)
    di = d + fpre + str(Trial[j]) + fsuf
    print(str(Trial[j]))
    fwithd = di + filename
    M = np.loadtxt(fwithd, dtype=float, delimiter=" ", usecols=16 + respindx)
    for i in range(len(respindx)):
        summ[j, i] = sum(M[:, i])
        act[j, :, i] = M[:, i]
    Numberofrates = np.shape(M)[1]
    ratedf = np.ndarray([Numberofrates, 51, 31])
    ratedf2 = np.ndarray([Numberofrates, 51, 31])
    Da = np.ndarray([Numberofrates, 51, 31])
    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    for i in range(Numberofrates):
        for k in range(51):
            for l in range(31):
                ratedf2[i, k, l] = M[(k) * 31 + l, i]
    for k in range(51):
        ratedf[:, 50 - k, :] = ratedf2[:, k, :]
    df = np.load(di + fpre + str(Trial[j]) + "_df.npy")
    for resp in range(Numberofrates):
        totalsum = (
            df[Bfo1, np.shape(df)[1] - 1, :, :]
            + df[Bifo1, np.shape(df)[1] - 1, :, :]
            + df[Bfn1, np.shape(df)[1] - 1, :, :]
            + df[Bifn1, np.shape(df)[1] - 1, :, :]
            + df[Bfa1, np.shape(df)[1] - 1, :, :]
            + df[Bifa1, np.shape(df)[1] - 1, :, :]
            + df[Bfs1, np.shape(df)[1] - 1, :, :]
            + df[Bifs1, np.shape(df)[1] - 1, :, :]
        )
        pop = df[vars[resp], np.shape(df)[1] - 1, :, :]
        subt = 1 / (np.exp((500 - totalsum) / (0.1 * 500)) + 1)
        Da[resp, :, :] = ((ratedf[resp, :, :] / pop - subt) / const) ** (1 / 0.58)
#        Da[resp,:,:] = (((ratedf[resp,:,:]/pop) - subt)/const)**(1/0.58)
