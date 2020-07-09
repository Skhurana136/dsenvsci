# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:24:47 2020

@author: khurana
"""
import pandas as pd
import data_reader.data_processing as proc
import numpy as np

# Saturated flow regime
Reg = "Fast"
fpre = "NS-A"
fsuf = r"/"
gw = 1
filename = "ratesAtFinish.dat"

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
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

# Default:
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())

Bfo1 = 8 - gw
Bfn1 = 15 - gw
Bfa1 = 25 - gw
Bmo1 = 9 - gw
Bmn1 = 16 - gw
Bma1 = 26 - gw

AFbiomassvars = [
    Bfo1,
    Bfa1,
    Bfn1,
    Bmo1,
    Bma1,
    Bmn1,
]

AFbiomassgvarnames = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
]

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

microvars = [
    Bfo1,
    Bfa1,
    Bfn1]

microbes = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
]
# Calculate bulk Damkohler numbers in the domain

biomass_path_data = r"Y:\Home\khurana\4. Publications\Restructuring\Paper1\Figurecodes\biomass_withbreakthrough_forMartin_v4_complete.csv"
biomass = pd.read_csv(biomass_path_data, sep = '\t')
biomass.columns
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
        M = np.loadtxt(filepath, dtype=float, delimiter=" ", usecols=16 + respindx)
        df = np.load(directory + fpre + t + fsuf + fpre + t + "_df.npy")
        for bioindx, bioname, i,c in zip(microvars, microbes, range(len(respindx)), gvarnames):
            biorow = biomass.loc[(biomass.Regime == Reg) & (biomass.Trial == t) & (biomass.Chem == bioname)]
            meanrate = np.mean(M[:, i])
            meanbiomass = np.mean(df[bioindx, - 1, :, :])
            traveltime = biorow.Breakthroughtime.iloc[0]
            da = traveltime*meanrate
            dabio = traveltime*meanrate/meanbiomass
            ka = traveltime*meanrate
            kabio = traveltime*meanrate/meanbiomass
            row.append([Reg, t, scdict[t]['Het'], scdict[t]['Anis'], gratenames[i], meanrate, bioname, meanbiomass, meanrate/meanbiomass, da, ka, dabio, kabio, c])

df = pd.DataFrame.from_records(row, columns = ['Regime', 'Trial', 'Variance', 'Anisotropy', 'Rate_type', 'Meanrate', 'Microbe', 'Meanbiomassconc', 'Rateperbio', 'Da', 'Ka', 'Dabio', 'Kabio', 'Chem'])


df.to_csv(r"Z:\Saturated_flow\diffusion_transient\Da_mean_ss.csv", sep = '\t')