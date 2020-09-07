# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:55:52 2020

@author: khurana
"""
import numpy as np
from analyses.saturated_transient import calcconcmasstime

def calcmassfluxnew(numpyarray, yin, yout, xleft, xright, gvarnames, flowregime):
    import data_reader.data_processing as proc
    reactivespecies = proc.masterdissolvedspecies(flowregime)
    microbialspecies = proc.mastermicrobialspecies(flowregime)
    mobilespecies = list(t for t in microbialspecies.keys() if microbialspecies[t]['Location'] == "Mobile")
    species = {**reactivespecies, **microbialspecies}

    vedge = 0.005
    velem = 0.01
    por = 0.2
    
    massfluxin = np.zeros([len(gvarnames)])
    massfluxout = np.zeros([len(gvarnames)])
    df = numpyarray
    veliredg = df[2, -1, yin, xright]
    veliledg = df[2, -1, yin, xleft]
    veloredg = df[2, -1, yout, xright]
    veloledg = df[2, -1, yout, xleft]
    veloelem = df[2, -1, yout, xleft + 1 : xright]
    velielem = df[2, -1, yin, xleft + 1 : xright]
#    velelem = df[2, -1, yin + 1 : yout, xleft + 1 : xright]
#    vellelem = df[2, -1, yin + 1 : yout, xleft]
#    velrelem = df[2, -1, yin + 1 : yout, xright]
    if flowregime == "Saturated":
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
    elif flowregime == "Unsaturated":
        satiredg = df[4, -1, yin, xright]
        satiledg = df[4, -1, yin, xleft]
        satoredg = df[4, -1, yout, xright]
        satoledg = df[4, -1, yout, xleft]
        satoelem = df[4, -1, yout, xleft + 1 : xright]
        satielem = df[4, -1, yin, xleft + 1 : xright]
#        satlelem = df[4, 1:, yin + 1 : yout, xleft]
#        satrelem = df[4, 1:, yin + 1 : yout, xright]
#        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in gvarnames:
        idx = gvarnames.index(i)
        if i == "Nitrogen":
            Nspecies = mobilespecies + ["Particulate organic matter"]
            ninlet = 0
            noutlet = 0
            for n in Nspecies:
                 ninlet = (
                    ninlet
                    + (
                        df[species[n]['TecIndex'], -1, yin, xleft] * satiledg * veliledg * vedge
                        + df[species[n]['TecIndex'], -1, yin, xright] * satiredg * veliredg * vedge
                        + sum(
                            df[species[n]['TecIndex'], -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem
                        )
                    )
                    / por
                    )
                 noutlet = (
                    noutlet
                    + (
                        df[species[n]['TecIndex'], -1, yout, xleft] * satoledg * veloledg * vedge
                        + df[species[n]['TecIndex'], -1, yout, xright] * satoredg * veloredg * vedge
                        + sum(
                            df[species[n]['TecIndex'], -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem
                        )
                    )
                    / por
                )
            sumin = 0
            sumout = 0
            for r in ["Ammonium", "Nitrate"]:                            
                rin = (
                        df[reactivespecies[r]['TecIndex'], -1, yin, xleft]*satiledg*veliledg*vedge
                + df[reactivespecies[r]['TecIndex'], -1, yin, xright]*satiredg*veliredg*vedge
                + np.sum(df[reactivespecies[r]['TecIndex'], -1, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1))/por
                rout = (
                        df[reactivespecies[r]['TecIndex'], -1, yout, xleft]*satoledg*veloledg*vedge
                + df[reactivespecies[r]['TecIndex'], -1, yout, xright]*satoredg*veloredg*vedge
                + np.sum(
                        df[reactivespecies[r]['TecIndex'], -1, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                        )
                ) / por
                sumin = sumin + rin
                sumout = sumout + rout
            massfluxin[idx] = ninlet/10 + sumin
            massfluxout[idx] = noutlet/10 + sumout
        elif i == "TOC":
            cinlet = 0
            coutlet = 0
            for c in list(mobilespecies + ["DOC"] + ["Particulate organic matter"]):
                cinlet = (
                        cinlet
                        + (
                        df[species[c]['TecIndex'], -1, yin, xleft]
                            * satiledg
                            * veliledg
                            * vedge
                        
                        + df[species[c]['TecIndex'], -1, yin, xright]
                                * satiredg
                                * veliredg
                                * vedge
                            
                        + np.sum(
                                    df[species[c]['TecIndex'], -1, yin, xleft + 1 : xright]
                                    * satielem
                                    * velielem
                                    * velem,
                                    axis=-1,
                                )
                            )/por
                        )
                    
                coutlet = (
                    coutlet
                    + (df[species[c]['TecIndex'], -1, yout, xleft]
                            * satoledg
                            * veloledg
                            * vedge
                            
                        + df[species[c]['TecIndex'], -1, yout, xright]
                                * satoredg
                                * veloredg
                                * vedge
                            
                        + np.sum(
                                    df[species[c]['TecIndex'], -1, yout, xleft + 1 : xright]
                                    * satoelem
                                    * veloelem
                                    * velem,
                                    axis=-1,
                                )
                            )/por
                        )
                        
            massfluxin[idx] = cinlet
            massfluxout[idx] = coutlet
        else:
            massfluxin[idx] = (
                df[species[i]['TecIndex'], -1, yin, xleft]
                        * satiledg
                        * veliledg
                        * vedge
                   
                +  df[species[i]['TecIndex'], -1, yin, xright]
                        * satiredg
                        * veliredg
                        * vedge
                   
                + np.sum(
                            df[species[i]['TecIndex'], -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem,
                            axis=-1,
                        )
                    ) / por
            massfluxout[idx] = (                    df[species[i]['TecIndex'], -1, yout, xleft]
                        * satoledg
                        * veloledg
                        * vedge
                   
                    + df[species[i]['TecIndex'], -1, yout, xright]
                        * satoredg
                        * veloredg
                        * vedge
                    + np.sum(                            df[species[i]['TecIndex'], -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem,
                                axis=-1,
                        )
                    ) / por
    return massfluxin, massfluxout

    
def calcmassflux(
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
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    por = 0.2
    doc1 = 10 - gw
    Bmo1 = 9 - gw
    Bmn1 = 16 - gw
    Bms1 = 21 - gw
    Bma1 = 26 - gw
    Bimo1 = 14 - gw
    Bimn1 = 19 - gw
    Bims1 = 24 - gw
    Bima1 = 28 - gw
    POM1 = 30 - gw
    Amm1 = 12 - gw
    nitra1 = 17 - gw
    Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    mf = np.zeros([len(Trial) * len(gvarnames), 9])
    for j in range(len(Trial)):
        di = directory + fpre + str(Trial[j]) + fsuf
        print(str(Trial[j]))
        #        di+fpre+str(Tforfpre[k])+str(Trial[j])+'_df'
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
            idx = j * len(gvarnames) + i
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
                        + (
                            df[n - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                            + df[n - 3, -1, yin, xright] * satiredg * veliredg * vedge
                            + sum(
                                df[n - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                        / por
                    )
                    noutlet = (
                        noutlet
                        + (
                            df[n - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                            + df[n - 3, -1, yout, xright] * satoredg * veloredg * vedge
                            + sum(
                                df[n - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                        / por
                    )
                ninlet = (
                    ninlet / 10
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
                    / por
                )
                noutlet = (
                    noutlet / 10
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
                    / por
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
                        + (
                            df[c - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                            + df[c - 3, -1, yin, xright] * satiredg * veliredg * vedge
                            + sum(
                                df[c - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                        / por
                    )
                    coutlet = (
                        coutlet
                        + (
                            df[c - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                            + df[c - 3, -1, yout, xright] * satoredg * veloredg * vedge
                            + sum(
                                df[c - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                        / por
                    )
                mf[idx, 4] = abs(
                    cinlet
                )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx, 5] = abs(
                    coutlet
                )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            else:
                mf[idx, 4] = abs(
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
                    / por
                )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx, 5] = abs(
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
                    / por
                )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            # calculating removal in mass
            mf[idx, 6] = mf[idx, 4] - mf[idx, 5]
            # normalizing removal in mass with incoming mass
            mf[idx, 7] = (mf[idx, 4] - mf[idx, 5]) / mf[idx, 4]
            # comparing removal with homogeneous case
            mf[idx, 8] = mf[idx, 6] / mf[i, 6]
    return mf


def calcsum(data,yin,yout,xleft,xright,gvarnames,flowregime):
    import data_reader.data_processing as proc
    species = proc.mastermicrobialspecies("Saturated")
    vedge = 0.005
    velem = 0.01
    por = 0.2
    df = data
    #    sumall = np.zeros([len(Trial)*len(vars),9])
    sumall = np.zeros([len(gvarnames)])
    if flowregime == "Saturated":
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
    for g in gvarnames:
        idx = gvarnames.index(g)
        sumall[idx] = (
                (
                    (
                        df[species[g]['TecIndex'], -1, yin, xleft] * satiledg
                        + df[species[g]['TecIndex'], -1, yout, xleft] * satoledg
                        + df[species[g]['TecIndex'], -1, yin, xright] * satiredg
                        + df[species[g]['TecIndex'], -1, yout, xright] * satoredg
                    )
                    * vedge
                    * vedge
                    + (
                        sum(
                            df[species[g]['TecIndex'],
                                -1,
                                yin,
                                xleft + 1 : xright,
                            ]
                            * satielem
                        )
                        + sum(df[species[g]['TecIndex'],
                                -1,
                                yout,
                                xleft + 1 : xright,
                            ]
                            * satoelem
                        )
                        + sum(
                            df[species[g]['TecIndex'], -1, yin + 1 : yout, xleft]
                            * satlelem
                        )
                        + sum(df[species[g]['TecIndex'], -1, yin + 1 : yout, xright]
                            * satrelem
                        )
                    )
                    * vedge
                    * velem
                )
                + sum(
                    sum(
                        df[species[g]['TecIndex'],
                            -1,
                            yin + 1 : yout,
                            xleft + 1 : xright,
                        ]
                        * satelem
                    )
                )
                * velem
                * velem
            ) * por
    return sumall


def calcaggres(mf):
    H = sorted(np.unique(mf[:, 1, ...]))
    A = sorted(np.unique(mf[..., 2]))
    C = sorted(np.unique(mf[..., 3]))
    meanresults = []
    stdresults = []
    covresults = []
    for i in H:
        remH = mf[np.where(mf[..., 1] == i)]
        for j in A:
            remHA = remH[np.where(remH[..., 2] == j)]
            for k in C:
                meanresults.append(
                    [
                        np.mean(remHA[np.where(remHA[..., 3] == k)][..., 6]),
                        np.mean(remHA[np.where(remHA[..., 3] == k)][..., 7]),
                        i,
                        j,
                        k,
                    ]
                )
                stdresults.append(
                    [
                        np.std(remHA[np.where(remHA[..., 3] == k)][..., 6]),
                        np.std(remHA[np.where(remHA[..., 3] == k)][..., 7]),
                        i,
                        j,
                        k,
                    ]
                )
                covresults.append(
                    [
                        np.cov(remHA[np.where(remHA[..., 3] == k)][..., 7]),
                        np.cov(remHA[np.where(remHA[..., 3] == k)][..., 7]),
                        i,
                        j,
                        k,
                    ]
                )
    cleanmresults = [x for x in meanresults if str(x[0]) != "nan"]
    cleansresults = [x for x in stdresults if str(x[0]) != "nan"]
    cleancresults = [x for x in covresults if str(x[0]) != "nan"]
    meanresultsarr = np.array(cleanmresults)
    stdresultsarr = np.array(cleansresults)
    covresultsarr = np.array(cleancresults)
    return meanresultsarr, stdresultsarr, covresultsarr


def calcoxiccells(
    limit,
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
):
    oxiccells = np.zeros([len(Trial), 51, 1])
    for j in range(len(Trial)):
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial[j],
            Het[j],
            Anis[j],
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
        c = []
        for k in range(51):
            c.append(
                np.count_nonzero(
                    df[vars[gvarnames.index("DO")] - 3, np.shape(df)[1] - 1, k, :]
                    > limit
                )
            )
        #        print (Trial[j], c)
        #        print ("Calculating number of oxic cells: ", c)
        oxiccells[j, :, 0] = c

    return oxiccells


# ATTENTION!!!!
# This function needs to be checked:
# def removalperbiomass():
#    removal = pd.read_csv("X:/massflux_withbreakthrough_revisedwithmobilebiomass.csv", delimiter="\t")
#    Biomass = pd.read_csv("X:/biomass_withbreakthrough.csv", delimiter = "\t")
#    newseries = removal[removal["Chem"]=="DO"]
#    newseries['Rem'] = 100*(newseries.loc[:,"Inlet_massflux"] - newseries.loc[:,"Outlet_massflux"])/newseries.loc[:,"Inlet_massflux"]
#    newb1series = Biomass[Biomass["Chem"]=="Active fixed Aerobes"]
#    newb2series = Biomass[Biomass["Chem"]=="Active mobile Aerobes"]
#    newbseries = newb1series.loc[:,'Total_biomass'].reset_index().add(newb2series.loc[:,'Total_biomass'].reset_index())
#    num = newseries.loc[:,'Rem'].reset_index().to_numpy
#    b = newbseries.loc[:,'Total_biomass'].reset_index().to_numpy
#    for i in range(len(num)):
#        num[:,'Rem'][i]/b[i]
#    newseries.loc[:,'Rem'].reset_index().div(newbseries.loc[:,'Total_biomass'].reset_index())

#    num.div(newbseries.loc[:,'Total_biomass'].reset_index())
#    for i in range(len(newseries)):
#        newseries.loc[:,'Rem'].reset_index().div(newbseries.loc[:,'Total_biomass'].reset_index())
#    return mf
#
