"""
# -*- coding: utf-8 -*-
# ======================================================================
# Created by : Swamini Khurana
# Created on : On Tues Dec 07 2021 at 00:32:06
# ======================================================================
# __author__ = Swamini Khurana
# __copyright__ = Copyright (c) 2021, Swamini Khurana, dsenvsci
# __credits__ = [Swamini Khurana]
# __license__ = MIT
# __version__ = 0.0.1
# __maintainer__ = Swamini Khurana
# __email__ = swamini.khurana@gmail.com
# __status__ = development
# ======================================================================
The file has been build for providing with configurations
about the reaction network solver
""" 

import numpy as np
from DS.analyses.transient import conc_time
import DS.data_reader.data_processing as proc

def massflux(numpyarray, yin, yout, xleft, xright, gvarnames, flowregime, vedge, velem, por):
    species = proc.speciesdict(flowregime)
    mobilespecies = list(t for t in species.keys() if (species[t]['Location'] == "Mobile") and (species[t]['State']!= "Dissolved"))
    
    massfluxin = np.zeros([len(gvarnames)])
    massfluxout = np.zeros([len(gvarnames)])
    df = numpyarray
    veliredg = df[2, -1, yin, xright]
    veliledg = df[2, -1, yin, xleft]
    veloredg = df[2, -1, yout, xright]
    veloledg = df[2, -1, yout, xleft]
    veloelem = df[2, -1, yout, xleft + 1 : xright]
    velielem = df[2, -1, yin, xleft + 1 : xright]
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
        def effsat(data):
            slope = 1/(0.8-0.2)
            constant = -1/3
            sat = slope*data + constant
            return sat
        satiredg = effsat(df[4, -1, yin, xright])
        satiledg = effsat(df[4, -1, yin, xleft])
        satoredg = effsat(df[4, -1, yout, xright])
        satoledg = effsat(df[4, -1, yout, xleft])
        satoelem = effsat(df[4, -1, yout, xleft + 1 : xright])
        satielem = effsat(df[4, -1, yin, xleft + 1 : xright])
    for i in gvarnames:
        idx = gvarnames.index(i)
        if i == "Nitrogen":
            Nspecies = mobilespecies
            ninlet = 0
            noutlet = 0
            for n in Nspecies:
                 ninlet = ninlet+(df[species[n]['TecIndex'], -1, yin, xleft]*satiledg*veliledg*vedge+
                            df[species[n]['TecIndex'], -1, yin, xright]*satiredg*veliredg*vedge+
                            sum(df[species[n]['TecIndex'], -1, yin, xleft + 1 : xright]*satielem*velielem*velem))/por
                 noutlet = noutlet+(df[species[n]['TecIndex'], -1, yout, xleft]*satoledg*veloledg*vedge+
                             df[species[n]['TecIndex'], -1, yout, xright]*satoredg*veloredg*vedge+
                        sum(df[species[n]['TecIndex'], -1, yout, xleft + 1 : xright]*satoelem*veloelem*velem))/por
            sumin = 0
            sumout = 0
            for r in ["Ammonium", "Nitrate"]:                            
                rin = (df[species[r]['TecIndex'], -1, yin, xleft]*satiledg*veliledg*vedge+
                       df[species[r]['TecIndex'], -1, yin, xright]*satiredg*veliredg*vedge+
                       np.sum(df[species[r]['TecIndex'], -1, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1))/por
                rout = (df[species[r]['TecIndex'], -1, yout, xleft]*satoledg*veloledg*vedge+
                        df[species[r]['TecIndex'], -1, yout, xright]*satoredg*veloredg*vedge+
                        np.sum(df[species[r]['TecIndex'], -1, yout, xleft + 1 : xright]*satoelem*veloelem*velem,axis=-1))/por
                sumin = sumin + rin
                sumout = sumout + rout
            massfluxin[idx] = ninlet/10 + sumin
            massfluxout[idx] = noutlet/10 + sumout
        elif i == "TOC":
            cinlet = 0
            coutlet = 0
            for c in list(mobilespecies + ["DOC"]):
                cinlet = cinlet+(df[species[c]['TecIndex'], -1, yin, xleft]*satiledg*veliledg*vedge+
                df[species[c]['TecIndex'], -1, yin, xright]*satiredg*veliredg*vedge+
                np.sum(df[species[c]['TecIndex'], -1, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1,))/por                    
                coutlet = coutlet+(df[species[c]['TecIndex'], -1, yout, xleft]*satoledg*veloledg*vedge+
                 df[species[c]['TecIndex'], -1, yout, xright]*satoredg*veloredg*vedge+
                 np.sum(df[species[c]['TecIndex'], -1, yout, xleft + 1 : xright]*satoelem*veloelem*velem,axis=-1,))/por
            massfluxin[idx] = cinlet
            massfluxout[idx] = coutlet
        else:
            massfluxin[idx] = (df[species[i]['TecIndex'], -1, yin, xleft]*satiledg*veliledg*vedge+
                               df[species[i]['TecIndex'], -1, yin, xright]*satiredg*veliredg*vedge+
                               np.sum(df[species[i]['TecIndex'], -1, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1,))/por
            massfluxout[idx] = (df[species[i]['TecIndex'], -1, yout, xleft]*satoledg*veloledg*vedge+
                                df[species[i]['TecIndex'], -1, yout, xright]*satoredg*veloredg*vedge+
                                np.sum(df[species[i]['TecIndex'], -1, yout, xleft + 1 : xright]*satoelem*veloelem*velem,axis=-1,))/por
    return massfluxin, massfluxout
    
def sum_biomass(data,yin,yout,xleft,xright,gvarnames,flowregime):
    species = proc.speciesdict(flowregime)
    vedge = 0.005
    velem = 0.01
    por = 0.2
    vbc = 0.5 * 0.3
    df = data
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
        def effsat(data):
            slope = 1/(0.8-0.2)
            constant = -1/3
            sat = slope*data + constant
            return sat
        satiredg = effsat(df[4, -1, yin, xright])
        satiledg = effsat(df[4, -1, yin, xleft])
        satoredg = effsat(df[4, -1, yout, xright])
        satoledg = effsat(df[4, -1, yout, xleft])
        satoelem = effsat(df[4, -1, yout, xleft + 1 : xright])
        satielem = effsat(df[4, -1, yin, xleft + 1 : xright])
        satlelem = effsat(df[4, -1, yin + 1 : yout, xleft])
        satrelem = effsat(df[4, -1, yin + 1 : yout, xright])
        satelem = effsat(df[4, -1, yin + 1 : yout, xleft + 1 : xright])
    for g in gvarnames:
        idx = gvarnames.index(g)
        sumall[idx] = (((df[species[g]['TecIndex'], -1, yin, xleft]*satiledg
                        + df[species[g]['TecIndex'], -1, yout, xleft]*satoledg
                        + df[species[g]['TecIndex'], -1, yin, xright]*satiredg
                        + df[species[g]['TecIndex'], -1, yout, xright]*satoredg)*vedge**2
                    + (sum(df[species[g]['TecIndex'],-1,yin,xleft + 1 : xright]*satielem)
                        + sum(df[species[g]['TecIndex'],-1,yout,xleft + 1 : xright]*satoelem)
                        + sum(df[species[g]['TecIndex'], -1, yin + 1 : yout, xleft]*satlelem)
                        + sum(df[species[g]['TecIndex'], -1, yin + 1 : yout, xright]*satrelem))*vedge*velem)
                    + sum(sum(df[species[g]['TecIndex'],-1,yin + 1 : yout,xleft + 1 : xright]*satelem))*velem**2)*por/vbc
    return sumall

def chem_stock(data,yin,yout,xleft,xright,gvarnames,flowregime, vedge, velem, por):
    species = proc.speciesdict(flowregime)
    mobilespecies = list(t for t in species.keys() if (species[t]['Location'] == "Mobile") and (species[t]['State']!= "Dissolved"))

    species = proc.speciesdict(flowregime)

    df = data
    stock_domain = np.zeros([len(gvarnames)])
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
        def effsat(data):
            slope = 1/(0.8-0.2)
            constant = -1/3
            sat = slope*data + constant
            return sat
        satiredg = effsat(df[4, -1, yin, xright])
        satiledg = effsat(df[4, -1, yin, xleft])
        satoredg = effsat(df[4, -1, yout, xright])
        satoledg = effsat(df[4, -1, yout, xleft])
        satoelem = effsat(df[4, -1, yout, xleft + 1 : xright])
        satielem = effsat(df[4, -1, yin, xleft + 1 : xright])
        satlelem = effsat(df[4, -1, yin + 1 : yout, xleft])
        satrelem = effsat(df[4, -1, yin + 1 : yout, xright])
        satelem = effsat(df[4, -1, yin + 1 : yout, xleft + 1 : xright])
    for g in gvarnames:
        idx = gvarnames.index(g)
        if g == "Nitrogen":
            Nspecies = mobilespecies
            n_microbes = 0
            for n in Nspecies:
                n_microbes = n_microbes+(((df[species[n]['TecIndex'], -1, yin, xleft]*satiledg
                                   + df[species[n]['TecIndex'], -1, yout, xleft]*satoledg
                                   + df[species[n]['TecIndex'], -1, yin, xright]*satiredg
                                   + df[species[n]['TecIndex'], -1, yout, xright]*satoredg)*vedge**2
                                  + (sum(df[species[n]['TecIndex'],-1,yin,xleft + 1 : xright]*satielem)
                                     + sum(df[species[n]['TecIndex'],-1,yout,xleft + 1 : xright]*satoelem)
                                     + sum(df[species[n]['TecIndex'], -1, yin + 1 : yout, xleft]*satlelem)
                                     + sum(df[species[n]['TecIndex'], -1, yin + 1 : yout, xright]*satrelem))*vedge*velem)
                                 + sum(sum(df[species[n]['TecIndex'],-1,yin + 1 : yout,xleft + 1 : xright]*satelem))*velem**2)*por
            pure_n = 0
            for r in ["Ammonium", "Nitrate"]:                            
                pure_n = pure_n + (((df[species[r]['TecIndex'], -1, yin, xleft]*satiledg
                                   + df[species[r]['TecIndex'], -1, yout, xleft]*satoledg
                                   + df[species[r]['TecIndex'], -1, yin, xright]*satiredg
                                   + df[species[r]['TecIndex'], -1, yout, xright]*satoredg)*vedge**2
                                  + (sum(df[species[r]['TecIndex'],-1,yin,xleft + 1 : xright]*satielem)
                                     + sum(df[species[r]['TecIndex'],-1,yout,xleft + 1 : xright]*satoelem)
                                     + sum(df[species[r]['TecIndex'], -1, yin + 1 : yout, xleft]*satlelem)
                                     + sum(df[species[r]['TecIndex'], -1, yin + 1 : yout, xright]*satrelem))*vedge*velem)
                                 + sum(sum(df[species[r]['TecIndex'],-1,yin + 1 : yout,xleft + 1 : xright]*satelem))*velem**2)*por
            stock_domain[idx] = n_microbes/10 + pure_n
        elif g == "TOC":
            c_domain = 0
            for c in list(mobilespecies + ["DOC"]):
                c_domain = c_domain+(((df[species[c]['TecIndex'], -1, yin, xleft]*satiledg
                                   + df[species[c]['TecIndex'], -1, yout, xleft]*satoledg
                                   + df[species[c]['TecIndex'], -1, yin, xright]*satiredg
                                   + df[species[c]['TecIndex'], -1, yout, xright]*satoredg)*vedge**2
                                  + (sum(df[species[c]['TecIndex'],-1,yin,xleft + 1 : xright]*satielem)
                                     + sum(df[species[c]['TecIndex'],-1,yout,xleft + 1 : xright]*satoelem)
                                     + sum(df[species[c]['TecIndex'], -1, yin + 1 : yout, xleft]*satlelem)
                                     + sum(df[species[c]['TecIndex'], -1, yin + 1 : yout, xright]*satrelem))*vedge*velem)
                                 + sum(sum(df[species[c]['TecIndex'],-1,yin + 1 : yout,xleft + 1 : xright]*satelem))*velem**2)*por
            stock_domain[idx] = c_domain
        else:
            stock_domain[idx] = (((df[species[g]['TecIndex'], -1, yin, xleft]*satiledg
                                   + df[species[g]['TecIndex'], -1, yout, xleft]*satoledg
                                   + df[species[g]['TecIndex'], -1, yin, xright]*satiredg
                                   + df[species[g]['TecIndex'], -1, yout, xright]*satoredg)*vedge**2
                                  + (sum(df[species[g]['TecIndex'],-1,yin,xleft + 1 : xright]*satielem)
                                     + sum(df[species[g]['TecIndex'],-1,yout,xleft + 1 : xright]*satoelem)
                                     + sum(df[species[g]['TecIndex'], -1, yin + 1 : yout, xleft]*satlelem)
                                     + sum(df[species[g]['TecIndex'], -1, yin + 1 : yout, xright]*satrelem))*vedge*velem)
                                 + sum(sum(df[species[g]['TecIndex'],-1,yin + 1 : yout,xleft + 1 : xright]*satelem))*velem**2)*por
    return stock_domain

def oxiccells(limit,Trial,Het,Anis,gw,d,fpre,fsuf,yin,yout,xleft,xright,vars,gvarnames):
    oxiccells = np.zeros([len(Trial), 51, 1])
    for j in range(len(Trial)):
        df, massendtime, masstime, conctime, Velocity, head = conc_time(Trial[j],Het[j],Anis[j],
            gw,d,fpre,fsuf,yin,yout,xleft,xright,vars,gvarnames)
        c = []
        for k in range(51):
            c.append(np.count_nonzero(df[vars[gvarnames.index("DO")] - 3, np.shape(df)[1] - 1, k, :]>limit))
        oxiccells[j, :, 0] = c

    return oxiccells

def diversity_carbon (data, carbon_num, bio_num):

    """Function to evaluate the Shannon diversity and total
    carbon stock in the domain.
    
    Parameter
    ---------
    data : Array, float.
        Array of concentration of concentration of carbon and biomass species.
        Each column is a different species. Each row is a new time point.
    carbon_num : int.
        Number of carbon species. Default value is 1.
    bio_num : int.
        Number of microbial species. Default value is 1.
    
    """

    C = data[:,:carbon_num]
    B = data[:,carbon_num:]
    
    ## SHANNON DIVERSITY AND CARBON STOCK
    proportion = B/np.sum(B,axis=0)
    Shannon = -np.sum(proportion*np.log(proportion), axis = 1)
    total_C_stock = np.sum(C,axis=1) + np.sum(B, axis=1)
    C_stock = np.sum(C,axis=1)

    return Shannon, C_stock, total_C_stock

def normalize_carbon (data, carbon_num, bio_num, carbon_initial, biomass_initial):

    """Function to normalize carbon concentration by initials values.
    
    Parameter
    ---------
    data : Array, float.
        Array of concentration of concentration of carbon and biomass species.
        Each column is a different species. Each row is a new time point.
    carbon_num : int.
        Number of carbon species. Default value is 1.
    bio_num : int.
        Number of microbial species. Default value is 1.
    carbon_initial : Array, float.
        Array of initial concentrations of each carbon species.
    biomass_initial : Array, float.
        Array of initial concentrations of each biomass species.

    """


    C = data[:,:carbon_num]
    C_norm = np.sum(C, axis =1)/sum(carbon_initial)
    B = data[:,carbon_num:]
    
    return C_norm, B