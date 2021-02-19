# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:30:37 2020

@author: khurana
"""

from ogs5py import OGS, by_id
from gstools import SRF, TPLExponential
import numpy as np
import csv
import matplotlib.pyplot as plt

conditions= {
        "Bfo1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 200, "label" : "immobilebiomass"},
        "Bmo1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}, 
        "doc1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 800, "ic" : 2000, "label" : "chemical"},
        "dox1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 250, "ic" : 250, "label" : "chemical"},
        "Amm1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 60, "ic" : 30, "label" : "chemical"},
        "Bifo1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 5.0, "label" : "immobilebiomass"},
        "Bimo1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"},
        "Bfn1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 250, "label" : "immobilebiomass"},
        "Bmn1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}, 
        "nitra1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 250, "ic" : 250, "label" : "chemical"},
        "Bifn1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 5, "label" : "immobilebiomass"},
        "Bimn1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 5, "label" : "mobilebiomass"}, 
        "Bfs1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 200, "label" : "immobilebiomass"},
        "Bms1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}, 
        "sulpha1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 1500, "ic" : 1000, "label" : "chemical"},
        "Bifs1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 5, "label" : "immobilebiomass"},
        "Bims1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}, 
        "Bfa1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 200, "label" : "immobilebiomass"},
        "POM1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 5, "ic" : 5, "label" : "mobilebiomass"},
        "tr1":{"mindex" : 1, "diffindx" : 1, "diff" : "1e-9", "bc" : 20, "ic" : 0, "label" : "chemical"},
        "Bma1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}, 
        "Bifa1": {"mindex" : 0, "diffindx" : 0, "diff" : "1e-9", "bc" : 0, "ic" : 5, "label" : "immobilebiomass"},
        "Bima1": {"mindex" : 1, "diffindx" : 1, "diff" : 0, "bc" : 2, "ic" : 2, "label" : "mobilebiomass"}
        }

def genogsfiles2(foldername, modelname, variance, anisotropy, domainx, cellsizex, domainy, cellsizey, seedhet):
    
    model = OGS(task_root = foldername, task_id = modelname, output_dir = "out")
    model.gli.generate("rectangular", dim = 2, ori = (0.0,0.0,0.0), size = (domainx, domainy,0.0))
    model.gli.add_polyline("PLY_BC_BOT", [[0, 0, 0], [domainx, 0, 0]])
    model.gli.add_polyline("PLY_BC_TOP", [[0, domainy, 0], [domainx, domainy, 0]])
    model.gli.add_polyline("PLY_BC_LEFT", [[0, 0, 0], [0, domainy, 0]])
    model.gli.add_polyline("PLY_BC_RIGHT", [[domainx, 0, 0], [domainx, domainy, 0]])
    model.gli.add_points([0.0, domainy/2, 0.0], "SOURCE")
    model.msh.generate("rectangular", dim = 2, mesh_origin = (0.0, 0.0, 0.0), element_no = (int(domainx/cellsizex), int(domainy/cellsizey), 0), element_size = (cellsizex, cellsizey, 0))
    
    # covariance model for conductivity field
    cov_model = TPLExponential(dim=2, var = variance,  # variance (C calculated internally, so that `var` is 1)
                               len_low=0.01,  # lower truncation of the power law
                               len_scale= domainx/20,  # length scale (a.k.a. range), len_up = len_low + len_scale
                               anis=anisotropy,  # anisotropy between main direction and transversal ones
                               hurst=0.25)
    srf = SRF(model=cov_model, mean=-14, seed = seedhet)
    
    cond = np.exp(srf.mesh(model.msh).T)
    model.mpd.add(name=str(seedhet))
    model.mpd.add_block(  # edit recent mpd file
            MSH_TYPE="GROUNDWATER_FLOW",
            MMP_TYPE="PERMEABILITY",
            DIS_TYPE="ELEMENT",
            DATA=by_id(cond),
            )
    model.mmp.add_block(
            GEOMETRY_DIMENSION = 2,
            POROSITY = [1, 0.2],
            TORTUOSITY = [1, 1],
            STORAGE = [1, 0.0],
            PERMEABILITY_TENSOR = ["ISOTROPIC", 2e-06],
            PERMEABILITY_DISTRIBUTION = model.mpd.file_name,
            MASS_DISPERSION = [1, 0.02, 0.0],
            DENSITY = [1, 1500.0])
    
    model.pcs.add_block(  # set the process type
            PCS_TYPE="GROUNDWATER_FLOW",
            TIM_TYPE = "STEADY",
            NUM_TYPE="NEW", RELOAD = [1, 5000])
    for r in range(len(conditions)):
        model.pcs.add_block(
                PCS_TYPE = "MASS_TRANSPORT",
                TIM_TYPE = "STEADY",
                NUM_TYPE = "NEW", RELOAD = [1, writingfreq])
    model.num.add_block(  # numerical solver
               PCS_TYPE="GROUNDWATER_FLOW",
               LINEAR_SOLVER=[2, 6, 1.0e-14, 5000, 1.0, 100, 2],
               ELE_GAUSS_POINTS = 3
               )
    model.num.add_block(  # numerical solver
                PCS_TYPE="MASS_TRANSPORT",
    LINEAR_SOLVER=[2, 6, 1.0e-14, 5000, 1.0, 100, 2],
               ELE_GAUSS_POINTS = 3
    )
    model.st.add_block(  
            PCS_TYPE="GROUNDWATER_FLOW",
            PRIMARY_VARIABLE="HEAD",
            GEO_TYPE=["POLYLINE", "PLY_BC_LEFT"],
            DIS_TYPE=["CONSTANT_NEUMANN", 0.001],
            )
    for p in ["GROUNDWATER_FLOW", "MASS_TRANSPORT"]:
        model.tim.add_block(PCS_TYPE= p,
                            TIME_START = 0,
                            TIME_END= 365,
                            TIME_STEPS = [365000, 0.001],
                            TIME_UNIT = "DAY"
                            )
    model.msp.add_block(
            DENSITY = [1, 1500.0])
    for t,values in conditions.items():
        model.mcp.add_block(
                NAME = t,
                MOBILE = conditions[t]["mindex"],
                DIFFUSION = [conditions[t]["diffindx"], conditions[t]["diff"]])
    model.bc.add_block(PCS_TYPE="GROUNDWATER_FLOW",
               PRIMARY_VARIABLE="HEAD",
               GEO_TYPE=["POLYLINE", "PLY_BC_RIGHT"],
               DIS_TYPE=["CONSTANT", 0.0],
               )
    model.bc.add_block(PCS_TYPE="GROUNDWATER_FLOW",
               PRIMARY_VARIABLE="HEAD",
               GEO_TYPE=["POLYLINE", "PLY_BC_LEFT"],
               DIS_TYPE=["CONSTANT", 0.001],
               )
    for mobsp in list(t for t,values in conditions.items() if conditions[t]["label"]!="immobilebiomass"):
        model.bc.add_block(PCS_TYPE="MASS_TRANSPORT",
                           PRIMARY_VARIABLE = mobsp,
                           GEO_TYPE=["POINT", "SOURCE"],
                           DIS_TYPE=["CONSTANT", conditions[mobsp]["bc"]],
                           )
    model.ic.add_block(PCS_TYPE="GROUNDWATER_FLOW",
                       PRIMARY_VARIABLE="HEAD",
                       GEO_TYPE= "DOMAIN",
                       DIS_TYPE=["CONSTANT", 0.0],
                       )
    for t,values in conditions.items():
        model.ic.add_block(PCS_TYPE="MASS_TRANSPORT",
                           PRIMARY_VARIABLE = t,
               GEO_TYPE = "DOMAIN",
               DIS_TYPE=["CONSTANT", conditions[t]["ic"]],
               )
    model.mfp.add_block(
            FLUID_TYPE = "WATER",
            DESNITY = [1, 1000.0],
            VISCOSITY = [1, 0.001],
            HEAT_CAPACITY = [1, 0.0],
            HEAT_CONDUCTIVITY = [1, 0.0],
            PCS_TYPE = "HEAD")
    model.write_input()
    
    f = foldername + "/" + modelname + ".out"
    csvfile = open(f, "w")
    writer = csv.writer(
            csvfile,
            delimiter=" ",
            quotechar="\t",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
            )
    writer.writerow(["#OUTPUT"])
    writer.writerow([" $NOD_VALUES"])
    writer.writerow(["  HEAD"])
    writer.writerow(["  VELOCITY_X1"])
    writer.writerow(["  VELOCITY_Y1"])
    writer.writerow(["  POM1"])
    writer.writerow(["  tr1"])
    writer.writerow([" $GEO_TYPE"])
    writer.writerow(["  DOMAIN"])
    writer.writerow([" $TIM_TYPE"])
    writer.writerow(["  STEPS " + str(writingfreq)])
    writer.writerow([" $DAT_TYPE"])
    writer.writerow(["  TECPLOT"])
    writer.writerow([" $PCS_TYPE"])
    writer.writerow(["  ALL"])
    writer.writerow(["#STOP"])
    csvfile.close()
    
    plt.figure()
    plt.imshow(cond.reshape(int(domainy/cellsizey),int(domainx/cellsizex)), origin="lower")
    plt.title ("Heterogeneity field: " + foldername)
    plt.show()
    plt.savefig(foldername + "/hetfield.png", dpi = 300)
    
    return plt

voptions = [0.1, 1, 10]
ls = [0.1]
anoptions = [2, 5, 10]
x = 0.01
y = 0.01
counter = 37
for va in voptions:
    for an in anoptions:
        foldername = str(counter)
        modelname = "2Dmodel"
        writingfreq = 1
        genogsfiles2(str(foldername), str(modelname), variance = va, anisotropy = an, domainx = 1.0, cellsizex = x, domainy = 0.5, cellsizey = y, seedhet = counter)
        counter += 1
        
def genogsfiles3(foldername, modelname, variance, anisotropy, domainx, cellsizex, domainy, cellsizey, domainz, cellsizez, seedhet):
    
    model = OGS(task_root = foldername, task_id = modelname, output_dir = "out")
    model.gli.generate("rectangular", dim = 3, ori = (0.0,0.0,0.0), size = (domainx, domainy, domainz))
    #model.gli.add_polyline("BOT_F", [[0, 0, 0], [domainx, 0, 0]])
    #model.gli.add_surface("RIGHT_S", ["RIGHT_F", "RIGHT_TFB", "RIGHT_B", "RIGHT_BFB"])
    model.msh.generate("rectangular", dim = 3, mesh_origin = (0.0, 0.0, 0.0), element_no = (int(domainx/cellsizex), int(domainy/cellsizey), int(domainz/cellsizez)), element_size = (cellsizex, cellsizey, cellsizez))
    
    # covariance model for conductivity field
    cov_model = TPLExponential(dim=3, var = variance,  # variance (C calculated internally, so that `var` is 1)
                               len_low=0.01,  # lower truncation of the power law
                               len_scale=1.0,  # length scale (a.k.a. range), len_up = len_low + len_scale
                               #anis=anisotropy,  # anisotropy between main direction and transversal ones
                               hurst=0.25)
    srf = SRF(model=cov_model, mean=-14, seed = seedhet)
    
    cond = np.exp(srf.mesh(model.msh))
    model.mpd.add(name=str(seedhet))
    model.mpd.add_block(  # edit recent mpd file
            MSH_TYPE="GROUNDWATER_FLOW",
            MMP_TYPE="PERMEABILITY",
            DIS_TYPE="ELEMENT",
            DATA=by_id(cond),
            )
    model.mmp.add_block(
            GEOMETRY_DIMENSION = 3,
            POROSITY = [1, 0.2],
            TORTUOSITY = [1, 1],
            STORAGE = [1, 0.0],
            PERMEABILITY_TENSOR = ["ISOTROPIC", 2],
            PERMEABILITY_DISTRIBUTION = model.mpd.file_name,
            MASS_DISPERSION = [1, 0.02, 0.0],
            DENSITY = [1, 1500.0])
    model.pcs.add_block(  # set the process type
            PCS_TYPE="GROUNDWATER_FLOW",
            TIM_TYPE = "STEADY",
            NUM_TYPE="NEW", RELOAD = [1, 5000])
    for r in range(len(conditions)):
        model.pcs.add_block(
                PCS_TYPE = "MASS_TRANSPORT",
                TIM_TYPE = "STEADY",
                NUM_TYPE = "NEW", RELOAD = [1, writingfreq])
    model.num.add_block(  # numerical solver
               PCS_TYPE="GROUNDWATER_FLOW",
               LINEAR_SOLVER=[2, 6, 1.0e-14, 5000, 1.0, 100, 2],
               ELE_GAUSS_POINTS = 3
               )
    model.num.add_block(  # numerical solver
                PCS_TYPE="MASS_TRANSPORT",
    LINEAR_SOLVER=[2, 6, 1.0e-14, 5000, 1.0, 100, 2],
               ELE_GAUSS_POINTS = 3
    )
    model.st.add_block(  
            PCS_TYPE="GROUNDWATER_FLOW",
            PRIMARY_VARIABLE="HEAD",
            GEO_TYPE=["SURFACE", "boundary_3"],
            DIS_TYPE=["CONSTANT_NEUMANN", 0.001],
            )
    for p in ["GROUNDWATER_FLOW", "MASS_TRANSPORT"]:
        model.tim.add_block(PCS_TYPE= p,
                            TIME_START = 0,
                            TIME_END= 365,
                            TIME_STEPS = [365000, 0.001],
                            TIME_UNIT = "DAY"
                            )
    model.msp.add_block(
            DENSITY = [1, 1500.0])
    for t,values in conditions.items():
        model.mcp.add_block(
                NAME = t,
                MOBILE = conditions[t]["mindex"],
                DIFFUSION = [conditions[t]["diffindx"], conditions[t]["diff"]])
    model.bc.add_block(PCS_TYPE="GROUNDWATER_FLOW",
               PRIMARY_VARIABLE="HEAD",
               GEO_TYPE=["SURFACE", "boundary_1"],
               DIS_TYPE=["CONSTANT", 0.0],
               )
    for mobsp in list(t for t,values in conditions.items() if conditions[t]["label"]!="immobilebiomass"):
        model.bc.add_block(PCS_TYPE="MASS_TRANSPORT",
                           PRIMARY_VARIABLE = mobsp,
                           GEO_TYPE=["SURFACE", "boundary_3"],
                           DIS_TYPE=["CONSTANT", conditions[mobsp]["bc"]],
                           )
    model.ic.add_block(PCS_TYPE="GROUNDWATER_FLOW",
                       PRIMARY_VARIABLE="HEAD",
                       GEO_TYPE= "DOMAIN",
                       DIS_TYPE=["CONSTANT", 0.0],
                       )
    for t,values in conditions.items():
        model.ic.add_block(PCS_TYPE="MASS_TRANSPORT",
                           PRIMARY_VARIABLE = t,
               GEO_TYPE = "DOMAIN",
               DIS_TYPE=["CONSTANT", conditions[t]["ic"]],
               )
    model.mfp.add_block(
            FLUID_TYPE = "WATER",
            DENSITY = [1, 1000.0],
            VISCOSITY = [1, 0.001],
            HEAT_CAPACITY = [1, 0.0],
            HEAT_CONDUCTIVITY = [1, 0.0],
            PCS_TYPE = "HEAD")
    model.write_input()
    
    f = foldername + "/" + modelname + ".out"
    csvfile = open(f, "w")
    writer = csv.writer(
            csvfile,
            delimiter=" ",
            quotechar="\t",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
            )
    writer.writerow(["#OUTPUT"])
    writer.writerow([" $NOD_VALUES"])
    writer.writerow(["  HEAD"])
    writer.writerow(["  VELOCITY_X1"])
    writer.writerow(["  VELOCITY_Y1"])
    writer.writerow(["  VELOCITY_Z1"])
    writer.writerow(["  POM1"])
    writer.writerow(["  tr1"])
    writer.writerow([" $GEO_TYPE"])
    writer.writerow(["  DOMAIN"])
    writer.writerow([" $TIM_TYPE"])
    writer.writerow(["  STEPS " + str(writingfreq)])
    writer.writerow([" $DAT_TYPE"])
    writer.writerow(["  TECPLOT"])
    writer.writerow([" $PCS_TYPE"])
    writer.writerow(["  ALL"])
    writer.writerow(["#STOP"])
    csvfile.close()
    
    plt.figure()
    plt.imshow(cond.reshape(int(domainz/cellsizez)*int(domainy/cellsizey),int(domainx/cellsizex)), origin="lower")
    plt.title ("Heterogeneity field: " + foldername)
    plt.show()
    plt.savefig(foldername + "/hetfield.png", dpi = 300)
    
    return plt

voptions = [0.1, 1, 10]
ls = [0.1]
anoptions = [2, 5, 10]
x = 0.01
y = 0.01
z = 0.01
counter = 37
for va in voptions:
    for an in anoptions:
        foldername = str(counter)
        modelname = "2dmodel"
        writingfreq = 2
        genogsfiles2(str(foldername), str(modelname), variance = va, anisotropy = an, domainx = 4.0, cellsizex = x, domainy = 1.0, cellsizey = 0.01, seedhet = counter)#, domainz = 2.0, cellsizez = z)
        counter += 1