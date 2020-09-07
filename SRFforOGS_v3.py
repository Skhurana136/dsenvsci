# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:59:49 2019

@author: khurana
"""

from gstools import SRF, TPLExponential
import matplotlib.pyplot as plt
import numpy as np
import csv

y = np.linspace(0, 0.5, 51)
x = np.linspace(0, 0.3, 31)
extnc = ".dat"
extnf = "png"
pref = "A"
randomkeyword = "-T"

v = [5, 5, 5]
# v2 = [2, 5, 10]
ls = [0.1]
an = [2, 5, 10]
counter = 37  # First realization series
counter = 46  # Second realization series
counter = 55  # Third realization series
counter = (
    64
)  # Fourth realization with an extra variance of 5 which was not included earlier
counter = 76  # three more realizations for variance of 5
for i in range(len(v)):
    for j in range(len(ls)):
        for k in range(len(an)):
            model = TPLExponential(
                dim=2,  # spatial dimension
                var=v[i],  # variance (C calculated internally, so that `var` is 1)
                len_low=0.01,  # lower truncation of the power law
                len_scale=ls[
                    j
                ],  # length scale (a.k.a. range), len_up = len_low + len_scale
                anis=an[k],  # anisotropy between main direction and transversal ones
                hurst=0.25,  # hurst coefficient from the power law
            )
            srf = SRF(model, mean=-14, mean_in=-14)
            field = srf((x, y), mesh_type="structured", seed=counter)
            #            fig = plt.figure()
            plt.imshow(field.T, origin="lower")
            plt.show()
            #            fig.savefig("A1.png")
            #            plt.close()
            perm = 2 * np.exp(field)
            print(
                pref + str(counter),
                v[i],
                ls[j],
                an[k],
                np.min(perm),
                np.mean(perm),
                np.max(perm),
                np.std(perm),
            )

            fname = pref + str(counter) + randomkeyword + extnc
            permfile = np.arange(perm.size).reshape(
                (x.size, y.size)
            )  # Contains f(x, y)
            assert permfile.shape == (x.size, y.size)
            csvfile = open(fname, "w")
            writer = csv.writer(
                csvfile,
                delimiter=" ",
                quotechar="\t",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            writer.writerow(["#MEDIUM_PROPERTIES_DISTRIBUTED"])
            writer.writerow(["$MSH_TYPE"])
            writer.writerow([" GROUNDWATER_FLOW"])
            writer.writerow(["$MMP_TYPE"])
            writer.writerow([" PERMEABILITY"])
            writer.writerow(["$DIS_TYPE"])
            writer.writerow([" NEAREST_VALUE"])
            writer.writerow([";hj"])
            writer.writerow(["$CONVERSION_FACTOR"])
            writer.writerow([" 1.0"])
            writer.writerow(["$DATA"])
            # writer.writerow([";x, y, Z, perm"])
            for x_index in range(x.size):
                for y_index in range(y.size):
                    writer.writerow([x[x_index], y[y_index], 0, perm[x_index, y_index]])
            writer.writerow(["#STOP"])
            csvfile.close()
            counter = counter + 1

#For generating a bigger domain
y = np.linspace(0, 0.5*2, 51*2-1)
x = np.linspace(0, 0.3, 31)
extnc = ".dat"
extnf = "png"
pref = "Double"
randomkeyword = "-T"

#v = [0.1, 1, 10]
ls = [0.1]
an = [2, 5, 10]
counter = 37  # First realization series
numberrealizations = 3
v = [0.1, 1, 5, 10]
numberrealizations = 1
counter = 64 #Second realization series
for number in range(numberrealizations):
    for i in range(len(v)):
        for j in range(len(ls)):
            for k in range(len(an)):
                model = TPLExponential(
                        dim=2,  # spatial dimension
                        var=v[i],  # variance (C calculated internally, so that `var` is 1)
                        len_low=0.01,  # lower truncation of the power law
                        len_scale=ls[
                                j
                                ],  # length scale (a.k.a. range), len_up = len_low + len_scale
                                anis=an[k],  # anisotropy between main direction and transversal ones
                                hurst=0.25,  # hurst coefficient from the power law
                                )
                srf = SRF(model, mean=-14, mean_in=-14)
                field = srf((x, y), mesh_type="structured", seed=counter)
                fig = plt.figure()
                plt.imshow(field.T, origin="lower")
                plt.show()
                fig.savefig(pref + str(counter) + ".png")
                plt.close()
                perm = 2 * np.exp(field)
                print(
                        pref + str(counter),
                        v[i],
                        ls[j],
                        an[k],
                        np.min(perm),
                        np.mean(perm),
                        np.max(perm),
                        np.std(perm),
                        )
                fname = pref + str(counter) + randomkeyword + extnc
                permfile = np.arange(perm.size).reshape(
                        (x.size, y.size)
                        )  # Contains f(x, y)
                assert permfile.shape == (x.size, y.size)
                csvfile = open(fname, "w")
                writer = csv.writer(
                        csvfile,
                        delimiter=" ",
                        quotechar="\t",
                        quoting=csv.QUOTE_MINIMAL,
                        lineterminator="\n",
                        )
                writer.writerow(["#MEDIUM_PROPERTIES_DISTRIBUTED"])
                writer.writerow(["$MSH_TYPE"])
                writer.writerow([" GROUNDWATER_FLOW"])
                writer.writerow(["$MMP_TYPE"])
                writer.writerow([" PERMEABILITY"])
                writer.writerow(["$DIS_TYPE"])
                writer.writerow([" NEAREST_VALUE"])
                writer.writerow([";hj"])
                writer.writerow(["$CONVERSION_FACTOR"])
                writer.writerow([" 1.0"])
                writer.writerow(["$DATA"])
                # writer.writerow([";x, y, Z, perm"])
                for x_index in range(x.size):
                    for y_index in range(y.size):
                        writer.writerow([x[x_index], y[y_index], 0, perm[x_index, y_index]])
                writer.writerow(["#STOP"])
                csvfile.close()
                counter = counter + 1
                
#For generating a smaller domain
y = np.linspace(0, 0.5/2, 26)
x = np.linspace(0, 0.3, 31)
extnc = ".dat"
extnf = "png"
pref = "Half"
randomkeyword = "-T"

v = [0.1, 1, 5, 10]
ls = [0.1]
an = [2, 5, 10]
counter = 37  # First realization series
numberrealizations = 3
for number in range(numberrealizations):
    for i in range(len(v)):
        for j in range(len(ls)):
            for k in range(len(an)):
                model = TPLExponential(
                        dim=2,  # spatial dimension
                        var=v[i],  # variance (C calculated internally, so that `var` is 1)
                        len_low=0.01,  # lower truncation of the power law
                        len_scale=ls[
                                j
                                ],  # length scale (a.k.a. range), len_up = len_low + len_scale
                                anis=an[k],  # anisotropy between main direction and transversal ones
                                hurst=0.25,  # hurst coefficient from the power law
                                )
                srf = SRF(model, mean=-14, mean_in=-14)
                field = srf((x, y), mesh_type="structured", seed=counter)
                fig = plt.figure()
                plt.imshow(field.T, origin="lower")
                plt.show()
                fig.savefig(pref + str(counter) + ".png")
                plt.close()
                perm = 2 * np.exp(field)
                print(
                        pref + str(counter),
                        v[i],
                        ls[j],
                        an[k],
                        np.min(perm),
                        np.mean(perm),
                        np.max(perm),
                        np.std(perm),
                        )
                fname = pref + str(counter) + randomkeyword + extnc
                permfile = np.arange(perm.size).reshape(
                        (x.size, y.size)
                        )  # Contains f(x, y)
                assert permfile.shape == (x.size, y.size)
                csvfile = open(fname, "w")
                writer = csv.writer(
                        csvfile,
                        delimiter=" ",
                        quotechar="\t",
                        quoting=csv.QUOTE_MINIMAL,
                        lineterminator="\n",
                        )
                writer.writerow(["#MEDIUM_PROPERTIES_DISTRIBUTED"])
                writer.writerow(["$MSH_TYPE"])
                writer.writerow([" GROUNDWATER_FLOW"])
                writer.writerow(["$MMP_TYPE"])
                writer.writerow([" PERMEABILITY"])
                writer.writerow(["$DIS_TYPE"])
                writer.writerow([" NEAREST_VALUE"])
                writer.writerow([";hj"])
                writer.writerow(["$CONVERSION_FACTOR"])
                writer.writerow([" 1.0"])
                writer.writerow(["$DATA"])
                # writer.writerow([";x, y, Z, perm"])
                for x_index in range(x.size):
                    for y_index in range(y.size):
                        writer.writerow([x[x_index], y[y_index], 0, perm[x_index, y_index]])
                writer.writerow(["#STOP"])
                csvfile.close()
                counter = counter + 1
                
#Domain length - 1.25 meters
y = np.linspace(0, 1.25, 126)
x = np.linspace(0, 0.3, 31)
extnc = ".dat"
extnf = "png"
pref = "Doublehalf"
randomkeyword = "-T"

v = [0.1, 1, 10]
ls = [0.1]
an = [2, 5, 10]
counter = 37  # First realization series
numberrealizations = 3
for number in range(numberrealizations):
    for i in range(len(v)):
        for j in range(len(ls)):
            for k in range(len(an)):
                model = TPLExponential(
                        dim=2,  # spatial dimension
                        var=v[i],  # variance (C calculated internally, so that `var` is 1)
                        len_low=0.01,  # lower truncation of the power law
                        len_scale=ls[
                                j
                                ],  # length scale (a.k.a. range), len_up = len_low + len_scale
                                anis=an[k],  # anisotropy between main direction and transversal ones
                                hurst=0.25,  # hurst coefficient from the power law
                                )
                srf = SRF(model, mean=-14, mean_in=-14)
                field = srf((x, y), mesh_type="structured", seed=counter)
                fig = plt.figure()
                plt.imshow(field.T, origin="lower")
                plt.show()
                fig.savefig(pref + str(counter) + ".png")
                plt.close()
                perm = 2 * np.exp(field)
                print(
                        pref + str(counter),
                        v[i],
                        ls[j],
                        an[k],
                        np.min(perm),
                        np.mean(perm),
                        np.max(perm),
                        np.std(perm),
                        )
                fname = pref + str(counter) + randomkeyword + extnc
                permfile = np.arange(perm.size).reshape(
                        (x.size, y.size)
                        )  # Contains f(x, y)
                assert permfile.shape == (x.size, y.size)
                csvfile = open(fname, "w")
                writer = csv.writer(
                        csvfile,
                        delimiter=" ",
                        quotechar="\t",
                        quoting=csv.QUOTE_MINIMAL,
                        lineterminator="\n",
                        )
                writer.writerow(["#MEDIUM_PROPERTIES_DISTRIBUTED"])
                writer.writerow(["$MSH_TYPE"])
                writer.writerow([" GROUNDWATER_FLOW"])
                writer.writerow(["$MMP_TYPE"])
                writer.writerow([" PERMEABILITY"])
                writer.writerow(["$DIS_TYPE"])
                writer.writerow([" NEAREST_VALUE"])
                writer.writerow([";hj"])
                writer.writerow(["$CONVERSION_FACTOR"])
                writer.writerow([" 1.0"])
                writer.writerow(["$DATA"])
                # writer.writerow([";x, y, Z, perm"])
                for x_index in range(x.size):
                    for y_index in range(y.size):
                        writer.writerow([x[x_index], y[y_index], 0, perm[x_index, y_index]])
                writer.writerow(["#STOP"])
                csvfile.close()
                counter = counter + 1