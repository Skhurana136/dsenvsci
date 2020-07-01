# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 11:17:20 2020

@author: khurana
"""
import seaborn as sns
import matplotlib.pyplot as plt
from data_reader.data_processing import tracerstudies

def plot_tracer(filename):
    combined_tracer = tracerstudies(filename)
    combined_tracer["%fraction"] = combined_tracer["fraction"] * 100
    combined_tracer["%fraction_withslow"] = combined_tracer["fraction_withslow"] * 100
    sns.set(rc={"figure.figsize": (7, 4)})
    sns.set_style("whitegrid")
    sns.boxplot(
        x="Xlabels",
        y="%fraction",
        hue="Regime",
        data=combined_tracer,
        hue_order=["Slow", "Medium", "Fast"],
        palette=["coral", "mediumseagreen", "steelblue"],
    )
    plt.xlabel("Variance:Anisotropy")
    plt.ylabel("% of homogeneous scenario")
    plt.title("Time taken for tracer breakthrough")

    return plt
