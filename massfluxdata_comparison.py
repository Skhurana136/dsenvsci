# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""

import pandas as pd

#Load data
path_data = "//msg-filer2/scratch_60_days/khurana/massflux_Original_complete.csv"
data = pd.read_csv(path_data, sep = "\t")
data.columns
data.dtypes

regimes = data.Regime.unique().tolist()
Time_series = data.Time_series.unique().tolist()
chem_series = data.Chem.unique().tolist()
trial_series = data.Trial.unique().tolist()
spatial_base = 'H'
temp_base = 0

for r in regimes:
    for t in trial_series:
        for c in chem_series:
            rel_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t) & (data.Time_series == 0)]['reldelmassflux'].values[0]
            print(rel_base)
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'relbase'] = rel_base
        
data['reldelmassflux_fraction'] = data['delmassflux']/data['relbase']

subdata = data[data['Trial']!=52]

x = subdata.groupby(["Regime", "Chem", "Variance", "Anisotropy"]).mean()