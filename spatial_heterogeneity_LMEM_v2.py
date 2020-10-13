# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

#temporal_variation_factor_analysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import statsmodels.formula.api as smf

head_path = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/headatinlet.csv"
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

def corrfunc (data):
    normdata = data - np.mean(data)
    autocorr_f = (np.correlate(normdata, normdata, mode='full'))/(np.std(data)**2 * np.shape(data)[0])
    return autocorr_f

cor1 = corrfunc(head["H1"])[5474:]
cor2 = corrfunc(head["H2"])[5474:]
cor3 = corrfunc(head["H3"])[5474:]

cor1tim = np.where(cor1 < 0.75)[0][0]
cor2tim = np.where(cor2 < 0.75)[0][0]
cor3tim = np.where(cor3 < 0.75)[0][0]
mincor = cor2.min()

cleanup_tim = {"cov" : {1 : cov1, 2 : cov2, 5 : cov3},
               "timtrace" : {1 : cor1tim, 2: cor2tim, 5 : cor3tim},
               "chem_factors" : {'DOC':0,'DO':1, 'Ammonium':2, 'Nitrate':3, 'Sulphate':4, 'Particulate organic matter':5,
                                 'Mobile active aerobic degraders':6, 'Mobile active nitrate reducers':7,
                                 'Mobile active sulphate reducers':8, 'Mobile active ammonia oxidizers':9,
                                 'Mobile inactive aerobic degraders':10, 'Mobile inactive nitrate reducers':11,
                                 'Mobile inactive sulphate reducers':12, 'Mobile inactive ammonia oxidizers':13,
                                 'Nitrogen':14, 'TOC':15},
                   "Pe" : {'Slow':2,'Medium':11, 'Fast':22}}

directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "Normalized_RMSamplitude_chem.csv"
data = pd.read_csv(directory + filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")
chemstoplot = ["DOC", "DO", "TOC", "Nitrogen"]
data["cov"] = data ["Time_series"]

data  = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")

data["timtrace"] = data ["Time_series"]
data["chem_factors"] = data ["Chem"]
data["Pe"] = data ["Regime"]
data.replace(cleanup_tim, inplace = True)

dadata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Conc_da_ss.csv", sep = "\t")

mdata = pd.merge(data, dadata[["Regime", "Trial", "Chem", "%reldelmassflux", "%fraction_rel_delmf"]], on = ["Regime", "Trial", "Chem"])
mdata.loc[mdata["PeDa"] > 40, "PeDamark"] = 3
mdata.loc[(mdata["PeDa"] > 20) & (mdata["PeDa"] < 40), "PeDamark"] = 2
mdata.loc[(mdata["PeDa"] > 1) & (mdata["PeDa"] < 20), "PeDamark"] = 1
mdata.loc[mdata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : "Da/Pe < 1",
          1 : "1 < Da/Pe < 15",
          2 : "15 < Da/Pe < 40",
          3 : "Da/Pe > 40"}

#Steady state
allfeatures = dadata.columns.tolist()
features = ['Variance', 'Anisotropy', 'Regime', 'Chem', 'Breakthroughtime', 'fraction', 'avgconc_in', 'avgconc_out', 'PeDa']
yfeature = ['%reldelmassflux', 'del4massflux']

mydata = dadata[features+yfeature]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata["Variance"]= mydata["Variance"].astype(str)
mydata["Anisotropy"]= mydata["Anisotropy"].astype(str)

colsymdict = {"Slow" : "indianred", "Medium" : "g", "Fast" : "steelblue",
              "DOC": 'D', "DO" : '^', "Nitrogen" : 's', "TOC" : 'o'}

#Plot yfeature with breakthrough time divided by regime, variance and anisotropy
import seaborn as sns
sns.pairplot(mydata, kind="scatter", hue="Regime", palette=colsymdict)

#Variance and Anisotropy and avg_conc_in don't vary with scenarios, or flow regimes
nrows = 4
ncols = 3
vartoplot = mydata.Variance.unique().tolist()
antoplot = mydata.Anisotropy.unique().tolist()
regtoplot = mydata.Regime.unique().tolist()

for xplot in ["Breakthroughtime", "fraction"]:
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, sharex = True, sharey = True)
    vr = 0
    for v in vartoplot[1:]:
        vdata = mydata[mydata.Variance.eq(v)]
        ar = 0
        for a in antoplot[1:]:
            adata = vdata[vdata.Anisotropy.eq(a)]
            idx = axes[vr][ar]
            for c in chemstoplot:
                for r in regtoplot:
                    cdata = adata[(adata.Chem.eq(c)) & (adata.Regime.eq(r))]
                    idx.scatter("fraction", "perreldelmassflux", alpha = 0.5, c = colsymdict[r],
                                marker = colsymdict[c], data = cdata)
            ar += 1
    vr += 1
    idx.set_xlabel("Breakthrough time (days)")
    idx.set_ylabel("Relative change in mass flux (%)")

#Breakthrough time varies with variance and anisotropy, and
#relative mass flux does not vary with breakthrough time
#Both variance and anisotropy result in a different breakthrough time
#Both variance and anisotropy result in a different impact on breakthrough time

#So fixed effect - Breakthrough time or fraction
#Random effect - Regime, Variance, Anisotropy
#1. Regime#2. Variance#3. Anisotropy
xplotlist = ["Breakthroughtime", "fraction"]
yplot = 'del4massflux'
y = mydata[yplot]

row = []
for xplot in xplotlist:
    print(xplot)
    f1 = yplot + ' ~ ' + xplot
    for g in ["Regime", "Anisotropy"]:
        print(g)
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g])
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        row.append([yplot, xplot,g,RMSE, "Fixed"])
    
#Using variance as a group variable results in a singular matrix, so I am not pursuing it.
#Anisotropy as a group variable is less powerful than regime, so I will not use it either.
#Can be used in combination with variance and anisotropy
#Interested in pursuing random effects of regime on breakthrough time
for xplot in xplotlist:
    print(xplot)
    f1 = yplot + ' ~ ' + xplot
    for g in ["Regime", "Anisotropy"]:
        print(g)
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~'+xplot)
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        row.append([yplot, xplot,g,RMSE, "RandomFixed"])

#Random fixed improved by a minuscule amount.
#Nested random effects will require data transformation and addition of new columns
mydata["R_V"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str)
mydata["R_A"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str)
mydata["V_A"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
mydata["R_V_A"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
for xplot in xplotlist:
    print(xplot)
    f1 = yplot + ' ~ ' + xplot
    for g in ["R_V", "R_A", "V_A", "R_V_A"]:
        print(g)
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~'+xplot)
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        row.append([yplot, xplot,g,RMSE, "RandomNestedFixed"])

results = pd.DataFrame.from_records(row, columns = ["response", "independent", "Details", "RMSE", "Model_Type"])
results.to_csv("X:/RMSE_results_withbasecase.csv", sep = ',')

mydata["R_V_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Chem"]
mydata["R_A_C"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["V_A_C"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["R_V_A_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata ["R_C"] = mydata["Regime"] + "_" + mydata["Chem"]

groups = ["Regime", "Anisotropy", "R_V", "R_A", "V_A", "R_V_A", "Chem", "R_V_C", "R_A_C", "V_A_C", "R_V_A_C", "R_C"]
summarycols = ["Group", 0, 1, 2, 3, 'Coef.', 'Std.Err.', 'z', 'P>|z|', '[0.025', '0.975]']
allsummaries = pd.DataFrame(columns = summarycols)
oneshot = []
for xplot in xplotlist:
    f1 = yplot + ' ~ ' + xplot
    for g in groups:
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~'+xplot)
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        oneshot.append([yplot, xplot,g,RMSE, "RandomFixed"])
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0) 

results = pd.DataFrame.from_records(oneshot, columns = ["response", "independent", "Details", "RMSE", "Model_Type"])
results.to_csv("X:/RMSE_results_morey.csv", sep = ',')
allsummaries.to_csv("X:/Grouped_summary_morey.csv", sep = ',')

#Introducing Chemical made a big impact on performance on the model.
#Chemical should be a dependent variable in conjunction with breakthrough time
xplotlist = ["Breakthroughtime", "fraction", "Breakthroughtime + Chem", "fraction + Chem"]
allsummaries = pd.DataFrame(columns = summarycols)
oneshot = []
for rex in ["Breakthroughtime", "fraction"]:
    for yplot in ['perreldelmassflux', 'del4massflux']:
        y = mydata[yplot]
        for xplot in xplotlist:
            f1 = yplot + ' ~ ' + xplot
            for g in groups[7:]:
                mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~'+ rex)
                mdf = mymd.fit()
                y_predict = mdf.fittedvalues
                RMSE = np.sqrt(((y-y_predict)**2).mean())
                oneshot.append([yplot, xplot,g,RMSE, "RandomFixed", rex])
                summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
                summary["Group"] = g
                allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0) 

results = pd.DataFrame.from_records(oneshot, columns = ["response", "independent", "Details", "RMSE", "Model_Type", "Recurse"])
results.to_csv("X:/RMSE_results_wthChem.csv", sep = ',')
allsummaries.to_csv("X:/Grouped_summary_withChem.csv", sep = ',')

#Trying the same but grouping with Da and Pe instead
mydata = dadata[["%reldelmassflux", "PeDa", "Breakthroughtime", "fraction", "del4massflux", "conc_Da", "Pe", "Chem"]]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
mydata.loc[mydata["PeDa"] > 40, "PeDamark"] = 3
mydata.loc[(mydata["PeDa"] > 20) & (mydata["PeDa"] < 40), "PeDamark"] = 2
mydata.loc[(mydata["PeDa"] > 1) & (mydata["PeDa"] < 20), "PeDamark"] = 1
mydata.loc[mydata["PeDa"] < 1, "PeDamark"] = 0
groups = ["conc_Da","Pe", "PeDamark"]
summarycols = ["Group", 0, 1, 2, 3, 'Coef.', 'Std.Err.', 'z', 'P>|z|', '[0.025', '0.975]']
allsummaries = pd.DataFrame(columns = summarycols)
oneshot = []
for yplot in ['perreldelmassflux', 'del4massflux']:
    y = mydata[yplot]
    for xplot in ['Breakthroughtime', 'fraction']:
        f1 = yplot + ' ~ ' + xplot
        for g in groups:
            mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g] , re_formula = '~'+xplot)
            mdf = mymd.fit()
            y_predict = mdf.fittedvalues
            RMSE = np.sqrt(((y-y_predict)**2).mean())
            oneshot.append([yplot, xplot,g,RMSE, "RandomFixed", xplot])
            summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
            summary["Group"] = g
            allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
            plt.figure()
            plt.scatter(y, y_predict, label = xplot + yplot)
            plt.ylabel (yplot)
            plt.xlabel (xplot)
    for xplot in ['Breakthroughtime', 'fraction', 'conc_Da', 'Pe']:
        f1 = yplot + ' ~ ' + xplot
        for g in groups[2:]:
            mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g] , re_formula = '~'+xplot)
            mdf = mymd.fit()
            y_predict = mdf.fittedvalues
            RMSE = np.sqrt(((y-y_predict)**2).mean())
            oneshot.append([yplot, xplot,g,RMSE, "RandomFixed", xplot])
            summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
            summary["Group"] = g
            allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
            plt.figure()
            plt.scatter(y, y_predict, label = xplot + yplot)
            plt.ylabel (yplot)
            plt.xlabel (xplot)
    for xplot in ['Breakthroughtime', 'fraction', 'Breakthroughtime + conc_Da', 'Breakthroughtime + Pe', 'fraction + Pe', 'fraction + conc_Da']:
        f1 = yplot + ' ~ ' + xplot
        for g in groups[2:]:
            mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g] , re_formula = '~'+xplot)
            mdf = mymd.fit()
            y_predict = mdf.fittedvalues
            RMSE = np.sqrt(((y-y_predict)**2).mean())
            oneshot.append([yplot, xplot,g,RMSE, "RandomFixed", xplot])
            summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
            summary["Group"] = g
            allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
            plt.figure()
            plt.scatter(y, y_predict, label = xplot + yplot)
            plt.ylabel (yplot)
            plt.xlabel (xplot)
    for xplot in ['Breakthroughtime', 'fraction', 'Breakthroughtime + conc_Da', 'Breakthroughtime + Pe', 'fraction + Pe', 'fraction + conc_Da']:
        f1 = yplot + ' ~ ' + xplot
        for g in groups:
            mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g])
            mdf = mymd.fit()
            y_predict = mdf.fittedvalues
            RMSE = np.sqrt(((y-y_predict)**2).mean())
            oneshot.append([yplot, xplot,g,RMSE, "RandomFixed", xplot])
            summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
            summary["Group"] = g
            allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
            plt.figure()
            plt.scatter(y, y_predict, label = xplot + yplot)
            plt.ylabel (yplot)
            plt.xlabel (xplot)
     
results = pd.DataFrame.from_records(oneshot, columns = ["response", "independent", "Details", "RMSE", "Model_Type", "Recurse"])
results.to_csv("X:/RMSE_results_wthChem_PeDa.csv", sep = ',')
allsummaries.to_csv("X:/Grouped_summary_withChem_PeDa.csv", sep = ',')
        

###Trying the same with random forests###
from sklearn.ensemble import RandomForestRegressor as rf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import metrics

allfeatures = dadata.columns.tolist()
features = ['Variance', 'Anisotropy', 'Pe', 'conc_Da', 'Chem', 'Breakthroughtime', 'fraction', 'PeDa']
yfeature = ['%reldelmassflux', 'del4massflux']

mydata = dadata[features+yfeature]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
mydata.loc[mydata["PeDa"] > 40, "PeDamark"] = 3
mydata.loc[(mydata["PeDa"] > 20) & (mydata["PeDa"] < 40), "PeDamark"] = 2
mydata.loc[(mydata["PeDa"] > 1) & (mydata["PeDa"] < 20), "PeDamark"] = 1
mydata.loc[mydata["PeDa"] < 1, "PeDamark"] = 0

X = mydata[["fraction", "Pe", "conc_Da", "Breakthroughtime"]]
y = mydata[["del4massflux"]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))

X = mydata[["Variance", "Anisotropy", "fraction", "Pe", "conc_Da", "Breakthroughtime"]]
y = mydata[["del4massflux"]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))

X = mydata[["fraction", "PeDa", "Breakthroughtime"]]
y = mydata[["del4massflux"]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))


#Final outputs for paper:
mydata = dadata[["Breakthroughtime", "fraction", "del4massflux", "Regime", "Chem", "Variance", "Anisotropy"]]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata["R_V"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str)
mydata["R_A"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str)
mydata["V_A"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
mydata["R_V_A"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
mydata["R_V_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Chem"]
mydata["R_A_C"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["V_A_C"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["R_V_A_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata ["R_C"] = mydata["Regime"] + "_" + mydata["Chem"]
yplot = 'del4massflux'
y = mydata[yplot]
xplotlist = ['fraction', 'fraction + Chem']
groups = ["Regime", "Variance", "R_V", "V_A", "R_V_A", "Chem", "R_V_C", "V_A_C", "R_V_A_C", "R_C"]
summarycols = ["Group", 0, 1, 2, 3, 'Coef.', 'Std.Err.', 'z', 'P>|z|', '[0.025', '0.975]']
allsummaries = pd.DataFrame(columns = summarycols)
oneshot = []
for xplot in xplotlist:
    if ('time' in (xplot)):
        rec = 'Breakthroughtime'
    elif ('fraction' in (xplot)):
        rec = 'fraction'
    print (xplot, rec)
    f1 = yplot + ' ~ ' + xplot
    for g in groups:
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~'+rec)
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        oneshot.append([yplot, xplot,g,RMSE, "RandomFixed"])
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0) 
    for g in groups:
        mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g])
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        oneshot.append([yplot, xplot,g,RMSE, "Norecurse"])
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
        
mydata = dadata[["%reldelmassflux", "PeDa", "Breakthroughtime", "fraction", "del4massflux", "conc_Da", "Pe", "Chem"]]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
mydata.loc[mydata["PeDa"] > 40, "PeDamark"] = 3
mydata.loc[(mydata["PeDa"] > 20) & (mydata["PeDa"] < 40), "PeDamark"] = 2
mydata.loc[(mydata["PeDa"] > 1) & (mydata["PeDa"] < 20), "PeDamark"] = 1
mydata.loc[mydata["PeDa"] < 1, "PeDamark"] = 0
groups = ["conc_Da","Pe", "PeDamark"]
for xplot in xplotlist:
    rec = 'fraction'
    f1 = yplot + ' ~ ' + xplot
    for g in groups:
        mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g] , re_formula = '~'+rec)
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        oneshot.append([yplot, xplot,g,RMSE, "RandomFixed"])
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)
    for g in groups:
        mymd = smf.mixedlm (formula = f1, data = mydata, groups = mydata[g])
        mdf = mymd.fit()
        y_predict = mdf.fittedvalues
        RMSE = np.sqrt(((y-y_predict)**2).mean())
        oneshot.append([yplot, xplot,g,RMSE, "norecurse"])
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0)

results = pd.DataFrame.from_records(oneshot, columns = ["response", "fixed_effect", "random_effect", "RMSE", "Model_Type"])
results.to_csv("X:/RMSE_results_forpaper_v2.csv", sep = ',')
allsummaries.to_csv("X:/Grouped_summary_forpaper_v2.csv", sep = ',')

#Trying sci-kit learn for feature selection
#Transient
allfeatures = dadata.columns.tolist()
features = ['Variance', 'Anisotropy', 'Pe', 'Chem','conc_Da', 'Breakthroughtime', 'fraction']
yfeature = ['del4massflux']

mydata = dadata[features+yfeature]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata["Variance"]= mydata["Variance"].astype(str)
mydata["Anisotropy"]= mydata["Anisotropy"].astype(str)

from sklearn.ensemble import RandomForestRegressor as rf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.svm import SVR
from sklearn.feature_selection import RFE, SelectFromModel

feattofilt = ['Variance', 'Anisotropy', 'Pe', 'conc_Da', 'Breakthroughtime', 'fraction']
X = mydata[feattofilt]
y = mydata[yfeature[0]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))

selector = RFE (regressor, n_features_to_select = 3, step = 1)
selector = selector.fit(X, y)
ranking = selector.ranking_

estimator = SVR(kernel = "linear")
estimator.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
selector = RFE (estimator, n_features_to_select = 3, step = 1)
selector = selector.fit(X, y)
ranking = selector.ranking_

feattofilt = ['Variance', 'Anisotropy', 'Pe', 'conc_Da', 'Breakthroughtime', 'fraction']
X = mydata[feattofilt]
y = mydata[yfeature[0]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
selector = SelectFromModel (regressor).fit(X,y)
thresh = selector.threshold
filt = selector.get_support()

estimator = SVR(kernel = "linear")
selector = SelectFromModel (estimator).fit(X,y)
co = selector.estimator_.coef_
filt = selector.get_support()
