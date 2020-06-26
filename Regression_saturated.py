# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 09:54:53 2020

@author: khurana
"""
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Regression or other methods of analyses
mfile = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/massflux_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)
bfile = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)

mfile["Pe"] = mfile["Regime"]
mfile["Pe"] = mfile["Pe"].replace(["Slow", "Medium", "Fast"], [2.02, 11.7, 22.45])
bfile["Pe"] = bfile["Regime"]
bfile["Pe"] = bfile["Pe"].replace(["Slow", "Medium", "Fast"], [2.02, 11.7, 22.45])

data = mfile[mfile["Chem"] == "TOC"]
y = data["del2massflux"]
X = data[["fraction", "Pe"]]

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

reg = LinearRegression().fit(X_train, y_train)
plt.plot(X_test, reg.predict(X_test), label="linear regression")
plt.plot(X_test, y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# Score
score = LinearRegression().fit(X_train, y_train).score(X_test, y_test)
print("Test score: {:.3f}".format(score))

# Transformation by taking log:
X_train["Pe_log"] = np.log10(X_train["Pe"])
X_test["Pe_log"] = np.log10(X_test["Pe"])
X_train_log = X_train[["fraction", "Pe_log"]]
X_test_log = X_test[["fraction", "Pe_log"]]
plt.hist(np.log(X_train_log), bins=5)
plt.ylabel("Number of appearances")
plt.xlabel("Value")
# Regression over the log transformed values
from sklearn.linear_model import Ridge

reg = Ridge().fit(X_train_log, y_train)
# SUBSTANTIAL IMPROVEMENT AS EXPECTED
plt.plot(X_test_log["fraction"], reg.predict(X_test_log), label="linear regression")
plt.plot(X_test_log["fraction"], y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# Score
score = Ridge().fit(X_train_log, y_train).score(X_test_log, y_test)
print("Test score: {:.3f}".format(score))
# no improvement in score
from sklearn.preprocessing import PolynomialFeatures

poly = PolynomialFeatures(degree=4).fit(X_train_log)
X_train_log_poly = poly.transform(X_train_log)
X_test_log_poly = poly.transform(X_test_log)

print("X_train.shape: {}".format(X_train_log.shape))
print("X_train_poly.shape: {}".format(X_train_log_poly.shape))
# names of features:
print("Polynomial feature names:\n{}".format(poly.get_feature_names()))
reg = Ridge().fit(X_train_log_poly, y_train)
score = reg.score(X_test_log_poly, y_test)
print("Test score: {:.3f}".format(score))

plt.plot(
    np.sort(X_test_log_poly[:, 1], axis=0),
    reg.predict(np.sort(X_test_log_poly, axis=0)),
    label="linear regression",
)
plt.plot(X_test_log_poly[:, 1], y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# score improved. but it is still at 0.345 ...

# Polynomial regression without taking log:
X_train = X_train[["fraction", "Pe"]]
X_test = X_test[["fraction", "Pe"]]
poly = PolynomialFeatures(degree=4).fit(X_train)
X_train_poly = poly.transform(X_train)
X_test_poly = poly.transform(X_test)
print("X_train.shape: {}".format(X_train.shape))
print("X_train_poly.shape: {}".format(X_train_poly.shape))
reg = Ridge().fit(X_train_poly, y_train)
score = reg.score(X_test_log_poly, y_test)
print("Test score: {:.3f}".format(score))
plt.plot(
    np.sort(X_test_poly[:, 1], axis=0),
    reg.predict(np.sort(X_test_poly, axis=0)),
    label="linear regression",
)
plt.plot(X_test_poly[:, 1], y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# score falls below 0 into negatives. this did not help

# Going back to transformation but doing an exponential transformation of PE this time
X_train["Pe_exp"] = np.exp(X_train["Pe"])
X_test["Pe_exp"] = np.exp(X_test["Pe"])
X_train_exp = X_train[["fraction", "Pe_exp"]]
X_test_exp = X_test[["fraction", "Pe_exp"]]
plt.hist(X_train_exp, bins=5)
plt.ylabel("Number of appearances")
plt.xlabel("Value")
# Regression over the log transformed values
reg = Ridge().fit(X_train_exp, y_train)
# SUBSTANTIAL IMPROVEMENT AS EXPECTED
plt.plot(
    np.sort(X_test_exp["fraction"]),
    reg.predict(np.sort(X_test_exp)),
    label="linear regression",
)
plt.plot(X_test_exp["fraction"], y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# Score
score = reg.score(X_test_exp, y_test)
print("Test score: {:.3f}".format(score))
# the score is bad at 0.167
from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler()
X_train = X_train[["fraction", "Pe"]]
X_test = X_test[["fraction", "Pe"]]
scaler.fit(X_train)
X_train_scaled = scaler.transform(X_train)
# print dataset properties before and after scaling
print("transformed shape: {}".format(X_train_scaled.shape))
print("per-feature minimum before scaling: \n{}".format(X_train.min(axis=0)))
print("per-feature maximum before scaling: \n{}".format(X_train.max(axis=0)))
print("per-feature minimum after scaling: \n{}".format(X_train_scaled.min(axis=0)))
print("per-feature maximum after scaling: \n{}".format(X_train_scaled.max(axis=0)))

# Apply SVM to the scaled data, transform the test set.
# transform test data
X_test_scaled = scaler.transform(X_test)
print("per-feature minimum after scaling: \n{}".format(X_test_scaled.min(axis=0)))
print("per-feature maximum after scaling: \n{}".format(X_test_scaled.max(axis=0)))

plt.hist(X_train_scaled, bins=5)
plt.ylabel("Number of appearances")
plt.xlabel("Value")
# Regression over the log transformed values
reg = LinearRegression().fit(X_train_scaled, y_train)
# SUBSTANTIAL IMPROVEMENT AS EXPECTED
plt.plot(
    np.sort(X_test_scaled[:, 0]),
    reg.predict(np.sort(X_test_scaled)),
    label="linear regression",
)
plt.plot(X_test_scaled, y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# Score
score = reg.score(X_test_scaled, y_test)
print("Test score: {:.3f}".format(score))

# Transformation
X_train_prod = np.array(X_train["fraction"] * X_train["Pe"])
X_test_prod = np.array(X_test["fraction"] * X_test["Pe"])
plt.hist(X_train_prod, bins=5)
plt.ylabel("Number of appearances")
plt.xlabel("Value")
# Regression over the log transformed values
reg = LinearRegression().fit(X_train_prod.reshape(-1, 1), y_train)
# SUBSTANTIAL IMPROVEMENT AS EXPECTED
plt.plot(
    np.sort(X_test_prod),
    reg.predict(np.sort(X_test_prod.reshape(-1, 1))),
    label="linear regression",
)
plt.plot(X_test_prod, y_test, "o", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
# Score
score = reg.score(X_test_prod.reshape(-1, 1), y_test)
print("Test score: {:.3f}".format(score))
# Regression is not working. Try auto-encoder
bins = np.linspace(0, 30, 4)
line = np.linspace(0, 30, 1000, endpoint=False).reshape(-1, 1)
print("bins: {}".format(bins))
which_bin = (np.digitize(X_train["Pe"], bins=bins)).reshape(-1, 1)
from sklearn.preprocessing import OneHotEncoder

encoder = OneHotEncoder(sparse=False)
encoder.fit(which_bin)
X_binned = encoder.transform(which_bin)
print(X_binned[:5])
line_binned = encoder.transform(np.digitize(line, bins=bins))

X_combined = np.hstack([X_train, X_binned])
print(X_combined.shape)

reg = LinearRegression().fit(X_combined, y_train)
line_combined = np.hstack([line, line_binned])
plt.plot(
    np.sort(X_combined),
    reg.predict(np.sort(X_combined)),
    label="linear regression combined",
)

for bin in bins:
    plt.plot([bin, bin], [0, 1], ":", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.plot(X[:, 0], y, "o", c="k")

X_product = (X_binned.T * np.array(X_train["fraction"])).T
reg0 = LinearRegression().fit(X_product[:, 0].reshape(-1, 1), y_train)
reg1 = LinearRegression().fit(X_product[:, 1].reshape(-1, 1), y_train)
reg2 = LinearRegression().fit(X_product[:, 2].reshape(-1, 1), y_train)

plt.plot(
    np.sort(X_test["fraction"]),
    reg0.predict(np.sort(np.array(X_test["fraction"]).reshape(-1, 1))),
    label="Reg0",
)
plt.plot(
    np.sort(X_test["fraction"]),
    reg1.predict(np.sort(np.array(X_test["fraction"]).reshape(-1, 1))),
    label="Reg1",
)
plt.plot(
    np.sort(X_test["fraction"]),
    reg2.predict(np.sort(np.array(X_test["fraction"]).reshape(-1, 1))),
    label="Reg2",
)
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.plot(X[:, 0], y, "o", c="k")

# Nothing worked. Going back and looking at the data afresh
mfmaster = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/massflux_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)
mfmaster["Pe"] = mfmaster["Regime"]
mfmaster["Pe"] = mfmaster["Pe"].replace(["Slow", "Medium", "Fast"], [2.02, 11.7, 22.45])
bmaster = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)

red_patch = mpatches.Patch(color="indianred", label="Slow flow")
green_patch = mpatches.Patch(color="g", label="Medium flow")
blue_patch = mpatches.Patch(color="steelblue", label="Fast flow")

Regimes = ["Fast", "Medium", "Slow"]
colors = ["steelblue", "g", "indianred"]
coefficients = []
intercepts = []
masterscores = []
res = []
chems = []
fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=[11, 8])
i = 0
for chemical in ["DOC", "DO", "Nitrogen", "TOC"]:
    print(chemical)
    mfile = mfmaster[mfmaster["Chem"] == chemical]
    paxis = axes.flat[i]
    handles = []
    scores = []
    for R in Regimes:
        print(R)
        X = mfile[mfile["Regime"] == R][["fraction"]]
        y = mfile[mfile["Regime"] == R][["del2massflux"]]
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
        reg = LinearRegression().fit(X_train, y_train)
        score = np.round(reg.score(X_test, y_test), 2)
        scores.append(score)
        coefficients.append(reg.coef_[0])
        intercepts.append(reg.intercept_)
        masterscores.append(score)
        res.append(R)
        chems.append(chemical)
        l1, = paxis.plot(
            np.sort(X), reg.predict(np.sort(X)), label=score, c=colors[Regimes.index(R)]
        )
        handles.append(l1)
        paxis.plot(X, y, ".", c=colors[Regimes.index(R)])
    paxis.set_title(chemical)
    paxis.legend(handles, scores, loc="best", title="Score", fontsize=12)
    i += 1
# plt.legend(title = "Score",loc='best')
fig.legend(
    handles=[red_patch, green_patch, blue_patch],
    loc="best",
    title="Flow regime",
    fontsize=12,
)
plt.annotate(
    "Fraction of removal of nutrient in base case",
    xy=(-1.2, 1),
    xytext=(-70, 0),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="left",
    va="center",
    rotation="vertical",
    fontsize=15,
)
plt.annotate(
    "Fraction of breakthrough time in base case",
    xy=(0, 0),
    xytext=(0, -40),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="center",
    va="center",
    fontsize=15,
)
plt.suptitle("Regression curve for nutrient removal in the domain", fontsize=15)
results = pd.DataFrame(
    np.hstack(
        (
            np.array(coefficients),
            np.array(intercepts),
            np.array(masterscores).reshape(-1, 1),
            np.array(chems).reshape(-1, 1),
            np.array(res).reshape(-1, 1),
        )
    ),
    columns=["Coefficient", "Intercept", "Score", "Chem", "Regime"],
)
plt.savefig(
    "Z:/Saturated_flow/diffusion_transient/steady_state_mass_flux_impact_regression.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.01,
)
plt.savefig(
    "Z:/Saturated_flow/diffusion_transient/steady_state_mass_flux_impact_regression.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.01,
)


bmaster = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)
bmaster["Pe"] = bmaster["Regime"]
bmaster["Pe"] = bmaster["Pe"].replace(["Slow", "Medium", "Fast"], [2.02, 11.7, 22.45])
Regimes = ["Fast", "Medium", "Slow"]
colors = ["steelblue", "g", "indianred"]
coefficients = []
intercepts = []
masterscores = []
res = []
chems = []
Chemicalspecies = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
]
species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
position = ["Immobile", "Mobile"]
fig, axes = plt.subplots(
    nrows=2,
    ncols=3,
    sharex=True,
    figsize=[len(Chemicalspecies) * 2, len(Chemicalspecies)],
)
i = 0
for chemical in Chemicalspecies:
    print(chemical)
    mfile = bmaster[bmaster["Chem"] == chemical]
    paxis = axes.flat[i]
    handles = []
    scores = []
    for R in Regimes:
        print(R)
        X = mfile[mfile["Regime"] == R][["fraction"]]
        y = mfile[mfile["Regime"] == R][["Change_umoles"]]
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
        reg = LinearRegression().fit(X_train, y_train)
        score = np.round(reg.score(X_test, y_test), 2)
        scores.append(score)
        coefficients.append(reg.coef_[0])
        intercepts.append(reg.intercept_)
        masterscores.append(score)
        res.append(R)
        chems.append(chemical)
        l1, = paxis.plot(
            np.sort(X), reg.predict(np.sort(X)), label=score, c=colors[Regimes.index(R)]
        )
        handles.append(l1)
        paxis.plot(X, y, ".", c=colors[Regimes.index(R)])
    paxis.legend(handles, scores, loc="best", title="Score", fontsize=12)
    i += 1
plt.legend(loc="best", fontsize=12)
for ax, typsp in zip(axes[0, :], species):
    ax.set_title(typsp, fontsize=15)
axes[0, -1].annotate(
    position[0],
    xy=(0, 0.5),
    xytext=(200, 0),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="left",
    va="center",
    rotation="vertical",
    fontsize=15,
)
axes[1, -1].annotate(
    position[1],
    xy=(0, 0.5),
    xytext=(200, 0),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="left",
    va="center",
    rotation="vertical",
    fontsize=15,
)
fig.legend(
    handles=[red_patch, green_patch, blue_patch],
    loc="best",
    title="Flow regime",
    fontsize=12,
)
plt.annotate(
    "Normalized biomass in the domain (uM C)",
    xy=(-2.5, 1),
    xytext=(-50, 0),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="left",
    va="center",
    rotation="vertical",
    fontsize=15,
)
plt.annotate(
    "Fraction of breakthrough time in base case",
    xy=(-0.4, -0.2),
    xytext=(0, 0),
    xycoords="axes fraction",
    textcoords="offset points",
    size="large",
    ha="center",
    va="center",
    fontsize=15,
)
# plt.suptitle("Regression curve for nutrient removal in the domain", fontsize = 15)
results = pd.DataFrame(
    np.hstack(
        (
            np.array(coefficients),
            np.array(intercepts),
            np.array(masterscores).reshape(-1, 1),
            np.array(chems).reshape(-1, 1),
            np.array(res).reshape(-1, 1),
        )
    ),
    columns=["Coefficient", "Intercept", "Score", "Chem", "Regime"],
)
plt.savefig(
    "Z:/Saturated_flow/diffusion_transient/steady_state_biomass_impact_regression.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.01,
)
plt.savefig(
    "Z:/Saturated_flow/diffusion_transient/steady_state_biomass_impact_regression.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.01,
)


# Tryingencoder again
mfmaster["product"] = mfmaster["fraction"] * mfmaster["Pe"]
mfmaster["product"] = np.log10(mfmaster["product"])
bins = np.linspace(0, 1.5, 11)
line = np.linspace(0, 1.5, 1000, endpoint=False).reshape(-1, 1)
print("bins: {}".format(bins))
X = mfmaster[mfmaster["Chem"] == "TOC"]
y = mfmaster[mfmaster["Chem"] == "TOC"]
X_train, X_test, y_train, y_test = train_test_split(
    X["product"], y["del2massflux"], random_state=0
)
which_bin = (np.digitize(X_train, bins=bins)).reshape(-1, 1)
from sklearn.preprocessing import OneHotEncoder

encoder = OneHotEncoder(sparse=False)
encoder.fit(which_bin)
X_binned = encoder.transform(which_bin)
print(X_binned[:5])
line_binned = encoder.transform(np.digitize(line, bins=bins))
X_combined = np.hstack([np.array(X_train).reshape(-1, 1), X_binned])
print(X_combined.shape)

reg = LinearRegression().fit(X_combined, y_train)
line_combined = np.hstack([line, line_binned])
plt.plot(
    line_combined[:, 0], reg.predict(line_combined), label="linear regression combined"
)
plt.plot(X_combined[:, 0], y_train, ".", c="k")
for bin in bins:
    plt.plot([bin, bin], [0, 1], ":", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.plot(X[:, 0], y, "o", c="k")
X_test_binned = encoder.transform((np.digitize(X_test, bins=bins)).reshape(-1, 1))
X_test_combined = np.hstack([np.array(X_test).reshape(-1, 1), X_test_binned])
reg.score(X_test_combined, y_test)

reg = LinearRegression().fit(np.array(X_train).reshape(-1, 1), y_train)
plt.plot(
    X_train,
    reg.predict(np.array(X_train).reshape(-1, 1)),
    label="linear regression combined",
)
plt.plot(X_train, y_train, ".", c="k")
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.plot(X_train, y_train, "o", c="k")
reg.score(np.array(X_test).reshape(-1, 1), y_test)

# Switch to SVM
mfmaster = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/massflux_withbreakthrough_forMartin_v4_complete.csv",
    sep="\t",
)
mfmaster["Pe"] = mfmaster["Regime"]
mfmaster["Pe"] = mfmaster["Pe"].replace(["Slow", "Medium", "Fast"], [2.02, 11.7, 22.45])
mfmaster["product"] = mfmaster["fraction"] / mfmaster["Pe"]
mfmaster["product"] = np.log10(mfmaster["product"])
X = mfmaster[mfmaster["Chem"] == "TOC"][["product"]]
y = mfmaster[mfmaster["Chem"] == "TOC"][["del2massflux"]]
X_train, X_test, y_train, y_test = train_test_split(
    np.array(X), np.array(y), random_state=0
)

line = np.linspace(np.min(X), np.max(X), 1000)
from sklearn.svm import SVR

for gamma in [0.5, 1, 5, 10]:
    svr = SVR(gamma=gamma).fit(np.array(X_train).reshape(-1, 1), y_train)
    print(svr.score(X_test.reshape(-1, 1), y_test))
    plt.plot(
        line.reshape(-1, 1),
        svr.predict(line.reshape(-1, 1)),
        label="SVR gamma = {}".format(gamma),
    )

plt.plot(X, y, ".", c="k")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.legend(loc="best")
