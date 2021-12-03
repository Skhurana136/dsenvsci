Created on Mon Nov 29 2021

# Creator: Swamini Khurana <swamini.khurana@gmail.com>
"""

"""
# %%
## Import libraries
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
## CONSOLIDATE SAMPLED PARAMETERS
def call_constants():
    params = [vparams, kparams, yparams, mparams]
    return params
# %%
def rates(x,t, c_n, b_n):
    
    # constants
    params = call_constants()
    v_params = params[0]
    k_params = params[1]
    y_params = params[2]
    m_params = params[3]
    c_mol_bio = 10 #carbon atoms in microbial necromass

    # assign each ODE to array elements
    C = x[:c_n]
    B = x[c_n:]

    # define rates for:
    # 1. respiration
    r_resp = v_params.reshape(c_n,b_n)*np.multiply(B,C[...,None])/(k_params.reshape(c_n, b_n)+C[...,None])
    # 2. growth
    r_growth = y_params*np.sum(r_resp,axis=0)
    # 3. mortality
    r_mort = m_params*B
    # 4. dead microbes back to C (distributed evenly across all C)
    r_hydrol = 0.1 * sum(r_mort)
    # 5. Complex DOM species added to simpler DOM species pools
    simple_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))

    # define each ODE
    # 1. rate of change of carbon species concentration
    C_rate = -1*np.sum(r_resp, axis = 1) + c_mol_bio*r_hydrol/c_n
    for ind in list(range(np.size(simple_C))):
        if ind < np.size(simple_C)-1:
            simp_c_ind = simple_C[ind] 
            comp_c_ind = simple_C[ind+1] 
            C_rate[simp_c_ind] = C_rate[simp_c_ind] + np.sum(r_resp[comp_c_ind])
    
    # 2. rate of change of biomass species concentration
    B_rate = r_growth - r_mort
    
    # Checking violation of conditions
    C_bad_ind = np.where(C<=0)
    B_bad_ind = np.where(B<=0)
    C_rate[C_bad_ind] = 0
    B_rate[B_bad_ind] = 0    

    all_rates = np.append(C_rate, B_rate, axis = 0)
    return all_rates

#%%
## CONSTANTS SAMPLING
#
# Number of DOM/Carbon species:

dom_n = np.random.randint(3,10,1)[0]
bio_n = np.random.randint(2,10,1)[0]
# Initialize the same number of parameters:

vparams = np.random.randint(2,10,bio_n*dom_n)
kparams = np.random.randint(1000,8000,bio_n*dom_n)
yparams = np.random.ranf(bio_n)
mparams = np.mean(vparams.reshape(dom_n, bio_n),
axis=0)/10

# initial conditions
dom_initial = np.random.randint(100,1000,dom_n)
biomass_initial = np.random.randint(5,100,bio_n)
x0 = np.append(dom_initial,biomass_initial)

# declare a time vector (time window)
t = np.arange(0,100,0.01)
# solve
x = odeint(rates, x0, t, args=(dom_n,bio_n))#, full_output=1)

C = x[:,:dom_n]
B = x[:,dom_n:]
plt.plot(t,C, linestyle = "-", label = "DOM")
plt.plot(t,B, linestyle = "--", label = "Biomass")
#plt.legend()

#%%
## SHANNON DIVERSITY AND CARBON STOCK
S = -np.sum(B/np.sum(B, axis = 0)*np.log(B/np.sum(B, axis = 0)), axis = 1)
Cresp = np.sum(C,axis=1)

plt.figure()
plt.plot(S, Cresp)
plt.xlabel ("Shannon")
plt.ylabel ("Carbon")

#%%
simC = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))
for ind in list(range(np.size(simC))):
        if ind < np.size(simC)-1:
            simp_c_ind = simC[ind] 
            comp_c_ind = simC[ind+1] 
            print(simp_c_ind, comp_c_ind)
            C[simp_c_ind] = C[simp_c_ind] + C[comp_c_ind]