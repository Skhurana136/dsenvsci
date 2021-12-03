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
c_mol_bio = 10 #carbon atoms in microbial necromass
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
    C_recycle = (1-y_params)*r_resp
    C_recycle = np.sum(C_recycle, axis = 1)

    # define each ODE
    # 1. rate of change of carbon species concentration
    C_rate = -1*np.sum(r_resp, axis = 1) + c_mol_bio*r_hydrol/c_n
    for ind in list(range(np.size(simple_C))):
        if ind < np.size(simple_C)-1:
            simp_c_ind = simple_C[ind] 
            comp_c_ind = simple_C[ind+1] 
            C_rate[simp_c_ind] = C_rate[simp_c_ind] + C_recycle[comp_c_ind]
    
    # 2. rate of change of biomass species concentration
    B_rate = r_growth - r_mort
    
    # Checking violation of conditions
    C_bad_ind = np.where(C<=0)
    B_bad_ind = np.where(B<=0)
    C_rate[C_bad_ind] = 0
    B_rate[B_bad_ind] = 0    

    all_rates = np.append(C_rate, B_rate, axis = 0)
    return all_rates

def divers_carbon (data, c_n, b_n):
    C = data[:,:c_n]
    B = data[:,c_n:]
    
    ## SHANNON DIVERSITY AND CARBON STOCK
    Shannon = -np.sum(B/np.sum(B, axis = 0)*np.log(B/np.sum(B, axis = 0)), axis = 1)
    C_stock = np.sum(C,axis=1) + 10*np.sum(B, axis=1)

    return Shannon, C_stock

#%%
## CONSTANTS SAMPLING
#
# Number of DOM/Carbon species:

dom_n = np.random.randint(10,1000,1)[0]
bio_n = np.random.randint(5,100,1)[0]
# Initialize the same number of parameters:

vparams = np.random.uniform(2,5,bio_n*dom_n)
kparams = np.random.randint(1000,8000,bio_n*dom_n)
yparams = np.random.ranf(bio_n)
mparams = np.mean(vparams.reshape(dom_n, bio_n),axis=0)/10

# initial conditions
dom_initial = np.random.randint(1000,10000,dom_n)
biomass_initial = np.random.randint(5,50,bio_n)
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

S, Cstock = divers_carbon(x, dom_n, bio_n)

plt.figure()
plt.plot(S, Cstock)
plt.xlabel ("Shannon")
plt.ylabel ("Carbon")