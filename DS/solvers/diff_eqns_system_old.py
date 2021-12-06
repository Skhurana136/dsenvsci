Created on Mon Nov 29 2021

# Creator: Swamini Khurana <swamini.khurana@gmail.com>
"""

"""
## Import libraries
from scipy.integrate import odeint
import numpy as np

from cfg import *

def rates(x,c_n, b_n):
    
    # constants
    #params, num_para_sets = call_constants()
    v_params = params[0]
    k_params = params[num_para_sets-3]
    y_params = params[num_para_sets-2]
    m_params = params[num_para_sets-1]

    # assign each ODE to array elements
    C = x[:c_n]
    B = x[c_n:]

    # define rates for:
    # 1. respiration
    # formula to implement for all B,C combinations:
    # r = v*B*C/(k+C). v,k are unique to each B,C combination
    r_resp = v_params.reshape(c_n,b_n)*np.multiply(B,C[...,None])/(k_params.reshape(c_n, b_n)+C[...,None])
    # 2. growth
    # formula to implement for all B
    # r = y* rate of respiration for each B for all C
    r_growth = y_params*np.sum(r_resp,axis=0)
    # 3. mortality
    # formula to implement for all B
    # r = m*B
    r_mort = m_params*B
    # 4. dead microbes back to C (distributed evenly across all C)
    # formula to implement for all dying B and add it to C
    # r = sum of dying biomass
    r_hydrol = 0.5*sum(r_mort)/c_n
    # 5. Complex DOM species added to simpler DOM species pools
    # formula to implement for all C
    # locate and order the C in decreasing order of simplicity
    simple_C = np.argsort(np.mean(k_params.reshape(c_n, b_n), axis=1))
    # Add non-absorbed C fraction (not used for growth) into recycling fraction
    # r = (1-y)*rate of respiration which will be unique for each B,C combination
    C_recycle = (1-y_params)*r_resp
    # In the end we are interested in the bulk unused C product to recycle it to other C pools:
    # r = sum(recycled pool for each C)
    C_recycle = np.sum(C_recycle, axis = 1)

    # define each ODE
    # 1. rate of change of carbon species concentration
    C_rate = -1*np.sum(r_resp, axis = 1) + r_hydrol
    for ind in list(range(np.size(simple_C))):
        if ind < np.size(simple_C)-1: # No action to be taken for the most complex compound
            # identify the index in C_rate for the simple C compound
            simp_c_ind = simple_C[ind] 
            # identify the next level complex compound and its index
            comp_c_ind = simple_C[ind+1]
            # Add the unused fraction of next level complex compound to the simpler C compound
            C_rate[simp_c_ind] += C_recycle[comp_c_ind]
    # save microbial necromass into the carbon pool that is midway simple:
    #mid_ind = simple_C[int(c_n/2)]
    #C_rate[mid_ind] += r_hydrol

    # 2. rate of change of biomass species concentration
    # Biomsdd grows and dies
    B_rate = r_growth - r_mort
    
    # Checking violation of conditions
    # Check if concentration of carbon and biomass is tending to negative. If yes, reset it to 0.
    C_bad_ind = np.where(C<=0)
    B_bad_ind = np.where(B<=0)
    C_rate[C_bad_ind] = 0
    B_rate[B_bad_ind] = 0    

    all_rates = np.append(C_rate, B_rate, axis = 0)
    return all_rates

def solve_network (rate_eqns, initial_conditions, time_space, other_args):
    # solve
    answer = odeint(rate_eqns, initial_conditions, time_space, args=other_args)

    return answer

def divers_carbon (data, c_n, b_n):
    C = data[:,:c_n]
    B = data[:,c_n:]
    
    ## SHANNON DIVERSITY AND CARBON STOCK
    proportion = B/np.sum(B,axis=0)
    Shannon = -np.sum(proportion*np.log(proportion), axis = 1)
    total_C_stock = np.sum(C,axis=1) + np.sum(B, axis=1)
    C_stock = np.sum(C,axis=1)

    return Shannon, C_stock, total_C_stock