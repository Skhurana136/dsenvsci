"""
# -*- coding: utf-8 -*-
# ======================================================================
# Created by : Swamini Khurana
# Created on : On Mon Dec 05 2021 at 00:32:06
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

from scipy.integrate import odeint, solve_bvp

class ReactionNetwork(object):
    """Class to define general information on the reaction network solver, such as
    rate constants, rate expressions."""

    def __init__(self,
        maximum_capacity = 300,
        carbon_mol_bio = 10,
        carbon_num = 1,
        bio_num = 1,
        carbon_input = 0,
        sigmoid_coeff_stolpovsky = 0.1,
        ):

        """Method to assign constants for the entire system.
        
        Parameter
        ---------
        maximum_capacity : float.
            Define the maximum capacity of the system.
            If no value is specified then default is 300.
        carbon_mol_bio : float.
            Define the carbon atoms in biomolecules.
            If no value is specified then default is 10.
        carbon_num : int.
            Number of carbon species. Default value is 1.
        bio_num : int.
            Number of microbial species. Default value is 1.
        carbon_input : float.
            Continuous addition of carbon to the most complex carbon pool. Default is 0, implying a closed system.
        sigmoid_coeff_stolpovsky : float.
            Coefficient to use in adapted sigmoid curve (see Stolpovsky et al., 2011). Can vary from 0.01 to 0.1.
        """
        
        self.max_cap = maximum_capacity
        self.c_mol_bio = carbon_mol_bio
        self.c_n = carbon_num
        self.b_n = bio_num
        self.c_bc = carbon_input
        self.st = sigmoid_coeff_stolpovsky
        
    def set_rate_constants(self, *user_input):
        """Method to assign constants to pass to other functions.
        
        Parameter
        ---------
        user_input : Nested list.
            Nested list to use for the reaction network.
            Each list must be one particular type of constant. 
            List 1: maximum rate constants.
            List 2: Half saturation rate constants.
            List 3: Yield coefficients.
            List 4: Mortality etc.
        """
        
        self.parameters = list(l for l in user_input)
        self.para_size = len(self.parameters)
        print("Reaction network constants have been defined")
    
    def solve_network (self, initial_conditions, time_space):
        """Method to set up the reaction network given the 
        number of dissolved organic matter species and microbial
        species.
        
        Parameter
        ---------
        initial_conditions : Array, float.
            Initial conditions to solve the system of ODEs.
            Array of the same size as variables to be solved.
        time_space : Array, float.
            Array of time points to solve the system of ODEs.
        """

        def rate_expressions(x , t):
            """Function to set up the reaction network given the 
            number of dissolved organic matter species and microbial
            species.
            """
    
            # constants
            v_params = self.parameters[0]
            k_params = self.parameters[self.para_size-3]
            y_params = self.parameters[self.para_size-2]
            m_params = self.parameters[self.para_size-1]
            lim_k_params = k_params/50

            # assign each ODE to array elements
            C = x[:self.c_n]
            B = x[self.c_n:]

            # define rates for:
            # 1. respiration
            # formula to implement for all B,C combinations:
            # r = v*B*C/(k+C). v,k are unique to each B,C combination
            r_resp = v_params.reshape(self.c_n,self.b_n)*np.multiply(B,C[...,None])/(k_params.reshape(self.c_n, self.b_n)+C[...,None])
            # 2. growth
            # formula to implement for all B
            # r = y* rate of respiration for each B for all C
            r_growth = y_params*np.sum(r_resp/np.exp((1-(r_resp/(lim_k_params.reshape(self.c_n, self.b_n)*B)))/self.st),axis=0)
            # 3. mortality
            # formula to implement for all B
            # r = m*B
            r_mort = m_params*B
            # 4. dead microbes back to C (distributed evenly across all C)
            # formula to implement for all dying B and add it to C
            # r = sum of dying biomass
            r_hydrol = sum(r_mort)/self.c_n
            # 5. Complex DOM species added to simpler DOM species pools
            # formula to implement for all C
            # locate and order the C in increasing order of complexity (ascending k_params)
            simple_C = np.argsort(np.mean(k_params.reshape(self.c_n, self.b_n), axis=1))
            # Add non-absorbed C fraction (not used for growth) into recycling fraction
            # r = (1-y)*rate of respiration which will be unique for each B,C combination
            C_recycle = r_resp - r_growth#(1-y_params)*np.sum(r_resp/np.exp((1-(r_resp/(lim_k_params.reshape(self.c_n, self.b_n)*B)))/self.st),axis=0)#(1-y_params)*r_resp
            # In the end we are interested in the bulk unused C product to recycle it to other C pools:
            # r = sum(recycled pool for each C)
            C_recycle = np.sum(C_recycle, axis = 1)

            # define each ODE
            # 1. rate of change of carbon species concentration
            C_rate = -1*np.sum(r_resp, axis = 1) + r_hydrol
            # add sequence of addition to the simpler carbon compounds
            for ind in list(range(np.size(simple_C))):
                if ind < np.size(simple_C)-1: # No action to be taken for the most complex compound
                    # identify the index in C_rate for the simple C compound
                    simp_c_ind = simple_C[ind] 
                    # identify the next level complex compound and its index
                    comp_c_ind = simple_C[ind-1]
                    # Add the unused fraction of next level complex compound to the simpler C compound
                    C_rate[simp_c_ind] += C_recycle[comp_c_ind]
            # save microbial necromass into the carbon pool that is midway simple:
            #mid_ind = simple_C[int(self.c_n/2)]
            #C_rate[mid_ind] += r_hydrol
            # add continuous input to the most complex carbon compound
            C_rate[simple_C[-1]] += self.c_bc

            # 2. rate of change of biomass species concentration
            # Biomsdd grows and dies
            B_rate = r_growth - r_mort
        
            # Checking violation of conditions
            # Check if concentration of carbon and biomass is tending to negative. If yes, reset it to 0.
            C_bad_ind = np.where(C<=0)
            B_bad_ind = np.where(B<=0)
            C_rate[C_bad_ind] = 0
            B_rate[B_bad_ind] = 0

            # Check if biomass is exceeding maximum capacity.
            # If yes, reset to maximum capacity and stop changing biomass.
            B_bad_ind = np.where(np.sum(B)>self.max_cap)
            B_rate[B_bad_ind] = 0

            rates = np.append(C_rate, B_rate, axis = 0)

            return rates
        
        ### WIP
        ## Defining boundary conditions/input conditions. This is not working right now.
        def boundary_conditions(x_bc0, x_bc1):
            """Method to set up the boundary conditions given the 
            time point and values.
            
            Parameter
            ---------
            x_bc0, x_bc1 : Array, float.
                x_bc0 : Boundary conditions 1.
                x_bc1 : Boundary conditions 2.
            """
            # Values at t=0:
            x_bc0[:self.c_n] = np.ones(self.c_n)*5
            x_bc0[self.c_n:] = np.ones(self.b_n)*1

            # Values at t=0:
            x_bc1[:self.c_n] = np.ones(self.c_n)*10
            x_bc1[self.c_n:] = np.ones(self.b_n)*5

            # These return values are what we want to be 0:
            return [x_bc0, x_bc1]
        ### End of WIP

        # solve
        self.initial_guess = odeint (rate_expressions, initial_conditions, time_space)
    

        #self.answer = solve_bvp(lambda t, x: rate_expressions(t, x, self.c_n, self.b_n),
        #            lambda x_bc0, x_bc1: boundary_conditions(x_bc0, x_bc1, self.c_n, self.b_n), time_space, self.initial_guess.T)

        return self.initial_guess#, self.answer

def diversity_carbon (data, c_n, b_n):

    """Method to evaluate the Shannon diversity and
    carbon stock in the domain.

    """

    C = data[:,:c_n]
    B = data[:,c_n:]
    
    ## SHANNON DIVERSITY AND CARBON STOCK
    proportion = B/np.sum(B,axis=0)
    Shannon = -np.sum(proportion*np.log(proportion), axis = 1)
    total_C_stock = np.sum(C,axis=1) + np.sum(B, axis=1)
    C_stock = np.sum(C,axis=1)

    return Shannon, C_stock, total_C_stock
