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

from scipy.integrate import odeint

class ReactionNetwork(object):
    """Class to define general information on the reaction network solver, such as
    rate constants, rate expressions."""

    @classmethod
    def set_constants(self, *user_input):
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
        #return parameterslist, num_para_sets
    
    @classmethod
    def solve_network (self, initial_conditions, time_space, c_n, b_n):
        """Method to set up the reaction network given the 
        number of dissolved organic matter species and microbial
        species.
        
        Parameter
        ---------
        c_n, b_n : int.
            c_n : number of carbon species.
            b_n : number of microbial species.
        """
        def rate_eqn(x , t, c_n, b_n):
            """Method to set up the reaction network given the 
            number of dissolved organic matter species and microbial
            species.
            
            Parameter
            ---------
            c_n, b_n : int.
                c_n : number of carbon species.
                b_n : number of microbial species.
            """
    
            # constants
            #params, num_para_sets = call_constants()
            v_params = self.parameters[0]
            k_params = self.parameters[self.para_size-3]
            y_params = self.parameters[self.para_size-2]
            m_params = self.parameters[self.para_size-1]

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

            rates = np.append(C_rate, B_rate, axis = 0)

            return rates

        # solve
        #rates_nw = form_rate_expressions(initial_conditions, time_space, c_n, b_n)
        self.answer = odeint (rate_eqn, initial_conditions, time_space, args=(c_n, b_n))

        return self.answer

    def diversity_carbon (self, c_n, b_n):
        """Method to evaluate the Shannong diversity and
        carbon stock in the domain.
            
            Parameters
            ---------
            c_n, b_n : int.
                c_n : number of carbon species.
                b_n : number of microbial species.
            """
        C = self.answer[:,:c_n]
        B = self.answer[:,c_n:]
        
        ## SHANNON DIVERSITY AND CARBON STOCK
        proportion = B/np.sum(B,axis=0)
        Shannon = -np.sum(proportion*np.log(proportion), axis = 1)
        total_C_stock = np.sum(C,axis=1) + np.sum(B, axis=1)
        C_stock = np.sum(C,axis=1)

        return Shannon, C_stock, total_C_stock
