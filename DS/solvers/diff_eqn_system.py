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

from scipy.integrate import solve_ivp, solve_bvp

class ReactionNetwork(object):
    """Class to define general information on the reaction network solver, such as
    rate constants, rate expressions."""

    def __init__(self,
        maximum_capacity = 300,
        carbon_mol_bio = 10,
        carbon_num = 1,
        bio_num = 3,
        fungal_num = 0,
        carbon_input = 0,
        sigmoid_coeff_stolpovsky = 0.1,
        necromass_distribution = "equally",
        enzyme_production_rate_constant = 0.7,
        efficiency_bio_uptake = 0.5,
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
            Number of mediators of carbon transformation in the system. Default value is 3.
        fungal_num : int.
            Number of fungal species. Default value is 2.
        carbon_input : float.
            Continuous addition of carbon to all carbon pools. Default is 0, implying a closed system.
        sigmoid_coeff_stolpovsky : float.
            Coefficient to use in adapted sigmoid curve (see Stolpovsky et al., 2011). Can vary from 0.01 to 0.1.
        necromass_distribution : string.
            Switch to specify if bacteria necromass distributes evenly among all carbon pools or
            "mid-labile" carbon pools. Default option is equal distribution among all carbon pools.
        enzyme_production_rate_constant : float.
            First order rate constant for production of extracellular enzymes by fungi. Default value is 0.7
        efficiency_bio_uptake : float.
            Fraction of depolymerized carbon pool that is used for microbial uptake (respiration + growth). Default value is 0.5.
        """
        
        self.max_cap = maximum_capacity
        self.c_mol_bio = carbon_mol_bio
        self.c_n = carbon_num
        self.f_n = fungal_num
        self.b_n = bio_num
        self.c_bc = carbon_input
        self.st = sigmoid_coeff_stolpovsky
        self.necromass_loop = necromass_distribution
        self.v_enz = enzyme_production_rate_constant
        self.z = efficiency_bio_uptake
        
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

    def identify_components_natures(self, recalcitrance_criterion):
        """Function to categorize the components of the reaction network.
        What is most recalcitrant (least labile) carbon pool?
        Which microbial species are fungi and the rest bacteria?
        
        Parameter
        ---------

        """
        
        # constants
        self.v_params = self.parameters[0].reshape(self.c_n, self.b_n)
        self.k_params = self.parameters[self.para_size-3].reshape(self.c_n, self.b_n)
        self.y_params = self.parameters[self.para_size-2]
        self.m_params = self.parameters[self.para_size-1]
        #self.lim_k_params = self.k_params/100

        if recalcitrance_criterion == 'oxidation_state':
            self.labile_c = np.argsort (np.mean(self.oxidation_state))
        elif recalcitrance_criterion == 'OC_ratio':
            self.labile_c = np.argsort (np.mean(self.oc_ratio))
        elif recalcitrance_criterion == 'half_saturation_constant':
            self.labile_c = np.argsort(np.mean(self.k_params, axis=1))
        elif recalcitrance_criterion == 'zakem_indicator':
            rho_i_j = self.v_params
            y_i_j = self.y_params
            l_j = self.m_params
            zakem_q_i_j = np.divide(rho_i_j, l_j)*(np.sum(np.multiply(y_i_j,rho_i_j)/rho_i_j) + y_i_j)
            self.labile_c = np.flipud(np.argsort(np.max(zakem_q_i_j, axis=0)))

        self.most_labile_c = self.labile_c[0]
        self.least_labile_c = self.labile_c[-1]
        self.middle_labile_group = self.labile_c[1:]
        
        ### Holding of the identification of fungal and bacterial groups for now
        # identify fungal groups in microbial biomass
        # top microbial species that are able to break down least_labile_c
        # activity_b = np.argsort(self.k_params[self.least_labile_c,:])
        #increasing order of k_params, so increasing order of resistance to use that C compound
        # self.fungi = activity_b[:2]
        # self.bacteria = activity_b[2:]

        print ("Most recalcitrant carbon pool is : ", self.least_labile_c)
        print ("Most labile carbon pool is : ", self.most_labile_c)
        #print ("Fungi pools are :", self.fungi)
        #print ("Bacteria pools are : ", self.bacteria)
    
    def solve_network (self, initial_conditions, time_span, time_space):
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

        def rate_expressions(t, x):
            """Function to set up the reaction network given the 
            number of dissolved organic matter species and microbial
            species.
            """
    
            # assign each ODE to array elements
            C = x[:self.c_n]
            B = x[self.c_n:]

            # define rates for:
            # 1. Production of exoenzymes to depolymerize labile C
            exoenzyme = self.v_enz * B
            
            # 2. Depolymerization of all C by exoenzymes
            c_depoly = self.v_params * C [...,None] * exoenzyme/(self.k_params + exoenzyme)#C[...,None])
            
            # 3. Respiration
            # formula to implement for all B,C combinations:
            b_uptake = self.z * c_depoly
            
            # 4. Growth
            # formula to implement for all B
            # r = y* rate of respiration for each B for all C
            b_growth = self.y_params * b_uptake
            
            # 5. mortality
            # formula to implement for all B
            # r = m*B
            b_mort = self.m_params*B
            
            # 6. dead microbes back to C (distributed evenly across all C)
            
            # formula to implement for all dying B and add it to C
            # r = sum of dying biomass
            b_necromass = np.sum(b_mort)
            
            # 7. Labile DOM species added to less labile DOM species pools
            # Add non-absorbed C fraction (not used for microbial uptake) into recycling fraction
            # In the end we are interested in the bulk unused C product to recycle it to other C pools:
            # r = sum(recycled pool for each C)
            c_recycle = np.sum((1-self.z) * c_depoly, axis = 1)

            # define each ODE
            C_rate = np.empty((self.c_n,0))
            B_rate = np.empty((self.b_n,0))

            # 1. rate of change of carbon species concentration
            # reduction in concentration of carbon due to depolymerization and microbial uptake
            C_rate = -1*np.sum(c_depoly,axis=1)

            # Add bacteria necromass to carbon pools
            if self.necromass_loop == "equally":
                C_rate += b_necromass/self.c_n
            else:
                C_rate[self.middle_labile_group] += b_necromass/(self.middle_labile_group.size)

            # add sequence of addition to the simpler carbon compounds
            for labile_c, less_labile_c in zip(self.labile_c[:-1], self.labile_c[1:]):
                C_rate[labile_c] += c_recycle[less_labile_c]                

            C_rate[self.least_labile_c] += c_recycle[self.least_labile_c]

            # add continuous input to all carbon pools
            C_rate += self.c_bc

            # 8. rate of change of biomass species concentration
            # Biomsdd grows and dies
            B_rate = np.sum(b_growth, axis = 0) - b_mort

            rates = np.append(C_rate, B_rate, axis = 0)

            return rates
        
        ### WIP
        ## Defining boundary conditions/input conditions. This is not working right now.
        #def boundary_conditions(x_bc0, x_bc1):
        #    """Method to set up the boundary conditions given the 
        #    time point and values.
        #    
        #    Parameter
        #    ---------
        #    x_bc0, x_bc1 : Array, float.
        #        x_bc0 : Boundary conditions 1.
        #        x_bc1 : Boundary conditions 2.
        #    """
        #    # Values at t=0:
        #    x_bc0[:self.c_n] = np.ones(self.c_n)*5
        #    x_bc0[self.c_n:] = np.ones(self.b_n)*1
        #
        #    # Values at t=0:
        #    x_bc1[:self.c_n] = np.ones(self.c_n)*10
        #    x_bc1[self.c_n:] = np.ones(self.b_n)*5
        #
        #    # These return values are what we want to be 0:
        #    return [x_bc0, x_bc1]
        ### End of WIP

        # solve
        #(self.initial_guess, d) = odeint (rate_expressions, initial_conditions, time_space, full_output=True)
        self.initial_guess = solve_ivp (rate_expressions, time_span, initial_conditions, t_eval = time_space, rtol = 1e-12, atol = 1e-12)

        print ("Your initial value problem is now solved.")
        #self.answer = solve_bvp(lambda t, x: rate_expressions(t, x, self.c_n, self.b_n),
        #            lambda x_bc0, x_bc1: boundary_conditions(x_bc0, x_bc1, self.c_n, self.b_n), time_space, self.initial_guess.T)

        return self.initial_guess

def generate_random_parameters(dom_n, bio_n, ratio_max_rate_mortality):

    """Function to generate random parameter sets to feed into the reaction network.
        
        Parameter
        ---------
        dom_n : int.
            Number of dissolved organic matter chemical species in the system.
        bio_n : int.
            Number of biomass species in the system.
        ratio_max_rate_mortality: float.
            Ratio to be adopted to modulate maximum rate of dom consumption and mortality
            of microbes to balance microbial die off and consumption of carbon.
        """

    # Initialize the same number of parameters as dom and biomass species:
    enzyme_prod_para = np.random.uniform(0.2, 0.4, bio_n)
    max_rate_para = np.random.normal(0.5, 0.001, bio_n*dom_n)#np.random.uniform(0.001,15,bio_n*dom_n)
    sat_conc_para = np.random.normal(600, 200, bio_n*dom_n)#np.random.randint(1000,8000,bio_n*dom_n)
    efficiency_para = np.random.normal(0.5, 0.01, bio_n)#np.random.uniform(0.1,0.99,bio_n)
    mortality_para = np.min(max_rate_para.reshape(dom_n, bio_n),axis=0)/ratio_max_rate_mortality

    return enzyme_prod_para, max_rate_para, sat_conc_para, efficiency_para, mortality_para

def generate_random_initial_conditions(dom_n, bio_n):

    """Function to generate random initial concentration distribution of carbon and biomass
    species in the system.
        
        Parameter
        ---------
        dom_n : int.
            Number of dissolved organic matter chemical species in the system.
        bio_n : int.
            Number of biomass species in the system.
        random_seed : int.
            Seed for generating 
        """

    # initial conditions
    dom_conc = np.random.normal(1000, 100, dom_n)#np.random.randint(500,1000,dom_n)
    biomass_conc = np.random.normal(60, 15, bio_n)#np.random.randint(100,1000,bio_n)

    return dom_conc, biomass_conc

def generate_random_boundary_conditions(dom_n, user_input, method_name = "uniform_constant"):
    """Function to generate continuous carbon input into the system.
        
        Parameter
        ---------
        method_name : str.

            Choose between different ways to enforce continuous input of carbon:
            uniform_constant : Default.
                            A random value that is added as a continuous source for all carbon pools.
            different_constant : A vector of random values added as a continuous source for the carbon pools.
            user_defined : A vector of values provided by the user as a list.

        """

    # To enforce steady state, identify the simplest carbon compound
        #simplest_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[-1]
        #complex_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[0]
    # Compute fresh carbon input commensurate to the mineralization of the simplest carbon compound
    if method_name == "uniform_constant":
        computed_carbon_input = np.random.normal(1000, 100, dom_n)
    elif method_name == "different_constant":
        computed_carbon_input = np.random.normal(1000, 100, dom_n)
    elif method_name == "user_defined":
        computed_carbon_input = user_input

    return computed_carbon_input
