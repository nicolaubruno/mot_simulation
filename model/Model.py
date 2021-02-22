#
# Libraries and modules
import numpy as np
import pandas as pd

#
class Model:

    ''' Attributes '''

    #
    # Constant
    #

    # Planck's constant
    @property
    def h(self):
        return 6.62607004e-34 # J * s

    # Elementary charge
    @property
    def e(self):
        return 1.60217662e-19 # C

    # Speed of light
    @property
    def c(self):
        return 2.99792458e8 # m / s

    # Boltzmann constant
    @property
    def k_B(self):
        return 1.38064852e-23 # J / K

    # Bohr magneton
    @property
    def mu_B(self):
        return 9.274009994e-24 # J / T
    
    #
    # Simulation parameters
    #

    #
    @property
    def atom(self):
        return self._atom
    
    #
    @property
    def transition(self):
        return self._transition

    #
    @property
    def conditions(self):
        return self._conditions
    
    

    ''' Methods '''

    #
    def __init__(self):
        self.__load_parameters()        

    #
    def __load_parameters(self):
        #
        # Active parameters
        #
        active_parameters = pd.read_csv('model/parameters/active.csv', header=0).iloc[0]

        #
        # Transition
        #
        all_transitions = pd.read_csv('model/parameters/transitions.csv', header=12, index_col=0)
        transition = all_transitions.loc[active_parameters['transition']].astype(object)
        self._transition = transition.iloc[1:]

        #
        # Atom
        #
        all_atoms = pd.read_csv('model/parameters/atoms.csv', header=8, index_col=0)
        self._atom = all_atoms.loc[transition['atom_id']].astype(object)

        #
        # Conditions
        #
        all_conditions = pd.read_csv('model/parameters/conditions.csv', header=10, index_col=0)
        self._condition = all_conditions.loc[active_parameters['conditions']]

        #
        # Environment
        #