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
    
    #
    @property
    def beams(self):
        return self._beams

    #
    @property
    def magnetic_field(self):
        return self._magnetic_field
    

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
        all_transitions = pd.read_csv('model/parameters/transitions.csv', header=12, index_col='id')
        self._transition = all_transitions.loc[active_parameters['transition']].astype(object)
        self._transition.index = [idx.strip() for idx in self._transition.index]

        #
        # Atom
        #
        all_atoms = pd.read_csv('model/parameters/atoms.csv', header=8, index_col='id')
        self._atom = all_atoms.loc[self._transition['atom_id']].astype(object)
        self._atom.index = [idx.strip() for idx in self._atom.index]

        #
        # Conditions
        #
        all_conditions = pd.read_csv('model/parameters/conditions.csv', header=10, index_col='id')
        self._conditions = all_conditions.loc[active_parameters['conditions']]
        self._conditions.index = [idx.strip() for idx in self._conditions.index]

        #
        # Environment
        #
        all_environments = pd.read_csv('model/parameters/environments.csv', header=8, index_col='id')
        env = all_environments.loc[active_parameters['environment']]
        env.index = [idx.strip() for idx in env.index]

        #
        # Beams
        #
        beams_setup = self.__string_to_array(env['beams_setup'], astype='int')
        all_beams = pd.read_csv('model/parameters/beams.csv', header=9, index_col='id')
        self._beams = all_beams.loc[beams_setup]
        self._beams.columns = [col.strip() for col in self._beams.columns]

        for idx in self._beams.index:
            elem = self._beams.at[idx, 'k_dic']
            self._beams.at[idx, 'k_dic'] = self.__string_to_array(elem, astype='float')

            elem = self._beams.at[idx, 'eps']
            self._beams.at[idx, 'eps'] = self.__string_to_array(elem, astype='float')

        #
        # Magnetic field
        #
        self._magnetic_field = env[['B_0', 'B_direction']]
        self._magnetic_field['B_0'] = float(self._magnetic_field['B_0'])
        self._magnetic_field['B_direction'] = self.__string_to_array(self._magnetic_field['B_direction'], astype='float')

    #
    def __string_to_array(self, str, astype='int'):
        return np.array(str.replace('[', '').replace(']', '').split(' ')).astype(astype)
