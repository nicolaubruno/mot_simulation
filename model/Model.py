#
# Libraries and modules
import numpy as np
import pandas as pd
import sys

#
class Model:

    ''' Attributes '''

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
    def settings(self):
        return self._settings
    
    #
    @property
    def beams(self):
        return self._beams

    # Magnetic Field Gradient
    @property
    def B_0(self):
        return self._B_0

    #
    # Initial temperature
    @property
    def T_0(self):
        return self._T_0

    #
    # Gravity status
    @property
    def g_bool(self):
        return self._g_bool
    
    #
    # Physical constants
    @property
    def constants(self):
        return self._constants
    

    ''' Methods '''

    #
    def __init__(self):
        self.__load_parameters()        

    #
    def __load_parameters(self):
        #
        # Active parameters
        active_parameters = pd.read_csv('model/parameters/active.csv', header=0, dtype='int').iloc[0]

        #
        # Transition
        path = 'model/parameters/transitions.csv'
        all_transitions = pd.read_csv(path, header=11, index_col='id')
        self. __check_id_in_CSV_file(all_transitions, path)

        self._transition = all_transitions.loc[active_parameters['transition']].astype(object)
        self._transition.index = [idx.strip() for idx in self._transition.index]

        #
        # Atom
        path = 'model/parameters/atoms.csv'
        all_atoms = pd.read_csv(path, header=8, index_col='id')
        self. __check_id_in_CSV_file(all_atoms, path)

        self._atom = all_atoms.loc[self._transition['atom_id']].astype(object)
        self._atom.index = [idx.strip() for idx in self._atom.index]

        #
        # Settings
        path = 'model/parameters/settings.csv'
        all_settings = pd.read_csv(path, header=8, index_col='id', dtype='int')
        self. __check_id_in_CSV_file(all_settings, path)

        self._settings = all_settings.loc[active_parameters['settings']]
        self._settings.index = [idx.strip() for idx in self._settings.index]

        #
        # Initial conditions
        path = 'model/parameters/initial_conditions.csv'
        all_conditions = pd.read_csv(path, header=6, index_col='id')
        self. __check_id_in_CSV_file(all_conditions, path)

        initial_conditions = all_conditions.loc[active_parameters['initial_conditions']]
        initial_conditions.index = [idx.strip() for idx in initial_conditions.index]

        self._T_0 = initial_conditions['T_0'] # Initial Temperature
        self._g_bool = initial_conditions['g_bool'] # Gravity status

        #
        # Environment
        path = 'model/parameters/environments.csv'
        all_environments = pd.read_csv(path, header=7, index_col='id')
        self. __check_id_in_CSV_file(all_environments, path)

        env = all_environments.loc[active_parameters['environment']]
        env.index = [idx.strip() for idx in env.index]

        #
        # Beams
        path = 'model/parameters/beams.csv'
        beams_setup = self.__string_to_array(env['beams_setup'], astype='int')
        all_beams = pd.read_csv(path, header=9, index_col='id')
        self. __check_id_in_CSV_file(all_beams, path)

        self._beams = all_beams.loc[beams_setup]
        self._beams.columns = [col.strip() for col in self._beams.columns]

        for idx in self._beams.index:
            elem = self._beams.at[idx, 'k_dic']
            self._beams.at[idx, 'k_dic'] = self.__string_to_array(elem, astype='float')

            elem = self._beams.at[idx, 'eps']
            self._beams.at[idx, 'eps'] = self.__string_to_array(elem, astype='float')

        #
        # Magnetic field gradient
        self._B_0 = float(env['B_0'])

        #
        # Constants
        path = 'model/parameters/constants.csv'
        self._constants = pd.read_csv(path, header=8, index_col='id')
        self. __check_id_in_CSV_file(self._constants, path)
        self._constants.columns = [idx.strip() for idx in self._constants.columns]

    #
    def run_simulation(self):
        pass       

    #
    # Utility methods
    #

    #
    def __check_id_in_CSV_file(self, dataframe, path):
        if True in dataframe.index.duplicated(keep='first'):
            sys.tracebacklimit = 0
            raise Exception('There are duplicate indexes in CSV file: \n\t%s\nPlease, replace them for unique indexes.' % path)

    #
    def __string_to_array(self, str, astype='int'):
        return np.array(str.replace('[', '').replace(']', '').split(' ')).astype(astype)
