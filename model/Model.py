#
# Libraries and modules
import sys, os, hashlib
import numpy as np
import pandas as pd
from datetime import datetime as dt
import mot_sim as C_ext

#
class Model:

    ''' Attributes '''

    #
    # Simulation parameters
    #

    #
    # Series
    @property
    def atom(self):
        return self._atom
    
    #
    # Series
    @property
    def transition(self):
        return self._transition

    #
    # Series
    @property
    def conditions(self):
        return self._conditions

    #
    # Dataframe
    @property
    def beams(self):
        return self._beams

    #
    # Dataframe
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
        # Atom
        path = 'model/parameters/atom.csv'
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Transition
        path = 'model/parameters/transition.csv'
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Beams
        path = 'model/parameters/beams.csv'
        self._beams = pd.read_csv(path, header=0)
        self._beams.index += 1

        #
        # Conditions
        path = 'model/parameters/conditions.csv'
        self._conditions = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Constants
        path = 'model/parameters/constants.csv'
        self._constants = pd.read_csv(path, header=0, index_col=0)

    #
    def run_simulation(self):
        #
        # Check if results directory exists
        results_dir_path = "model/results/"
        if not os.path.exists(results_dir_path):
            os.mkdir(results_dir_path)

        #
        # Create a directory to store the results of the simulations
        code = str(int(dt.now().timestamp()))
        results_dir_path += code + '/'
        os.mkdir(results_dir_path)

        #
        # Save the parameters of the simulation
        parameters_dir = results_dir_path + 'parameters/'
        os.mkdir(parameters_dir)

        self.atom.to_csv(parameters_dir + 'atom.csv')
        self.transition.to_csv(parameters_dir + 'transition.csv')
        self.beams.to_csv(parameters_dir + 'beams.csv')
        self.conditions.to_csv(parameters_dir + 'conditions.csv')
        self.constants.to_csv(parameters_dir + 'constants.csv')

        C_ext.simulate_atom(code)

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
