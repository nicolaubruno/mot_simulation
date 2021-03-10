#
# Libraries and modules
import sys, os, gc
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
    # (Series)
    @property
    def atom(self):
        return self._atom

    #
    # (Series)
    @property
    def transition(self):
        return self._transition

    #
    # (Series)
    @property
    def conds(self):
        return self._conds

    #
    # (Dataframe)
    @property
    def beams(self):
        return self._beams

    #
    # Constants
    @property
    def cts(self):
        return self._cts

    #
    # Simulation identification code
    @property
    def sim_code(self):
        return self._sim_code

    #
    # Simulation short name
    @property
    def sim_name(self):
        return self._sim_name

    #
    # (Array) Position histogram array
    @property
    def pos_freqs_arr(self):
        return self._pos_freqs_arr
    
    #
    # Atoms simulated
    @property
    def atoms_simulated(self):
        return self._atoms_simulated
    
    #
    # Iterations
    @property
    def iters(self):
        return self._iters


    ''' Methods '''

    #
    def __init__(self):
        #
        # Constants
        self._cts = {
            'h' :   6.62607004,  # Planck constant [10^{-34} J s]
            'e' :   1.60217662,  # Elementary charge [10^{-19} C]s
            'c' :   2.99792458,  # Speed of light [10^{8} m / s]
            'k_B':  1.38064852,  # Boltzmann constant [10^{-23} J / K]
            'mu_B': 9.274009994, # Bohr magneton [10^{-24} J / T]
            'u':    1.660539040 # Atomic mass unit [10^{-27} kg]
        }

        #
        # Get parameters
        self.__load_parameters()
        self._atoms_simulated = -1
        self._iters = -1
        self._sim_code = -1
        self._sim_name = ''

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
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds['num_sim'] = int(self._conds['num_sim'])
        self._conds['num_bins'] = int(self._conds['num_bins'])

    #
    def start_simulation(self, shortname = ''):
        #
        # Check if results directory exists
        results_dir_path = "model/results/"
        if not os.path.exists(results_dir_path):
            os.mkdir(results_dir_path)

        #
        # Create a directory to store the results of the simulations
        self._sim_code = int(dt.now().timestamp())
        self._sim_name = shortname
        results_dir_path += str(self.sim_code)
        if len(self.sim_name) > 0: results_dir_path += '_' +  self.sim_name + '/'
        os.mkdir(results_dir_path)

        #
        # Save the parameters of the simulation
        parameters_dir = results_dir_path + 'parameters/'
        os.mkdir(parameters_dir)

        self.atom.to_csv(parameters_dir + 'atom.csv')
        self.transition.to_csv(parameters_dir + 'transition.csv')
        self.beams.to_csv(parameters_dir + 'beams.csv')
        self.conds.to_csv(parameters_dir + 'conditions.csv')

        self._pos_freqs_arr = np.zeros((3, self.conds['num_bins']))
        self._atoms_simulated = 0
        self._iters = 0

    #
    def simulate_atom(self):
        #
        # Simulate atoms
        x_bins, y_bins, z_bins, iters, time = C_ext.simulate_atom()

        self._pos_freqs_arr[0] += np.array(x_bins)
        self._pos_freqs_arr[1] += np.array(y_bins)
        self._pos_freqs_arr[2] += np.array(z_bins)

        self._atoms_simulated += 1
        self._iters += iters

        del x_bins
        del y_bins
        del z_bins
        gc.collect()

    #
    def finish_simulation(self):
        #
        # Check file
        path = "model/results/" + str(self.sim_code) 
        if len(self.sim_name) > 0: path += '_' + self.sim_name + '/'
        path += "/positions.csv"

        if os.path.exists(path):
            os.remove(path)

        #
        # Create position histogram data
        columns_name = ["x", "y", "z"]
        indexes = [i+1 for i in range(self.conds['num_bins'])]
        pos_freqs = pd.DataFrame(columns=columns_name, index=indexes)
        pos_freqs.fillna(0, inplace=True)

        pos_freqs["x"] += self.pos_freqs_arr[0]
        pos_freqs["y"] += self.pos_freqs_arr[1]
        pos_freqs["z"] += self.pos_freqs_arr[2]

        pos_freqs["x"] = pos_freqs["x"].astype("int32")
        pos_freqs["y"] = pos_freqs["y"].astype("int32")
        pos_freqs["z"] = pos_freqs["z"].astype("int32")

        pos_freqs.to_csv(path)

    #
    def get_results(self, num=10):
        #
        # Variables
        dir_path = "model/results/"
        res = []
        obj_scandir = os.scandir(dir_path)
        i = 0

        for path in obj_scandir:
            str_splited = path.name.split("_")

            code = int(str_splited[0])
            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            res.append([code, name])
            i += 1

            if i > num - 1:
                break

        return sorted(res, key=lambda x: x[0])[::-1]

    #
    # Check if simulation code exists
    def check_sim_code(self, code):
        #
        # Variables
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")

            sim_code = int(str_splited[0])
            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if int(sim_code) == code:
                ret = True
                break

        return ret  

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
