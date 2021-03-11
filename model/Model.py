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
    # End of simulation
    @property
    def sim_end(self):
        return self._sim_end

    #
    # Looping counter
    @property
    def loop_counter(self):
        return self._loop_counter

    #
    # Looping values
    @property
    def loop_values(self):
        return self._loop_values

    #
    # Looping variable
    @property
    def loop_var(self):
        return self._loop_var

    #
    # Results dirs
    @property
    def results_dir(self):
        return self._results_dir
    

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

        self._atoms_simulated = -1
        self._iters = -1
        self._sim_code = -1
        self._sim_name = ''
        self._pos_freqs_arr = []
        self._sim_end = False
        self._loop_counter = 0
        self._loop_values = []

        #
        # Get parameters
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
        path = "model/parameters/conditions.csv"
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds["num_sim"] = int(self._conds["num_sim"])
        self._conds["num_bins"] = int(self._conds["num_bins"])

        #
        # Check loop variables
        values = self.__get_loop_values(str(self.conds["delta"]))
        if len(values) > 0:
            self._loop_var = "delta"
            self._loop_values = values.copy()

    #
    def start_simulation(self, shortname = ''):
        #
        # Check if results directory exists
        results_dir = "model/results/"
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)

        #
        # Create a directory to store the results of the simulations
        self._sim_code = int(dt.now().timestamp())
        self._sim_name = shortname
        results_dir += str(self.sim_code)
        if len(self.sim_name) > 0: results_dir += '_' +  self.sim_name + '/'
        else: results_dir += str(self.sim_code) + '/'
        os.mkdir(results_dir)
        self._results_dir = results_dir

        #
        # Loop values
        if len(self.loop_values) > 0:
            for i, loop_val in enumerate(self.loop_values):
                loop_dir = results_dir + "loop" + str(i) + '/'
                os.mkdir(loop_dir)

                #
                # Save the parameters of the simulation
                params_dir = loop_dir + "parameters/"
                os.mkdir(params_dir)

                if self.loop_var == "delta":
                    self.conds["delta"] = float(loop_val)

                    self.atom.to_csv(params_dir + "atom.csv")
                    self.transition.to_csv(params_dir + "transition.csv")
                    self.beams.to_csv(params_dir + "beams.csv", index=False)
                    self.conds.to_csv(params_dir + "conditions.csv")

            #
            # Enable simulation
            self._sim_end = False

        else:
            #
            # Save the parameters of the simulation
            loop_dir = results_dir + "loop0/"
            os.mkdir(loop_dir)

            params_dir = loop_dir + "parameters/"
            os.mkdir(params_dir)

            self.atom.to_csv(params_dir + "atom.csv")
            self.transition.to_csv(params_dir + "transition.csv")
            self.beams.to_csv(params_dir + "beams.csv", index=False)
            self.conds.to_csv(params_dir + "conditions.csv")

        self._pos_freqs_arr = np.zeros(self.conds['num_bins']**3)
        self._atoms_simulated = 0

    #
    def simulate(self):
        #
        # Simulate atoms
        if not self.sim_end:
            params_dir = self.results_dir + "loop" + str(self.loop_counter) + "/parameters/"
            freqs = C_ext.simulate_atom(params_dir, int((1000 * (dt.now().timestamp())) % 1000 + 1e3*self.atoms_simulated))

            self._atoms_simulated += 1
            self._pos_freqs_arr += np.array(freqs)

            if self.atoms_simulated == self.conds["num_sim"]:
                self.finish_simulation()

                if self.loop_counter == (len(self.loop_values) - 1):
                    self._sim_end = True

                else:
                    self._loop_counter += 1
                    self._atoms_simulated = 0
                    self._pos_freqs_arr = np.zeros(self.conds['num_bins']**3)

            del freqs
            gc.collect()

    #
    def finish_simulation(self):
        #
        # Check file
        path = self.results_dir + "loop" + str(self.loop_counter) + "/positions.csv"

        if os.path.exists(path):
            os.remove(path)

        #
        # Create position histogram data
        columns_name = ["freq"]
        indexes = []
        values = []

        for i in range(self.conds['num_bins']):
            for j in range(self.conds['num_bins']):
                for k in range(self.conds['num_bins']):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(self.pos_freqs_arr[self.conds['num_bins']**2 * i + self.conds['num_bins']*j + k])

        pos_freqs = pd.Series(np.array(values), index=indexes).astype("int32")
        pos_freqs.fillna(0, inplace=True)

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

    #
    def __get_loop_values(self, loop_str):
        #
        # Return variable
        values = []

        #
        # Check string
        if len(loop_str) > 4 and loop_str[0:4] == 'loop':
            #
            # Loop option 1
            if loop_str[4] == '[' and loop_str[-1] == ']':
                opts = loop_str[5:-1].split(' ')

                val = float(opts[0])
                end = float(opts[1])
                step = float(opts[2])
                
                values = []
                while val <= end:
                    values.append(val)
                    val += step

            elif loop_str[4] == '{' and loop_str[-1] == '}':
                values = loop_str[5:-1].split(' ')

            else:
                raise ValueError('Invalid loop variable')

        return values
