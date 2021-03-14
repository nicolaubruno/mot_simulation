#
# Libraries and modules
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext
from datetime import datetime as dt
from model.Result import Result

#
class Simulation:

    ''' Attributes '''

    __slots__ = ["_atom", "_transition", "_conds", "_beams", "_cts", "_code", "_name", "_pos_freqs_arr", "_atoms_simulated","_loop", "_results_dir", "_only_marginals"]

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
    def code(self):
        return self._code

    #
    # Simulation short name
    @property
    def name(self):
        return self._name

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
    # Looping status
    @property
    def loop(self):
        return self._loop

    #
    # Results dirs
    @property
    def results_dir(self):
        return self._results_dir

    #
    # Only marginal histograms
    @property
    def only_marginals(self):
        return self._only_marginals
    

    ''' Methods '''

    #
    def __init__(self, code=None, loop_idx=0, only_marginals=0):
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

        if only_marginals == 0:
            self._only_marginals = 0

        elif only_marginals == 1:
            self._only_marginals = 1

        self._atoms_simulated = 0
        self._code = code
        self._name = None

        # Loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "dirs": [],
            "active": 0
        }

        #
        # Load existing simulation
        if  self.code is not None:
            self.__get_name()
            self.__get_loop()

            if loop_idx < len(self.loop["dirs"]): 
                self._loop["active"] = loop_idx

        #
        # Get parameters
        self.__get_attr()

    #
    # Atoms simulated
    @atoms_simulated.setter
    def atoms_simulated(self, num):
        self._atoms_simulated = num

    #
    def __get_attr(self):
        #
        # Parameters directory
        if self.code is not None: params_dir = self.loop["dirs"][self.loop["active"]] + "parameters/"
        else: params_dir = "model/parameters/"

        #
        # Atom
        path = params_dir + 'atom.csv'
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Transition
        path = params_dir + 'transition.csv'
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Beams
        path = params_dir + 'beams.csv'
        self._beams = pd.read_csv(path, header=0)
        self._beams.index += 1

        #
        # Conditions
        path = params_dir + "conditions.csv"
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds["num_sim"] = int(self._conds["num_sim"])
        self._conds["num_bins"] = int(self._conds["num_bins"])

        if self.only_marginals == 0:
            self._pos_freqs_arr = np.zeros(self.conds["num_bins"]**3)

        elif self.only_marginals == 1:
            self._pos_freqs_arr = np.zeros((3, self.conds["num_bins"]))

        #
        # Check delta looping
        if self.code is None:
            values = self.__get_loop_values(str(self.conds["delta"]))

            if len(values) > 0:
                self._loop["var"] = "delta"
                self._loop["values"] = values.copy()

    #
    def start(self, shortname = '', only_marginals=0):
        #
        self._only_marginals = only_marginals

        if self.only_marginals == 0:
            self._pos_freqs_arr = np.zeros(self.conds["num_bins"]**3)

        elif self.only_marginals == 1:
            self._pos_freqs_arr = np.zeros((3, self.conds["num_bins"]))            

        #
        # Check simulation option
        if self.code is None:
            #
            # Check if results directory exists
            results_dir = "model/results/"
            if not os.path.exists(results_dir):
                os.mkdir(results_dir)

            #
            # Create a directory to store the results of the simulations
            self._code = int(dt.now().timestamp())
            self._name = shortname
            results_dir += str(self.code)
            if len(self.name) > 0: results_dir += '_' +  self.name + '/'
            else: results_dir += str(self.code) + '/'
            os.mkdir(results_dir)
            self._results_dir = results_dir

            #
            # Number of atoms simulated
            self._atoms_simulated = 0

            #
            # Loop values
            if len(self.loop["values"]) > 0:
                for i, loop_val in enumerate(self.loop["values"]):
                    #
                    # Create looping directory
                    self._loop["dirs"].append(results_dir + "loop" + str(i+1) + '_' + self.loop["var"] + '/')
                    os.mkdir(self.loop["dirs"][-1])

                    #
                    # Save the parameters of the simulation
                    params_dir = self.loop["dirs"][-1] + "parameters/"
                    os.mkdir(params_dir)

                    #
                    # Add loop variable 
                    if self.loop["var"] == "delta":
                        self.conds["delta"] = float(loop_val)

                    self.atom.to_csv(params_dir + "atom.csv")
                    self.transition.to_csv(params_dir + "transition.csv")
                    self.beams.to_csv(params_dir + "beams.csv", index=False)
                    self.conds.to_csv(params_dir + "conditions.csv")

            else:
                #
                # Save the parameters of the simulation
                self._loop["dirs"].append(results_dir + "loop1/")
                os.mkdir(self.loop["dirs"][-1])

                params_dir = self.loop["dirs"][-1] + "parameters/"
                os.mkdir(params_dir)

                self.atom.to_csv(params_dir + "atom.csv")
                self.transition.to_csv(params_dir + "transition.csv")
                self.beams.to_csv(params_dir + "beams.csv", index=False)
                self.conds.to_csv(params_dir + "conditions.csv")

    #
    def simulate(self, only_marginals):
        #
        # Simulate atoms
        params_dir = self.loop["dirs"][self.loop["active"]] + "parameters/"
        seed =  int((1000 * (dt.now().timestamp())) % 1000 + 1e3*self.atoms_simulated)
        freqs = C_ext.simulate_atom(params_dir, only_marginals, seed)

        return np.array(freqs)

    #
    def save(self):
        #
        # Saved results
        res = Result(self.code, self.loop["active"])

        #
        # Check file
        if self.only_marginals == 0:
            indexes = []
            values = []

            for i in range(self.conds['num_bins']):
                for j in range(self.conds['num_bins']):
                    for k in range(self.conds['num_bins']):
                        indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                        values.append(self.pos_freqs_arr[self.conds['num_bins']**2 * i + self.conds['num_bins']*j + k])

            values = np.array(values)

            #
            # Save file
            path = self.loop["dirs"][self.loop["active"]] + "/positions.csv"
            pos_freqs = pd.Series(values, index=indexes).astype("int32")
            pos_freqs.fillna(0, inplace=True)
            pos_freqs.to_csv(path)

        elif self.only_marginals == 1:
            indexes = [i+1 for i in range(res.conds['num_bins'])]
            columns = ["x", "y", "z"]

            data = {
                'x': self._pos_freqs_arr[0],\
                'y': self._pos_freqs_arr[1],\
                'z': self._pos_freqs_arr[2]
            }

            path = self.loop["dirs"][self.loop["active"]] + "marginals.csv"
            pos_freqs = pd.DataFrame(data).astype("int32")
            pos_freqs.fillna(0, inplace=True)
            pos_freqs.to_csv(path)

        #
        # Release memory
        gc.collect()

    #
    def get_results(self, num=10):
        #
        # Variables
        res = []
        obj_scandir = os.scandir("model/results/")
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
    def check_code(self, code):
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

            if sim_code == int(code):
                ret = True
                break

        return ret  

    #
    def __get_name(self):
        #
        # Get short name
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        self._name = ''

        for path in obj_scandir:
            str_splited = path.name.split("_")

            code = int(str_splited[0])
            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if code == self.code:
                self._name = name  
                break

    #
    def update_pos_freqs(self, freqs):
        self._pos_freqs_arr += np.array(freqs)

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

    #
    def __get_loop(self):
        # Variables
        i = 0

        # Scan results directory
        results_dir = "model/results/" + str(self.code)
        if len(self.name) > 0: results_dir += '_' + self.name
        results_dir += '/'
        obj_scandir = os.scandir(results_dir)

        for obj_dir in obj_scandir:
            if i == 0: 
                var = obj_dir.name.split("_")

                for j in range(1, len(var)):
                    if j == 1: self._loop["var"] += var[j]
                    else: self._loop["var"] += '_' + var[j]

            self._loop["dirs"].append(obj_dir.path + '/')

            if self.loop["var"] == "delta":
                conds = pd.read_csv(self.loop["dirs"][-1] + "parameters/conditions.csv", header=0, index_col=0, squeeze=True).astype(object)
                self._loop["values"].append(float(conds["delta"]))

            i += 1

