#
# Libraries and modules
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext

from multiprocessing import cpu_count
from pathos.multiprocessing import ProcessPool as Pool
from datetime import datetime as dt
from model.Result import Result

#
class Simulation:

    ''' Attributes '''

    __slots__ = ["_atom", "_transition", "_conds", "_beams", "_cts", "_code", "_name", "_pos_freqs_arr", "_atoms_simulated","_loop", "_results_dir", "_option"]

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
    # Simulation option
    @property
    def option(self):
        return self._option
    

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
        self._code = None
        self._name = None

        # Loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "active": -1
        }

    #
    def __get_params(self):
        #
        # Parameters directory
        params_dir = "model/parameters/"

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

        #
        # Check delta looping
        values = self.__get_loop_values(str(self.conds["delta"]))
        if len(values) > 0:
            self._loop["var"] = "delta"
            self._loop["values"] = values

        # Release memory
        del params_dir
        del path
        del values

    #
    def new(self, shortname = '', opt = 0):
        #
        # Create a directory to save the results
        #

        # Simulation identification
        self._code = int(dt.now().timestamp())
        self._name = shortname.strip()

        # Check if results directory exists
        self._results_dir = "model/results/"
        if not os.path.exists(self._results_dir):
            os.mkdir(self._results_dir)

        # Create directory
        self._results_dir += str(self.code)
        if self.name: self._results_dir += '_' + self.name
        self._results_dir += '/'
        os.mkdir(self.results_dir)

        #
        # Read parameters
        self.__get_params()

        #
        # Create directories for each result (looping)
        #
        
        # Looping
        num_res = len(self.loop["values"]) if len(self.loop["values"]) > 0 else 1
        for i in range(num_res):
            # Result directory
            if self.loop["var"]:
                res_dir = self.results_dir + "res" + str(i+1) + '_' + self.loop["var"] + '/'

            else:
                res_dir = self.results_dir + "res" + str(i+1) + '/'

            # Create directory
            os.mkdir(res_dir)

            #
            # Save parameters of the simulation
            #

            params_dir = res_dir + "parameters/"
            os.mkdir(params_dir)

            #
            # Add loop variable 
            if self.loop["var"] == "delta":
                self.conds["delta"] = float(self.loop["values"][i])

            self.atom.to_csv(params_dir + "atom.csv")
            self.transition.to_csv(params_dir + "transition.csv")
            self.beams.to_csv(params_dir + "beams.csv", index=False)
            self.conds.to_csv(params_dir + "conditions.csv")

            # Release memory
            del res_dir, params_dir

        # Release memory
        del num_res

        # Set results directory
        self._results_dir = self.results_dir + "res1"
        if self.loop["var"]: self._results_dir += '_' + self.loop["var"]
        self._results_dir += '/'

        #
        # Check simulation option
        # 

        available_opts = {
            0 : "3D distribution",\
            1 : "Marginal distributions"
        }

        if opt in available_opts.keys():
            self._option = opt

        # Release memory
        del available_opts

        #
        # Simulate 3D distribution (Heavy option)
        if self.option == 0:
            #
            # Frequencies of positions (1D-array)
            self._pos_freqs_arr = np.zeros(self.conds["num_bins"]**3)

        #
        # Simulate marginal distributions (Lightweight option)
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            self._pos_freqs_arr = np.zeros((3, self.conds["num_bins"]))

        #
        # Setup relevant variables
        self._atoms_simulated = 0
        if self.loop["var"]: self._loop["active"] = 0

        #
        # Release memory
        gc.collect()

    #
    def run(self, times):
        #
        # Check simulation status
        if self.atoms_simulated < self.conds["num_sim"]:
            #
            # Check number of executions
            if (self.atoms_simulated + times) > self.conds["num_sim"]:
                times = 1

            # 
            # Arguments to the pool 
            args_pool = []
            for i in range(times):
                seed =  int((1000 * (dt.now().timestamp())) % 1000 + 1e3*(self.atoms_simulated + i))
                args_pool.append((self.results_dir + "parameters/", self.option, seed))

            #
            # Wrapper function for pool
            def simulate_atom(args):
                params_dir, opt, s = args
                return C_ext.simulate_atom(params_dir, opt, s)

            #
            # Parallel execution
            with Pool(cpu_count()) as pool:
                freqs = pool.map(simulate_atom, args_pool)

                #
                # Add frequencies
                for k in range(times):
                    if self.option == 0:
                        for i in range(self.conds["num_bins"]**3):
                            self._pos_freqs_arr[i] += freqs[k][i]

                    elif self.option == 1:                
                        for i in range(3):
                            for j in range(self.conds["num_bins"]):
                                self._pos_freqs_arr[i][j] += freqs[k][i][j]

                del freqs
                del args_pool
                del simulate_atom
                del seed

            # Update atoms simulated
            self._atoms_simulated += times

            # Release memory
            gc.collect()

            return times

    #
    def open(self, code, loop_idx = 0, opt = 0):
        # Set loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "active": loop_idx
        }

        # Identification
        self._code = code
        self.__get_name()
        self.__get_loop()

        #
        # Set result directory
        #

        self._results_dir = "model/results/" + str(self.code)
        if self.name: self._results_dir += '_' + self.name
        self._results_dir += '/'

        if self.loop["var"]: 
            self._results_dir += "res" + str(self.loop["active"] + 1) + '_' + self.loop["var"] + '/'
        else:
            self._results_dir += "res1/"

        #
        # Check simulation option
        # 

        available_opts = {
            0 : "3D distribution",\
            1 : "Marginal distributions"
        }

        if opt in available_opts.keys():
            self._option = opt

        # Release memory
        del available_opts

        #
        # Simulate 3D distribution (Heavy option)
        if self.option == 0:
            #
            # Frequencies of positions (1D-array)
            del self._pos_freqs_arr
            self._pos_freqs_arr = np.zeros(self.conds["num_bins"]**3)

        #
        # Simulate marginal distributions (Lightweight option)
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            del self._pos_freqs_arr
            self._pos_freqs_arr = np.zeros((3, self.conds["num_bins"]))

        #
        # Setup relevant variables
        self._atoms_simulated = 0

        # Release memory
        gc.collect()

    #
    def save(self):
        #
        # Check file
        if self.option == 0:
            indexes = []
            values = []

            for i in range(self.conds['num_bins']):
                for j in range(self.conds['num_bins']):
                    for k in range(self.conds['num_bins']):
                        indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                        values.append(self.pos_freqs_arr[self.conds['num_bins']**2 * i + self.conds['num_bins']*j + k])

            values = np.array(values)

            # Save file
            path = self.results_dir + "/positions.csv"
            pos_freqs = pd.Series(values, index=indexes).astype("int32")
            pos_freqs.fillna(0, inplace=True)
            pos_freqs.to_csv(path)

            # Release memory
            del values
            del indexes
            del pos_freqs
            del path

        elif self.option == 1:
            indexes = [i+1 for i in range(self.conds['num_bins'])]
            columns = ["x", "y", "z"]

            data = {
                'x': self._pos_freqs_arr[0],\
                'y': self._pos_freqs_arr[1],\
                'z': self._pos_freqs_arr[2]
            }

            path = self.results_dir + "marginals.csv"
            pos_freqs = pd.DataFrame(data).astype("int32")
            pos_freqs.fillna(0, inplace=True)
            pos_freqs.to_csv(path)

            # Release memory
            indexes.clear()
            columns.clear()

            del indexes
            del columns
            del pos_freqs
            del data
            del path

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

        return res

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

            act_code = int(str_splited[0])
            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if act_code == self.code:
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

        return sorted(values)[::-1]

    #
    def __get_loop(self):
        # Variables
        i = 0

        #
        # Scan results directory
        #

        results_dir = "model/results/" + str(self.code)
        if len(self.name) > 0: results_dir += '_' + self.name
        results_dir += '/'
        obj_scandir = os.scandir(results_dir)

        for i, obj_dir in enumerate(obj_scandir):
            var = obj_dir.name.split("_")

            if i == 0: 
                for j in range(1, len(var)):
                    if j == 1: self._loop["var"] += var[j]
                    else: self._loop["var"] += '_' + var[j]

            #
            # Check looping
            #

            # Delta
            if self.loop["var"] == "delta":
                conds = pd.read_csv(obj_dir.path + "/parameters/conditions.csv" , header=0, index_col=0, squeeze=True).astype(object)
                order = int(var[0][3:])
                self._loop["values"].append((order, float(conds["delta"])))

        if len(self.loop["values"]) > 0:
            self._loop["values"] = np.array(sorted(self.loop["values"], key=(lambda x: x[0])))
            self._loop["values"] = (np.delete(self.loop["values"], 0, 1).T)[0]

