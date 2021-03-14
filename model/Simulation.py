#
# Libraries and modules
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext

from multiprocessing import cpu_count
from pathos.multiprocessing import ProcessPool as Pool
from datetime import datetime as dt
from model.Results import Results

#
class Simulation:

    ''' Attributes '''

    __slots__ = [
        "_pos_freqs_arr", "_atoms_simulated", "_option", "_results"
    ]

    #
    # Simulation parameters
    #

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
    # Simulation option
    @property
    def option(self):
        return self._option
    
    #
    # (Object) Results to be built in the simulation
    @property
    def results(self):
        return self._results
    

    ''' Methods '''

    #
    def __init__(self):
        #
        # Set-up initial values
        self._atoms_simulated = -1
        self._results = None

    #
    def new(self, shortname = '', opt = 0):
        #
        # Create a directory to save the results
        #

        # Empty results
        self._results = Results(int(dt.now().timestamp()), shortname.strip())

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
            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)

        #
        # Simulate marginal distributions (Lightweight option)
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))

        #
        # Setup relevant variables
        self._atoms_simulated = 0

        #
        # Release memory
        gc.collect()

    #
    def run(self, times):
        #
        # Check simulation status
        if self.atoms_simulated < self.results.conds["num_sim"]:
            #
            # Check number of executions
            if (self.atoms_simulated + times) > self.results.conds["num_sim"]:
                times = 1

            # 
            # Arguments to the pool 
            args_pool = []
            for i in range(times):
                seed =  int((1000 * (dt.now().timestamp())) % 1000 + 1e3*(self.atoms_simulated + i))
                args_pool.append((self.results.directory + "parameters/", self.option, seed))

            #
            # Wrapper function for pool
            def simulate_atom(args):
                params_dir, opt, s = args
                return C_ext.simulate_atom(params_dir, opt, s)

            #
            # Parallel execution
            with Pool(cpu_count(), maxtasksperchild=1000) as pool:
                freqs = pool.map(simulate_atom, args_pool)

                #
                # Add frequencies
                for k in range(times):
                    if self.option == 0:
                        for i in range(self.results.conds["num_bins"]**3):
                            self._pos_freqs_arr[i] += freqs[k][i]

                    elif self.option == 1:                
                        for i in range(3):
                            for j in range(self.results.conds["num_bins"]):
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
        #
        # Get results
        self._results = Results(code, loop_idx=loop_idx)

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
            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)

        #
        # Simulate marginal distributions (Lightweight option)
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            del self._pos_freqs_arr
            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))

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

            for i in range(self.results.conds['num_bins']):
                for j in range(self.results.conds['num_bins']):
                    for k in range(self.results.conds['num_bins']):
                        indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                        values.append(self.pos_freqs_arr[self.results.conds['num_bins']**2 * i + self.results.conds['num_bins']*j + k])

            values = np.array(values)

            # Save file
            path = self.results.directory + "/positions.csv"
            pos_freqs = pd.Series(values, index=indexes).astype("int32")
            pos_freqs.fillna(0, inplace=True)
            pos_freqs.to_csv(path)

            # Release memory
            del values
            del indexes
            del pos_freqs
            del path

        elif self.option == 1:
            indexes = [i+1 for i in range(self.results.conds['num_bins'])]
            columns = ["x", "y", "z"]

            data = {
                'x': self._pos_freqs_arr[0],\
                'y': self._pos_freqs_arr[1],\
                'z': self._pos_freqs_arr[2]
            }

            path = self.results.directory + "marginals.csv"
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
    # Check if results code exists
    def check_results_code(self, code):
        #
        # Variables
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")
            check_code = int(str_splited[0])

            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if check_code == int(code):
                ret = True
                break

        return ret  

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
    # Utility methods
    #

