#
# Libraries and modules
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext

from multiprocessing import cpu_count
from multiprocessing import Pool
from datetime import datetime as dt
from model.Results import Results
from random import randint

#
# Wrapper function for pool
def simulate_atom(args):
    params_dir, opt, s = args
    return C_ext.simulate_atom(params_dir, opt, s)

#
class Simulation:

    ''' Attributes '''

    __slots__ = [
        "_pos_freqs_arr", "_atoms_simulated", "_option", "_results", "_transitions", "_parallel_tasks"
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
    # Occurred transitions
    @property
    def transitions(self):
        return self._transitions

    #
    # Simulation option
    @property
    def option(self):
        return self._option
    
    #
    # Parallel tasks
    @property
    def parallel_tasks(self):
        return self._parallel_tasks
    

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
        self._transitions = np.zeros(3);
        self._results = None
        self._parallel_tasks = 0

    #
    def new(self, shortname = '', opt = 0):
        #
        # Create a directory to save the results
        #

        # Empty results
        self._results = Results(int(dt.now().timestamp()), shortname.strip())
        
        #
        # Check simulation option
        #--
        available_opts = {
            0 : "3D distribution",\
            1 : "Marginal distributions"
        }

        if opt in available_opts.keys():
            self._option = opt

        # Release memory
        del available_opts
        #--

        #
        # Simulate 3D distribution (Heavy option)
        #--
        if self.option == 0:
            #
            # Frequencies of positions (1D-array)
            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
        #--

        #
        # Simulate marginal distributions (Lightweight option)
        #--
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
        #--

        # Setup relevant variables
        self._atoms_simulated = 0

        # Transitions
        self._transitions = np.zeros(3)

        # Parallel tasks
        self._parallel_tasks = int(self.results.conds["parallel_tasks"])

        #
        # Release memory
        gc.collect()

    #
    def run(self):
        #
        # Check simulation status
        if self.atoms_simulated < self.results.conds["num_sim"]:
            #
            # Check number of executions
            if (self.atoms_simulated + self.parallel_tasks) > self.results.conds["num_sim"]:
                times = 1
            else:
                times = self.parallel_tasks

            # 
            # Arguments to the pool 
            args_pool = []
            for i in range(times):
                seed = int(randint(0, 1e10))
                args_pool.append((self.results.directory + "parameters/", self.option, seed))

            #
            # Parallel execution
            with Pool(cpu_count(), maxtasksperchild=100) as pool:
                res = pool.map(simulate_atom, args_pool)

                #
                # Add frequencies
                for k in range(times):
                    freqs = res[k][0]
                    time = res[k][1]
                    trans = res[k][2]

                    if self.option == 0:
                        for i in range(self.results.conds["num_bins"]**3):
                            self._pos_freqs_arr[i] += freqs[i]

                    elif self.option == 1:                
                        for i in range(3):
                            for j in range(self.results.conds["num_bins"]):
                                self._pos_freqs_arr[i][j] += freqs[i][j]

                #
                # Add transitions
                #--
                for i in range(3):
                    self._transitions[i] += trans[i]
                #--

                #
                # Release memory
                del freqs
                del args_pool
                del seed

                #
                # Finish pool and release memory
                pool.terminate()

            #
            # Update atoms simulated
            self._atoms_simulated += times

            # Release memory
            gc.collect()

            return times

    #
    def open(self, code, loop_idx = 0, opt = 0):
        # Get results
        self._results = Results(code, loop_idx=loop_idx)

        #
        # Check simulation option
        # --
        available_opts = {
            0 : "3D distribution",\
            1 : "Marginal distributions"
        }

        if opt in available_opts.keys():
            self._option = opt

        # Release memory
        del available_opts
        #--

        #
        # Simulate 3D distribution (Heavy option)
        #--
        if self.option == 0:
            #
            # Frequencies of positions (1D-array)
            del self._pos_freqs_arr
            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
        #--

        #
        # Simulate marginal distributions (Lightweight option)
        #--
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            del self._pos_freqs_arr
            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
        #--

        # Setup relevant variables
        self._atoms_simulated = 0

        # Parallel tasks
        self._paralell_tasks = int(self.results.conds["parallel_tasks"])

        # Release memory
        gc.collect()

    #
    def save(self):
        #
        # Check file
        if self.option == 0:
            self.results.add_positions(self.pos_freqs_arr)

        elif self.option == 1:
            self.results.add_marginal_positions(self.pos_freqs_arr)

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
    def available_results(self, max_num_results=10):
        #
        # Variables
        res = []
        check = True

        #
        # List all results directories
        results_dir = os.scandir("model/results/")
        for res_dir_item in results_dir:
            #
            # Check if the results is valid
            if res_dir_item.is_dir():
                res_dir = os.scandir(res_dir_item.path)

                check = True
                for res_item in res_dir:
                    if not check: break
                    elif res_item.is_dir():
                        res_loop = os.scandir(res_item.path)
                        check  = False

                        for res_loop_item in res_loop:
                            if res_loop_item.name == "marginals.csv":
                                check = True
                                break

                if check:
                    str_splited = res_dir_item.name.split("_")

                    code = int(str_splited[0])
                    name = ""
                    for j in range(1, len(str_splited)):
                        if j == 1: name += str_splited[j]
                        else: name += '_' + str_splited[j]

                    res.append([code, name])

        return sorted(res, key=lambda res: -res[0])[:max_num_results]

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
            sim_code = int(str_splited[0])

            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if sim_code == int(code):
                ret = True
                break

        return ret  