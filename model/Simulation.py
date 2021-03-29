#
# Libraries and modules
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext

from multiprocessing import Pool, cpu_count
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
        "_pos_freqs_arr", "_vel_freqs_arr", "_speed_freqs_arr", "_atoms_simulated", "_option", "_results", "_transitions", "_parallel_tasks"
    ]

    #
    # (Array) Position histogram array
    @property
    def pos_freqs_arr(self):
        return self._pos_freqs_arr

    #
    # (Array) Velocity histogram array
    @property
    def vel_freqs_arr(self):
        return self._vel_freqs_arr

    #
    # (Array) Speed histogram array
    @property
    def speed_freqs_arr(self):
        return self._speed_freqs_arr
    
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
    def new(self, shortname = '', opt = 0, results_group = None):
        #
        # Create a directory to save the results
        #

        # Empty results
        if results_group is not None and results_group > 1:
            results_group = self.available_results_groups()[results_group]

        self._results = Results(int(dt.now().timestamp()), shortname.strip(), results_group=results_group)
        
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
            # Frequencies (1D-array)
            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
            self._vel_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
        #--

        #
        # Simulate marginal distributions (Lightweight option)
        #--
        elif self.option == 1:
            #
            # Frequencies of positions (3D-array)
            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
            self._vel_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
        #--

        # Speed distribution
        self._speed_freqs_arr = np.zeros(self.results.conds["num_bins"])

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
                    pos_freqs = res[k][0]
                    vel_freqs = res[k][1]
                    speed_freqs = res[k][2]
                    time = res[k][3]
                    trans = res[k][4]

                    #
                    # Add position and velocity frequencies
                    #--
                    if self.option == 0:
                        for i in range(self.results.conds["num_bins"]**3):
                            self._pos_freqs_arr[i] += pos_freqs[i]
                            self._vel_freqs_arr[i] += vel_freqs[i]

                    elif self.option == 1:                
                        for i in range(3):
                            for j in range(self.results.conds["num_bins"]):
                                self._pos_freqs_arr[i][j] += pos_freqs[i][j]
                                self._vel_freqs_arr[i][j] += vel_freqs[i][j]
                    #--

                    #
                    # Add speed frequencies
                    #--
                    for i in range(self.results.conds["num_bins"]):
                        self._speed_freqs_arr[i] += speed_freqs[i]
                    #--

                #
                # Add transitions
                #--
                for i in range(3):
                    self._transitions[i] += trans[i]
                #--

                #
                # Release memory
                del pos_freqs
                del vel_freqs
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
            # Frequencies (1D-array)
            del self._pos_freqs_arr
            del self._vel_freqs_arr

            self._pos_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
            self._vel_freqs_arr = np.zeros(self.results.conds["num_bins"]**3)
        #--

        #
        # Simulate marginal distributions (Lightweight option)
        #--
        elif self.option == 1:
            #
            # Frequencies (3D-array)
            del self._pos_freqs_arr
            del self._vel_freqs_arr

            self._pos_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
            self._vel_freqs_arr = np.zeros((3, self.results.conds["num_bins"]))
        #--

        # Speed distribution
        self._speed_freqs_arr = np.zeros(self.results.conds["num_bins"])

        # Setup relevant variables
        self._atoms_simulated = 0

        # Parallel tasks
        self._parallel_tasks = int(self.results.conds["parallel_tasks"])

        # Release memory
        gc.collect()

    #
    def save(self):
        #
        # Add position and velocities
        if self.option == 0:
            self.results.add_positions(self.pos_freqs_arr)
            self.results.add_velocities(self.vel_freqs_arr)

        elif self.option == 1:
            self.results.add_marginals(self.pos_freqs_arr, self.vel_freqs_arr)

        #
        # Add speed frequencies
        self.results.add_speeds(self.speed_freqs_arr)

        #
        # Release memory
        gc.collect()

    #
    # Get available result groups
    def available_results_groups(self):
        #
        # Variables
        i = 2
        res = {1:"root"}

        #
        # List all results groups directories
        #--
        groups_dir = os.scandir("model/results/")
        for group_dir in groups_dir:
            # Check if the results group is valid
            if group_dir.is_dir():
                str_splited = group_dir.name.split("_")

                if(str_splited[0] == "group"):
                    name = ""
                    for j in range(1, len(str_splited)):
                        if j == 1: name += str_splited[j]
                        else: name += '_' + str_splited[j]

                    res[i] = name
                    i += 1
        #--

        return res

    #
    # Get available results
    def available_results(self, results_group):
        #
        # Variables
        res = {}
        available_groups = self.available_results_groups()
        check = True

        #
        # Get results dir
        #--
        if results_group == 1:
            path = "model/results/"

        else:
            path = "model/results/group_" + available_groups[results_group]
        #--

        #
        # List all results directories
        #--
        results_dir = os.scandir(path)
        for res_dir_item in results_dir:
            #
            # Check if the results is valid
            if res_dir_item.is_dir() and (len(res_dir_item.name) > 6) and not (res_dir_item.name[:5] == "group"):
                res_dir = os.scandir(res_dir_item.path)
                print(res_dir_item.name)

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

                    res[code] = name
        #--

        return {key: value for key, value in sorted(res.items(), key=(lambda x: -x[0]))}

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