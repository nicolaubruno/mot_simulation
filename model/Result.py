#
# Libraries and modules
import sys, os
import numpy as np
import pandas as pd

#
class Result:

    ''' Attributes '''

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
    # (Series) Conditions
    @property
    def conds(self):
        return self._conds

    #
    # (Dataframe)
    @property
    def beams(self):
        return self._beams

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
    # 3D-Histogram of the position distribution
    @property
    def pos_hist(self):
        return self._pos_hist
    
    #
    # Looping status
    @property
    def loop(self):
        return self._loop

    #
    # Results directory
    @property
    def results_dir(self):
        return self._results_dir


    ''' Methods '''

    #
    def __init__(self, sim_code, loop_idx=1):
        #
        # Simulation code
        self._sim_code = sim_code

        #
        # Get simulation name
        self.__get_sim_name()

        #
        # Results directory
        self._results_dir = "model/results/" + str(self.sim_code)
        if len(self.sim_name) > 0: self._results_dir += '_' + self.sim_name
        self._results_dir += '/'

        #
        # Looping status
        self._loop = {
            "var": '',\
            "values": [],\
            "dirs":[],\
            "active":0
        }

        #
        # Get looping status
        self.__get_loop()

        if loop_idx <= (len(self.loop["values"])+1): 
            self._loop["active"] = loop_idx-1

        #
        # Get attributes
        self.__get_attr(self.loop["active"])

    #
    def __get_attr(self, loop_idx):
        #
        # Directory of the result
        dir_path = self.loop["dirs"][loop_idx]

        #
        # Read parameters
        #

        #
        # Atom
        path = dir_path + "parameters/atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Transition
        path = dir_path + 'parameters/transition.csv'
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Beams
        path = dir_path + 'parameters/beams.csv'
        self._beams = pd.read_csv(path, header=0)
        self._beams.index += 1

        #
        # Conditions
        path = dir_path + 'parameters/conditions.csv'
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds['num_sim'] = int(self._conds['num_sim'])
        self._conds['num_bins'] = int(self._conds['num_bins'])

        #
        # Positions histogram
        #

        self._pos_hist = {
            'freqs':[],\
            'dens':[],\
            'bins':[],\
            'margs':[]
        }

        #
        # Positions frequencies
        path = dir_path + 'positions.csv'
        self._pos_hist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

        #
        # Positions densities
        self._pos_hist["margs"] = [{"freqs":[], "dens":[], "bins":[]} for i in range(3)]


        # Frequencies
        self._pos_hist["margs"][0]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(1, 2))
        self._pos_hist["margs"][1]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(0, 2))
        self._pos_hist["margs"][2]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(0, 1))

        for i in range(3):
            # Densities
            self._pos_hist["margs"][i]["dens"] = self._pos_hist["margs"][i]["freqs"] / np.sum(self._pos_hist["margs"][i]["freqs"])
            self._pos_hist["margs"][i]["dens"] = np.array(self._pos_hist["margs"][i]["dens"])

            #
            # Bins
            self._pos_hist["margs"][i]["bins"] = np.zeros(self.conds['num_bins'])
            self._pos_hist["margs"][i]["bins"] = self._pos_hist["margs"][i]["bins"] - self.conds['r_max']
            delta = 2*self.conds['r_max'] / self.conds['num_bins']

            for j in range(self.conds['num_bins']):
                self._pos_hist["margs"][i]["bins"][j] += j*delta

    #
    def __get_loop(self):
        # Variables
        i = 0

        # Scan results directory
        obj_scandir = os.scandir(self.results_dir)

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

    #
    def __get_sim_name(self):
        #
        # Get short name
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        self._sim_name = ''

        for path in obj_scandir:
            str_splited = path.name.split("_")

            code = int(str_splited[0])
            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if code == self.sim_code:
                self._sim_name = name  
                break

    #
    def loop_idx(self, idx):
        self.__get_attr(idx-1)