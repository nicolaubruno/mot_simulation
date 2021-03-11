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
    # Histogram of the positions marginals distributions
    @property
    def position_marginal_histogram(self):
        return self._position_marginal_histogram
    

    ''' Methods '''

    #
    def __init__(self, sim_code):
        #
        # Simulation code
        self._sim_code = sim_code

        #
        # Get attributes
        self.__get_attr()

    #
    def __get_attr(self):
        #
        # Get short name
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)

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
        # Directory of the result
        dir_path += str(self.sim_code) + '_' + self.sim_name + '/'

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
        }

        #
        # Positions frequencies
        path = dir_path + 'frequencies_positions.csv'
        self._pos_hist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

        #
        # Positions densities
        self._position_marginal_histogram = [{"freqs":[], "dens":[], "bins":[]} for i in range(3)]


        # Frequencies
        self._position_marginal_histogram[0]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(1, 2))
        self._position_marginal_histogram[1]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(0, 2))
        self._position_marginal_histogram[2]["freqs"] = np.sum(self.pos_hist["freqs"], axis=(0, 1))

        for i in range(3):
            # Densities
            self._position_marginal_histogram[i]["dens"] = self._position_marginal_histogram[i]["freqs"] / np.sum(self._position_marginal_histogram[i]["freqs"])
            self._position_marginal_histogram[i]["dens"] = np.array(self._position_marginal_histogram[i]["dens"])

            #
            # Positions bins
            self._position_marginal_histogram[i]["bins"] = np.zeros(self.conds['num_bins'])
            self._position_marginal_histogram[i]["bins"] = self._position_marginal_histogram[i]["bins"] - self.conds['r_max']
            delta = 2*self.conds['r_max'] / self.conds['num_bins']

            for j in range(self.conds['num_bins']):
                self._position_marginal_histogram[i]["bins"][j] += j*delta

    #
    def show_positions(self):
        pass