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
    # Positions histogram
    @property
    def pos_hist(self):
        return self._pos_hist


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
        path = dir_path + 'positions.csv'
        self._pos_hist['freqs'] = pd.read_csv(path, index_col=0)[['x', 'y', 'z']].to_numpy()
        self._pos_hist['freqs'] = self._pos_hist['freqs'].T

        #
        # Positions densities
        self._pos_hist['dens'] = []
        
        for i in range(3):
            self._pos_hist['dens'].append(self._pos_hist['freqs'][i] / sum(self._pos_hist['freqs'][i]))

        self._pos_hist['dens'] = np.array(self._pos_hist['dens'])

        #
        # Positions bins
        self._pos_hist["bins"] = np.zeros((3, self.conds['num_bins']))
        self._pos_hist["bins"] = self._pos_hist["bins"] - self.conds['r_max']
        delta = 2*self.conds['r_max'] / self.conds['num_bins']

        for i in range(3):
            for j in range(self.conds['num_bins']):
                self._pos_hist["bins"][i][j] += j*delta

    #
    def show_positions(self):
        pass