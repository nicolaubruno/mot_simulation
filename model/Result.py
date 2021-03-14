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
    def code(self):
        return self._code

    #
    # Simulation short name
    @property
    def name(self):
        return self._name

    #
    # 3D-Histogram of the positions
    @property
    def pos_3Dhist(self):
        return self._pos_3Dhist
    
    #
    # Marginal histograms of the positions
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
    def __init__(self, sim_code, loop_idx=0):
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

        if loop_idx < len(self.loop["dirs"]): 
            self._loop["active"] = loop_idx

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

        # 3D-Histograms of the positions
        self._pos_3Dhist = {
            'freqs':[],\
            'dens':[],\
            'bins':[]
        }

        # Marginal histograms of the positions
        self._pos_hist = [{"freqs":[], "dens":[], "bins":[]} for i in range(3)]

        #
        # Frequencies of the 3D-Histograms of the positions
        path = dir_path + 'positions.csv'
        if os.path.exists(path):
            self._pos_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

            # Frequencies
            self._pos_hist[0]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(1, 2))
            self._pos_hist[1]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 2))
            self._pos_hist[2]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 1))

            for i in range(3):
                # Densities
                self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"] / np.sum(self._pos_hist[i]["freqs"])
                self._pos_hist[i]["dens"] = np.array(self._pos_hist[i]["dens"])

                #
                # Bins
                self._pos_hist[i]["bins"] = - np.ones(self.conds['num_bins']) * float(self.conds['r_max'])
                delta = 2*float(self.conds['r_max']) / float(self.conds['num_bins'])

                for j in range(self.conds['num_bins']):
                    self._pos_hist[i]["bins"][j] += j*delta

        #
        # Frequencies of the marginal histograms of the positions
        #

        path = dir_path + 'marginals.csv'
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)
            '''
            for i in range(3):
                self._pos_hist[i]["freqs"] = df[i].to_numpy()
            '''

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
    def loop_idx(self, idx):
        self.__get_attr(idx)

    #
    def mass_centre(self):
        # Variables
        if not self.loop["var"]:
            r_c = [0,0,0]
            std_r_c = [0,0,0]

            for i in range(3):
                x = self.pos_3Dhist["margs"][i]["bins"]
                p = self.pos_3Dhist["margs"][i]["dens"]
                r_c[i] = np.sum(x*p)

            for i in range(3):
                x2 = self.pos_3Dhist["margs"][i]["bins"]*self.pos_3Dhist["margs"][i]["bins"]
                p = self.pos_3Dhist["margs"][i]["dens"]
                std_r_c[i] = np.sqrt(np.sum(x2*p) - r_c[i]**2)

        else:
            r_c = np.zeros((3,len(self.loop["values"])))
            r_c_square = np.zeros((3,len(self.loop["values"])))
            std_r_c = np.zeros((3,len(self.loop["values"])))

            for i in range(len(self.loop["values"])):
                self.loop_idx(i+1)

                for j in range(3):
                    x = self.pos_3Dhist["margs"][j]["bins"]
                    p = self.pos_3Dhist["margs"][j]["dens"]

                    r_c[j][i] = np.sum(x*p)
                    r_c_square[j][i] = np.sum(x*x*p)

            std_r_c = np.sqrt(r_c_square - r_c*r_c)

        return r_c, std_r_c

