#
# Libraries and modules
import sys, os
import numpy as np
import pandas as pd

#
class Results:

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
    # Identification code
    @property
    def code(self):
        return self._code

    #
    # Identification short name
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
    def directory(self):
        return self._directory


    ''' Methods '''

    #
    def __init__(self, code, name = '', loop_idx = 0):
        # Set loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "active": int(loop_idx)
        }

        #
        # Identification
        #

        self._code = None
        self._name = name.strip()

        # Get existent results
        if self.__check_code(code):
            # Get code
            self._code = code

            # Get name
            if len(self.name) > 0:
                if not self.__check_name(code):
                    raise ValueError('Name is not exists')

            else:
                self.__get_name()

            # Get loop
            self.__get_loop()

            # Get parameters
            self.__get_attr()

            # Get distributions
            self.__get_dists()

        # Create a new results
        else: self.__new(code, name)

    #
    # Get attributes
    def __get_attr(self):
        #
        # Change directory
        #

        self._directory = "model/results/" + str(self.code)

        if self.name:  
            self._directory += '_' + self._name

        self._directory += '/'

        if len(self.loop["var"]) > 0:
            self._directory += "res" + str(self.loop["active"] + 1) \
                + '_' + self.loop["var"] + '/'
        else:
            self._directory += "res1/"

        # Parameters directory
        params_dir =  self._directory + "parameters/"

        #
        # Read parameters
        #

        #
        # Atom
        path = params_dir + "atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Transition
        path = params_dir + "transition.csv"
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Beams
        path = params_dir + "beams.csv"
        self._beams = pd.read_csv(path, header=0)
        self._beams.index += 1

        #
        # Conditions
        path = params_dir + "conditions.csv"
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds['num_sim'] = int(self._conds['num_sim'])
        self._conds['num_bins'] = int(self._conds['num_bins'])

    #
    # Get distributions
    def __get_dists(self):
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
        path = self.directory + 'positions.csv'
        if os.path.exists(path):
            self._pos_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

            # Frequencies
            self._pos_hist[0]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(1, 2))
            self._pos_hist[1]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 2))
            self._pos_hist[2]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 1))

        #
        # Frequencies of the marginal histograms of the positions
        path = self.directory + 'marginals.csv'
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)

            # Frequencies
            self._pos_hist[0]["freqs"] = df['x'].to_numpy()
            self._pos_hist[1]["freqs"] = df['y'].to_numpy()
            self._pos_hist[2]["freqs"] = df['z'].to_numpy()

        #
        # Densities and bins
        if len(self._pos_hist[0]["freqs"]) > 0:
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
        #
        # Get root directory            
        root_dir = "model/results/" + str(self.code)

        if self.name:  
            root_dir += '_' + self._name

        root_dir += '/'

        # Variables
        i = 0

        # Scan results directory
        obj_scandir = os.scandir(root_dir)

        for obj_dir in obj_scandir:
            if i == 0: 
                var = obj_dir.name.split("_")

                for j in range(1, len(var)):
                    if j == 1: self._loop["var"] += var[j]
                    else: self._loop["var"] += '_' + var[j]

            if self.loop["var"] == "delta":
                conds = pd.read_csv(obj_dir.path + "/parameters/conditions.csv", header=0, index_col=0, squeeze=True).astype(object)
                self._loop["values"].append(float(conds["delta"]))

            i += 1

        self.loop["values"] = sorted(self.loop["values"])[::-1]

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
    def __new(self, code, name):
        self._code = code
        self._name = name

        # Check if results directory exists
        self._directory = "model/results/"
        if not os.path.exists(self._directory):
            os.mkdir(self._directory)

        # Create directory
        self._directory += str(self.code)
        if self.name: self._directory += '_' + self.name
        self._directory += '/'
        os.mkdir(self.directory)

        #
        # Create directories for each result (looping)
        #

        # Create new attributes
        self.__create_attr()
        
        # Looping
        num_res = len(self.loop["values"]) if len(self.loop["values"]) > 0 else 1
        for i in range(num_res):
            # Result directory
            if len(self.loop["var"]) > 0:
                res_dir = self.directory + "res" + str(i+1) + '_' + self.loop["var"] + '/'

            else:
                res_dir = self.directory + "res1/"

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

        # Set results directory
        self._directory = self.directory + "res"

        if len(self.loop["var"]) > 0: 
            self._directory += str(self.loop['active'] + 1) + '_' + self.loop["var"]

        else:
            self._directory += '1'

        self._directory += '/'

        # Release memory
        del num_res

    #
    # Create attributes
    def __create_attr(self):
        # Parameters directory
        params_dir = "model/parameters/"

        #
        # Atom
        path = params_dir + "atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Transition
        path = params_dir + "transition.csv"
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Beams
        path = params_dir + "beams.csv"
        self._beams = pd.read_csv(path, header=0)
        self._beams.index += 1

        #
        # Conditions
        path = params_dir + "conditions.csv"
        self._conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._conds['num_sim'] = int(self._conds['num_sim'])
        self._conds['num_bins'] = int(self._conds['num_bins'])

        #
        # Check delta looping
        values = self.__get_loop_values(str(self.conds["delta"]))
        if len(values) > 0:
            self._loop["var"] = "delta"
            self._loop["values"] = values
            self.conds["delta"] = self.loop["values"][self.loop["active"]]

    #
    def mass_centre(self):
        # Variables
        if len(self.loop["var"]) == 0:
            r_c = [0,0,0]
            std_r_c = [0,0,0]

            for i in range(3):
                x = self.pos_hist[i]["bins"]
                p = self.pos_hist[i]["dens"]
                r_c[i] = np.sum(x*p)

            for i in range(3):
                x2 = self.pos_hist[i]["bins"]*self.pos_hist[i]["bins"]
                p = self.pos_hist[i]["dens"]
                std_r_c[i] = np.sqrt(np.sum(x2*p) - r_c[i]**2)

        else:
            r_c = np.zeros((3,len(self.loop["values"])))
            r_c_square = np.zeros((3,len(self.loop["values"])))
            std_r_c = np.zeros((3,len(self.loop["values"])))

            for i in range(len(self.loop["values"])):
                self.loop_idx(i)

                for j in range(3):
                    x = self.pos_hist[j]["bins"]
                    p = self.pos_hist[j]["dens"]

                    r_c[j][i] = np.sum(x*p)
                    r_c_square[j][i] = np.sum(x*x*p)

            std_r_c = np.sqrt(r_c_square - r_c*r_c)

        return r_c, std_r_c

    #
    # Check if code exists
    def __check_code(self, code):
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
    # Check if name exists
    def __check_name(self, name):
        #
        # Variables
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")
            code = int(str_splited[0])

            check_name = ""
            for j in range(1, len(str_splited)):
                if j == 1: check_name += str_splited[j]
                else: check_name += '_' + str_splited[j]

            if code == self.code:
                if check_name == name:
                    ret = True
                    break

        return ret  

    #
    def loop_idx(self, idx):
        self._loop["active"] = idx
        self.__get_attr()
        self.__get_dists()
