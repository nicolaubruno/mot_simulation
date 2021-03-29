#
# Libraries and modules
#--
import sys, os, gc
import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
#--

#
class Results:
    ''' Attributes '''

    #
    # Constants in the international system
    @property
    def ctes(self):
        return self._ctes
 
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
    # (Series) Environment
    @property
    def env(self):
        return self._env

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
    # 3D-Histogram of positions
    @property
    def pos_3Dhist(self):
        return self._pos_3Dhist
    
    #
    # Marginal histograms of positions
    @property
    def pos_hist(self):
        return self._pos_hist

    #
    # 3D-Histogram of velocities
    @property
    def vel_3Dhist(self):
        return self._vel_3Dhist
    
    #
    # Marginal histograms of velocities
    @property
    def vel_hist(self):
        return self._vel_hist
    
    #
    # Speed histogram
    @property
    def speed_hist(self):
        return self._speed_hist

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

    #
    # Root directory
    @property
    def root_dir(self):
        return self._root_dir
    

    ''' Methods '''

    #
    def __init__(self, code, name = '', loop_idx = 0, results_group = None):
        #
        # Constants in the international system
        #--
        self._ctes = {
            'u': 1.660539040e-27,\
            'k_B': 1.38064852e-23,\
            'hbar': 1.0544718e-34,\
            'h': 6.626070040e-34
        }
        #--

        #
        # Set loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "active": int(loop_idx)
        }

        #
        # Root dir
        self._root_dir = "model/results/"
        if results_group is not None and results_group != 1: 
            self._root_dir += "group_" + results_group + "/"

        #
        # Identification
        self._code = None
        self._name = name.strip()

        #
        # Get existent results
        #--
        if self.__check_code(code):
            # Get code
            self._code = code

            # Get name
            if len(self.name) > 0:
                if not self.__check_name(code):
                    raise ValueError('Name is not exists')

            else:
                self.__get_name()

            # Get parameters
            self.__get_attr()

            # Get loop
            self.__get_loop()

            # Check loop
            if len(self.loop["var"]) > 0:
                self.__get_attr()

            # Get distributions
            self.__get_dists()

        # Create a new results
        else: self.__new(code, name)
        #--

        #
        # Cast values of the parameters
        self.__cast_params_values()

    #
    # Cast values of the parameters
    def __cast_params_values(self):
        # Atom
        self._atom['Z'] = int(self.atom['Z'])
        self._atom['mass'] = float(self.atom['mass'])

        # Transition
        self._transition['gamma'] = float(self.transition['gamma'])
        self._transition['lambda'] = float(self.transition['lambda'])
        self._transition['g_gnd'] = float(self.transition['g_gnd'])
        self._transition['g_exc'] = float(self.transition['g_exc'])
        self._transition['J_gnd'] = int(self.transition['J_gnd'])
        self._transition['J_exc'] = int(self.transition['J_exc'])

        # Environment
        self._env['B_0'] = float(self.env['B_0'])
        self._env['local_B'] = float(self.env['local_B'])
        self._env['delta'] = float(self.env['delta'])
        self._env['s_0'] = float(self.env['s_0'])
        self._env['w'] = float(self.env['w'])
        self._env['g_bool'] = int(self.env['g_bool'])

        # Conditions
        self._conds['T_0'] = float(self.conds['T_0'])
        self._conds['max_time'] = float(self.conds['max_time'])
        self._conds['wait_time'] = float(self.conds['wait_time'])
        self._conds['dt'] = float(self.conds['dt'])
        self._conds['max_r'] = float(self.conds['max_r'])
        self._conds['max_v'] = float(self.conds['max_v'])
        self._conds['num_sim'] = int(self.conds['num_sim'])
        self._conds['num_bins'] = int(self.conds['num_bins'])
        self._conds['parallel_tasks'] = int(self.conds['parallel_tasks'])

        #
        # Cast loop
        #--
        int_params = ['J_gnd', 'J_exc', 'g_bool', 'num_sim', 'num_bins', 'parallel_tasks']

        if len(self.loop["var"]) > 0:
            # Integer values
            if self.loop["var"] in int_params:
                for i, val in enumerate(self.loop["values"]):
                    self.loop["values"][i] = int(val)

            # Float values
            else:
                for i, val in enumerate(self.loop["values"]):
                    self.loop["values"][i] = float(val)
        #--

    #
    # Get attributes
    def __get_attr(self):
        #
        # Change directory
        #--
        self._directory = self.root_dir + str(self.code)

        if self.name:  
            self._directory += '_' + self._name

        self._directory += '/'

        if len(self.loop['var']) == 0:
            res_dir = os.scandir(self.directory)

            # Parameters directory
            self._directory = res_dir.__next__().path + '/'
            params_dir = self.directory + "parameters/"

        else:
            self._directory += 'res' + str(self.loop['active'] + 1) + '_' + self.loop['var'] + '/'
            params_dir = self._directory + "parameters/"
        #--

        #
        # Read parameters
        #

        #
        # Atom
        path = params_dir + "atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        #self._atom['mass'] = float(float(self.atom['mass']))

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
        #self._conds['num_sim'] = int(self._conds['num_sim'])
        #self._conds['num_bins'] = int(self._conds['num_bins'])

        #
        # Environment
        path = params_dir + "environment.csv"
        self._env = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        #self._env["s_0"] = float(self._env["s_0"])
        #self._env["w"] = float(self._env["w"])

    #
    # Get distributions
    def __get_dists(self):
        #
        # 3D-Histograms
        #--
        self._pos_3Dhist = {
            'freqs': None,\
            'dens': None,\
            'bins': None
        }


        self._vel_3Dhist = {
            'freqs': None,\
            'dens': None,\
            'bins': None
        }
        #--

        # Marginal histograms
        self._pos_hist = [{"freqs":[], "dens":[], "bins":[]} for i in range(3)]
        self._vel_hist = [{"freqs":[], "dens":[], "bins":[]} for i in range(3)]

        # Histogram of speeds
        self._speed_hist = {"freqs":[], "dens":[], "bins":[]}

        #
        # 3D-Histograms of positions
        #--
        path = self.directory + 'positions.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._pos_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((int(self.conds['num_bins']), int(self.conds['num_bins']), int(self.conds['num_bins'])))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._pos_3Dhist["dens"] = self.pos_3Dhist["freqs"] / np.sum(self.pos_3Dhist["freqs"])
            
            #
            # Bins
            self._pos_3Dhist["bins"] = np.zeros((3, int(self.conds['num_bins']))) - float(self.conds['max_r'])
            
            for i in range(3):
                for j in range(int(self.conds['num_bins'])):
                    delta = 2*float(self.conds['max_r']) / float(int(self.conds['num_bins']))
                    self._pos_3Dhist["bins"][i][j] += j*delta

            #
            # Marginal frequencies
            self._pos_hist[0]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(1, 2))
            self._pos_hist[1]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 2))
            self._pos_hist[2]["freqs"] = np.sum(self.pos_3Dhist["freqs"], axis=(0, 1))

            #
            # Defined marginals
            for i in range(3):
                #
                # Marginal densities
                self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"] / np.sum(self._pos_hist[i]["freqs"])

                #
                # Marginal bins
                self._pos_hist[i]["bins"] = - np.ones(int(self.conds['num_bins'])) * float(self.conds['max_r'])
                delta = 2*float(self.conds['max_r']) / float(int(self.conds['num_bins']))

                for j in range(int(self.conds['num_bins'])):
                    self._pos_hist[i]["bins"][j] += j*delta
        #--

        #
        # 3D-Histograms of velocities
        #--
        path = self.directory + 'velocities.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._vel_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((int(self.conds['num_bins']), int(self.conds['num_bins']), int(self.conds['num_bins'])))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._vel_3Dhist["dens"] = self.vel_3Dhist["freqs"] / np.sum(self.vel_3Dhist["freqs"])
            
            #
            # Bins
            self._vel_3Dhist["bins"] = np.zeros((3, int(self.conds['num_bins']))) - float(self.conds['max_v'])
            
            for i in range(3):
                for j in range(int(self.conds['num_bins'])):
                    delta = 2*float(self.conds['max_r']) / float(int(self.conds['num_bins']))
                    self._vel_3Dhist["bins"][i][j] += j*delta

            #
            # Marginal frequencies
            self._vel_hist[0]["freqs"] = np.sum(self.vel_3Dhist["freqs"], axis=(1, 2))
            self._vel_hist[1]["freqs"] = np.sum(self.vel_3Dhist["freqs"], axis=(0, 2))
            self._vel_hist[2]["freqs"] = np.sum(self.vel_3Dhist["freqs"], axis=(0, 1))

            #
            # Defined marginals
            for i in range(3):
                #
                # Marginal densities
                self._vel_hist[i]["dens"] = self._vel_hist[i]["freqs"] / np.sum(self._vel_hist[i]["freqs"])

                #
                # Marginal bins
                self._vel_hist[i]["bins"] = - np.ones(int(self.conds['num_bins'])) * float(self.conds['max_v'])
                delta = 2*float(self.conds['max_v']) / float(int(self.conds['num_bins']))

                for j in range(int(self.conds['num_bins'])):
                    self._vel_hist[i]["bins"][j] += j*delta
        #--

        #
        # Marginal histograms
        #--
        path = self.directory + 'marginals.csv'
        if os.path.exists(path):
            #
            # Read file
            df = pd.read_csv(path, index_col=0)

            # Check if velocity exists
            check_vel = (('vx' in df.columns) and ('vy' in df.columns) and ('vz' in df.columns))

            #
            # Frequencies
            self._pos_hist[0]["freqs"] = np.array(df['x'])
            self._pos_hist[1]["freqs"] = np.array(df['y'])
            self._pos_hist[2]["freqs"] = np.array(df['z'])

            if check_vel:
                self._vel_hist[0]["freqs"] = np.array(df['vx'])
                self._vel_hist[1]["freqs"] = np.array(df['vy'])
                self._vel_hist[2]["freqs"] = np.array(df['vz'])

            #
            # Densities and bins of marginal histograms
            for i in range(3):
                # Densities
                self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"] / np.sum(self._pos_hist[i]["freqs"])
                if check_vel: self._vel_hist[i]["dens"] = self._vel_hist[i]["freqs"] / np.sum(self._vel_hist[i]["freqs"])

                #
                # Bins
                self._pos_hist[i]["bins"] = - np.ones(int(int(self.conds['num_bins']))) * float(self.conds['max_r'])
                if check_vel: self._vel_hist[i]["bins"] = - np.ones(int(int(self.conds['num_bins']))) * float(self.conds['max_v'])
                pos_delta = 2*float(self.conds['max_r']) / float(int(self.conds['num_bins']))
                if check_vel: vel_delta = 2*float(self.conds['max_v']) / float(int(self.conds['num_bins']))

                for j in range(int(self.conds['num_bins'])):
                    self._pos_hist[i]["bins"][j] += j*pos_delta
                    if check_vel: self._vel_hist[i]["bins"][j] += j*vel_delta
        #--

        #
        # Histogram of speeds
        #--
        path = self.directory + 'speeds.csv'
        if os.path.exists(path):
            #
            # Frequencies
            self._speed_hist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy()

            # Densities
            self._speed_hist["dens"] = self._speed_hist["freqs"] / np.sum(self._speed_hist["freqs"])

            #
            # Bins
            #--
            self._speed_hist["bins"] = np.zeros(int(self.conds['num_bins']))
            speed_delta = float(self.conds['max_v']) / float(int(self.conds['num_bins']))

            for j in range(int(self.conds['num_bins'])):
                self._speed_hist["bins"][j] += j*speed_delta
            #--
        #--
    
    #
    def __get_name(self):
        #
        # Get short name
        obj_scandir = os.scandir(self.root_dir)
        self._name = ''

        for path in obj_scandir:
            str_splited = path.name.split("_")

            if str_splited[0] != 'group':
                act_code = int(str_splited[0])
                name = ""
                for j in range(1, len(str_splited)):
                    if j == 1: name += str_splited[j]
                    else: name += '_' + str_splited[j]

                if act_code == self.code:
                    self._name = name  
                    self._directory = path.path + "/"
                    break

    #
    def __get_loop(self):
        # Variables
        i = 0

        #
        # Directory
        #--
        res_dir = self.root_dir + str(self.code)

        if self.name:  
            res_dir += '_' + self._name

        res_dir += '/'

        # Scan results directory
        obj_scandir = os.scandir(res_dir)

        for obj_dir in obj_scandir:
            if i == 0: 
                var = obj_dir.name.split("_")

                for j in range(1, len(var)):
                    if j == 1: self._loop["var"] += var[j]
                    else: self._loop["var"] += '_' + var[j]

            #
            # Check environment loopings
            prohibited_variables = ["B_axial"]

            if self.loop["var"] in prohibited_variables:
                self._loop["var"] = ''

            elif self.loop["var"] in self.env.index:
                param = pd.read_csv(obj_dir.path + "/parameters/environment.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Check conditions loopings
            if self.loop["var"] in self.conds.index:
                param = pd.read_csv(obj_dir.path + "/parameters/conditions.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Check atom loopings
            prohibited_variables = ["symbol"]

            if self.loop["var"] in prohibited_variables:
                self._loop["var"] = ''

            elif self.loop["var"] in self.atom.index:
                param = pd.read_csv(obj_dir.path + "/parameters/atom.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Check transition loopings
            if self.loop["var"] in self.transition.index:
                param = pd.read_csv(obj_dir.path + "/parameters/transition.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            i += 1

        if self.loop["var"] == "delta":
            self.loop["values"] = sorted(self.loop["values"], key=(lambda x: abs(float(x))))

        else:
            self.loop["values"] = sorted(self.loop["values"], key=(lambda x: float(x)))

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
                
                if ((end - val) < 0 and step < 0) or ((end - val) > 0 and step > 0):
                    values = []
                    while val <= end:
                        values.append(val)
                        val += step
                else:
                    raise ValueError('Incorrect looping in the parameters')

            elif loop_str[4] == '{' and loop_str[-1] == '}':
                values = np.array(loop_str[5:-1].split(' '), dtype=float)

            else:
                raise ValueError('Invalid loop variable')

        return sorted(values, key=(lambda x: abs(x)))

    #
    def __new(self, code, name):
        self._code = code
        self._name = name

        # Check if results directory exists
        self._directory = self.root_dir
        if not os.path.exists(self._directory):
            os.mkdir(self._directory)

        # Create directory
        self._directory += str(self.code)
        if self.name: self._directory += '_' + self.name
        self._directory += '/'
        os.mkdir(self.directory)

        #
        # Create directories for each result (looping)

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
            if len(self.loop["var"]) > 0:
                #
                # Environment
                #--
                prohibited_variables = ["B_axial"]
                if (self.loop["var"] in self.env.index) and not (self.loop["var"] in prohibited_variables):
                    self.env[self.loop["var"]] = self.loop["values"][i]
                #--

                #
                # Conditions
                #--
                prohibited_variables = []
                if (self.loop["var"] in self.conds.index) and not (self.loop["var"] in prohibited_variables):
                    self.conds[self.loop["var"]] = self.loop["values"][i]
                #--

                #
                # Atom
                #--
                prohibited_variables = ['symbol']
                if (self.loop["var"] in self.atom.index) and not (self.loop["var"] in prohibited_variables):
                    self.atom[self.loop["var"]] = self.loop["values"][i]
                #--

                #
                # Transition
                #--
                prohibited_variables = []
                if (self.loop["var"] in self.transition.index) and not (self.loop["var"] in prohibited_variables):
                    self.transition[self.loop["var"]] = self.loop["values"][i]
                #--

            self.atom.to_csv(params_dir + "atom.csv", header="atom")
            self.transition.to_csv(params_dir + "transition.csv", header="transition")
            self.beams.to_csv(params_dir + "beams.csv", index=False)
            self.conds.to_csv(params_dir + "conditions.csv", header="conditions")
            self.env.to_csv(params_dir + "environment.csv", header="environment")

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
        self._conds['num_sim'] = int(float(self._conds['num_sim']))
        self._conds['num_bins'] = int(float(self._conds['num_bins']))
        
        #
        # Environment
        path = params_dir + "environment.csv"
        self._env = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Check looping values
        self.__set_loop()
        
    #
    # Set looping values
    def __set_loop(self):
        #
        # Check environment looping
        prohibited_variables = ["B_axial"]

        for idx in self.env.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.env[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._env[idx] = self.loop["values"][self.loop["active"]]

        #
        # Check conditions looping
        prohibited_variables = []

        for idx in self.conds.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.conds[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._conds[idx] = self.loop["values"][self.loop["active"]]

        #
        # Check atom looping
        prohibited_variables = ["symbol"]

        for idx in self.atom.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.atom[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._atom[idx] = self.loop["values"][self.loop["active"]]

        #
        # Check transition looping
        prohibited_variables = []

        for idx in self.transition.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.transition[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._transition[idx] = self.loop["values"][self.loop["active"]]

    #
    def mass_centre(self, axis=[0,1,2], fixed_loop_idx = False):
        #
        # Returns the best parameters to fit a Gaussian function
        def fit_gaussian(x, y):
            #
            # Gaussian function
            def gaussian(x, mean, std_dev): \
                return np.max(y) * np.exp(-((x - mean)/std_dev)**2 / 2)

            #
            # Guess values
            #--
            guess_values = np.zeros(2)
            guess_values[0] = np.sum(x*y) / np.sum(y)
            guess_values[1] = np.sqrt(np.sum(y * (x - guess_values[0])**2) / np.sum(y))
            #--
            #--

            popt, pcov = curve_fit(gaussian, x, y, guess_values)

            return popt

        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            # Variables
            if len(axis) == 1:
                r_c, std_r_c = fit_gaussian(self.pos_hist[axis[0]]["bins"], self.pos_hist[axis[0]]["dens"])
            
            else:
                r_c = np.zeros(len(axis))
                std_r_c = np.zeros(len(axis))

                #
                # Fit a normal Gaussian function
                for idx, val in enumerate(axis):
                    r_c[idx], std_r_c[idx] = fit_gaussian(self.pos_hist[val]["bins"], self.pos_hist[val]["dens"])
        #--

        #
        # With looping
        #--
        else:
            #
            # Set initial variables
            #--
            if len(axis) == 1:
                r_c = np.zeros(len(self.loop["values"]))
                std_r_c = np.zeros(len(self.loop["values"]))
            
            else:
                r_c = np.zeros((len(axis),len(self.loop["values"])))
                std_r_c = np.zeros((len(axis),len(self.loop["values"])))
            #--

            # Mass centre for each looping value
            for i in range(len(self.loop["values"])):
                #
                # Get looping values
                self.loop_idx(i)

                if len(axis) > 1:
                    for j in range(len(axis)):
                        r_c[j][i], std_r_c[j][i] = fit_gaussian(self.pos_hist[axis[j]]["bins"], self.pos_hist[axis[j]]["dens"])

                else:
                    x = self.pos_hist[axis[0]]["bins"]
                    p = self.pos_hist[axis[0]]["dens"]

                    r_c[i], std_r_c[i] = fit_gaussian(self.pos_hist[axis[0]]["bins"], self.pos_hist[axis[0]]["dens"])
        #--

        return r_c, std_r_c      

    #
    def average_velocity(self, axis=[0,1,2], fixed_loop_idx = False):
        #
        # Returns the best parameters to fit a Gaussian function
        def fit_gaussian(x, y):
            #
            # Gaussian function
            def gaussian(x, mean, std_dev): \
                return np.max(y) * np.exp(-((x - mean)/std_dev)**2 / 2)

            #
            # Guess values
            #--
            guess_values = np.zeros(2)
            guess_values[0] = np.sum(x*y) / np.sum(y)
            guess_values[1] = np.sqrt(np.sum(y * (x - guess_values[0])**2) / np.sum(y))
            #--
            #--

            popt, pcov = curve_fit(gaussian, x, y, guess_values)

            return popt

        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            # Variables
            if len(axis) == 1:
                vel_c, std_vel_c = fit_gaussian(self.vel_hist[axis[0]]["bins"], self.vel_hist[axis[0]]["dens"])
            
            else:
                vel_c = np.zeros(len(axis))
                std_vel_c = np.zeros(len(axis))

                #
                # Fit a normal Gaussian function
                for idx, val in enumerate(axis):
                    vel_c[idx], std_vel_c[idx] = fit_gaussian(self.vel_hist[val]["bins"], self.vel_hist[val]["dens"])
        #--

        #
        # With looping
        #--
        else:
            #
            # Set initial variables
            #--
            if len(axis) == 1:
                vel_c = np.zeros(len(self.loop["values"]))
                std_vel_c = np.zeros(len(self.loop["values"]))
            
            else:
                vel_c = np.zeros((len(axis),len(self.loop["values"])))
                std_vel_c = np.zeros((len(axis),len(self.loop["values"])))
            #--

            # Mass centre for each looping value
            for i in range(len(self.loop["values"])):
                #
                # Get looping values
                self.loop_idx(i)

                if len(axis) > 1:
                    for j in axis:
                        vel_c[j][i], std_vel_c[j][i] = fit_gaussian(self.vel_hist[j]["bins"], self.vel_hist[j]["dens"])

                else:
                    x = self.vel_hist[axis[0]]["bins"]
                    p = self.vel_hist[axis[0]]["dens"]

                    vel_c[i], std_vel_c[i] = fit_gaussian(self.vel_hist[axis[0]]["bins"], self.vel_hist[axis[0]]["dens"])
        #--

        return vel_c, std_vel_c      

    #
    # Get temperatures in uK
    def temperature(self, fixed_loop_idx = False, method=0):
        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            v_av, v_dev = self.average_velocity(fixed_loop_idx=fixed_loop_idx)
            
            if method == 0:
                temp = ((np.sum(v_dev*1e-2)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])

            elif method == 2:
                x = self.speed_hist['bins']*1e-2
                p = self.speed_hist['dens']

                temp = (np.sum(x*p)**2 * np.pi * float(self.atom['mass']) * self.ctes['u']) / (8*self.ctes['k_B'])

        #--

        #
        # With looping
        #--
        else:
            v_av, v_dev = self.average_velocity()

            if method == 0:
                temp = ((np.sum(v_dev*1e-2, axis=0)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square, axis=0) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])

            elif method == 2:
                temp = np.zeros(len(self.loop['values']))
                
                for i in range(len(self.loop['values'])):
                    self.loop_idx(i)
                    x = self.speed_hist['bins']*1e-2
                    p = self.speed_hist['dens']

                    temp[i] = (np.sum(x*p)**2 * np.pi * float(self.atom['mass']) * self.ctes['u']) / (8*self.ctes['k_B'])

        return temp   

    #
    # Doppler temperature
    def doppler_temperature(self, fixed_loop_idx = False):
        #
        # Check looping
        if self.loop["var"] == "gamma" and fixed_loop_idx:
            temp = np.zeros(len(self.loop["values"]))
            for i, gamma in enumerate(self.loop["values"]):
                temp[i] = 1e9*(self.ctes['hbar'] * gamma) / (2 * self.ctes['k_B']) # uK

        else:
            temp = 1e9*(self.ctes['h'] * self.transition['gamma']) / (2 * self.ctes['k_B']) # uK

        return temp

    #
    # Get 2D-histogram of positions removing an axis
    def pos_2Dhist(self, axis = 0, val = 0):
        #
        # Get bin index
        bin_idx = 0
        for idx, bin_size in enumerate(self.pos_3Dhist["bins"][axis]):
            if idx > 0 and (float(val) <= bin_size):
                bin_idx = idx - 1
                break

        #
        # Get bins
        if axis == 0:
            axis_label = {'y', 'z'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][1]), len(self.pos_3Dhist["bins"][2])))

        elif axis == 1:
            axis_label = {'x', 'z'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][0]), len(self.pos_3Dhist["bins"][2])))

        elif axis == 2:
            axis_label = {'x', 'y'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][0]), len(self.pos_3Dhist["bins"][1])))

        #
        # Get densities
        for i in range(len(hist)):
            for j in range(len(hist[i])):
                if axis == 0:
                    hist[i][j] = self.pos_3Dhist["dens"][bin_idx,i,j]

                elif axis == 1:
                    hist[i][j] = self.pos_3Dhist["dens"][i,bin_idx,j]

                elif axis == 2:
                    hist[i][j] = self.pos_3Dhist["dens"][i,j,bin_idx]

        return hist

    #
    # Check if code exists
    def __check_code(self, code):
        #
        # Variables
        obj_scandir = os.scandir(self.root_dir)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")

            if(str_splited[0] == "group"):
                group_dir = os.scandir(path.path)

                for res_dir in group_dir:
                    res_name = res_dir.name.split("_")
                    sim_code = int(res_name[0])

                    name = ""
                    for j in range(1, len(res_name)):
                        if j == 1: name += res_name[j]
                        else: name += '_' + res_name[j]

                    if sim_code == int(code):
                        self._root_dir = path.path + "/"
                        ret = True
                        break
            else:
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
        obj_scandir = os.scandir(self.root_dir)
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
        self.__cast_params_values()

    #
    # Add frequencies in the 3D histogram of positions
    def add_positions(self, pos_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array

        indexes = []
        values = []

        for i in range(int(self.conds['num_bins'])):
            for j in range(int(self.conds['num_bins'])):
                for k in range(int(self.conds['num_bins'])):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(pos_freqs_arr[int(self.conds['num_bins'])**2 * i + int(self.conds['num_bins'])*j + k])

        values = np.array(values)

        #
        # Save file

        path = self.directory + "/positions.csv"
        pos_freqs = pd.Series(values, index=indexes).astype("int32")
        pos_freqs.fillna(0, inplace=True)
        pos_freqs.to_csv(path)

        #
        # Update distributions
        self.__get_dists()

        #
        # Add marginal distribution files
        pos_freqs_arr = [self.pos_hist[i]["freqs"] for i in range(3)]
        self.add_marginals(pos_freqs_arr, self.vel_freqs_arr)

        #
        # Release memory

        del values
        del indexes
        del pos_freqs
        del path

        gc.collect()

    #
    # Add frequencies in the 3D histogram of velocities
    def add_velocities(self, vel_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array
        indexes = []
        values = []

        for i in range(int(self.conds['num_bins'])):
            for j in range(int(self.conds['num_bins'])):
                for k in range(int(self.conds['num_bins'])):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(vel_freqs_arr[int(self.conds['num_bins'])**2 * i + int(self.conds['num_bins'])*j + k])

        values = np.array(values)

        #
        # Save file
        path = self.directory + "/velocities.csv"
        vel_freqs = pd.Series(values, index=indexes).astype("int32")
        vel_freqs.fillna(0, inplace=True)
        vel_freqs.to_csv(path)

        #
        # Update distributions
        self.__get_dists()

        #
        # Add marginal distribution files
        vel_freqs_arr = [self.vel_hist[i]["freqs"] for i in range(3)]
        self.add_marginals(self.pos_freqs_arr, vel_freqs_arr)

        #
        # Release memory

        del values
        del indexes
        del vel_freqs
        del path

        gc.collect()

    #
    # Add frequencies in the marginal histograms of positions
    def add_marginals(self, pos_freqs_arr, vel_freqs_arr):
        data = {
            'x': pos_freqs_arr[0],\
            'y': pos_freqs_arr[1],\
            'z': pos_freqs_arr[2],\
            'vx': vel_freqs_arr[0],\
            'vy': vel_freqs_arr[1],\
            'vz': vel_freqs_arr[2]
        }

        path = self.directory + "marginals.csv"
        freqs = pd.DataFrame(data).astype("int32")
        freqs.fillna(0, inplace=True)
        freqs.to_csv(path)

        #
        # Release memory
        del freqs
        del data
        del path

        gc.collect()

    #
    # Add frequencies in the histogram of speeds
    def add_speeds(self, speed_freqs_arr):
        #
        # Save file

        path = self.directory + "/speeds.csv"
        speed_freqs = pd.Series(speed_freqs_arr).astype("int32")
        speed_freqs.fillna(0, inplace=True)
        speed_freqs.to_csv(path)

        #
        # Update distributions
        self.__get_dists()

        # Release memory
        del speed_freqs
        del path

        gc.collect()
