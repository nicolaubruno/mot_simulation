#
# Libraries and modules
#--
import sys, os, gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.special import erf
#--

#
class Results:

    ''' Attributes '''

    # Constants in the international system
    @property
    def ctes(self):
        return self._ctes
 
    # (Series)
    @property
    def atom(self):
        return self._atom

    # (Series)
    @property
    def transition(self):
        return self._transition

    # (Series) Initial conditions
    @property
    def ini_conds(self):
        return self._ini_conds

    # (Series) Performance
    @property
    def perform(self):
        return self._perform

    # (Series) Magnetic field
    @property
    def B_params(self):
        return self._B_params

    # (Dictionary)
    @property
    def beams(self):
        return self._beams

    # (Dataframe) Information about the parameters
    @property
    def info(self):
        return self._info
    
    # Identification code
    @property
    def code(self):
        return self._code

    # Identification short name
    @property
    def name(self):
        return self._name

    # 3D-Histogram of positions
    @property
    def pos_3Dhist(self):
        return self._pos_3Dhist
    
    # Marginal histograms of positions
    @property
    def pos_hist(self):
        return self._pos_hist

    # 3D-Histogram of velocities
    @property
    def vel_3Dhist(self):
        return self._vel_3Dhist
    
    # Marginal histograms of velocities
    @property
    def vel_hist(self):
        return self._vel_hist
    
    # Trapped atoms
    @property
    def trapped_atoms(self):
        return self._trapped_atoms

    # Looping status
    @property
    def loop(self):
        return self._loop

    # Results directory
    @property
    def directory(self):
        return self._directory

    # Root directory
    @property
    def root_dir(self):
        return self._root_dir
    

    ''' Methods '''

    #
    def __init__(self, code, name = '', loop_idx = 0, results_group = None):
        # Constants in the international system
        #--
        self._ctes = {
            'u': 1.660539040e-27,\
            'k_B': 1.38064852e-23,\
            'hbar': 1.0544718e-34,\
            'h': 6.626070040e-34
        }
        #--

        # Set loop variable
        self._loop = {
            "var": '',\
            "values": [],\
            "active": int(loop_idx)
        }

        # Root dir
        self._model_dir = "model/"
        self._root_dir = self._model_dir + "results/"
        if results_group is not None and results_group != 1: 
            self._root_dir += "group_" + results_group + "/"

        # Identification
        self._code = None
        self._name = name.strip()

        # Null trapped atoms
        self._trapped_atoms = -1

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
            self.__get_log()

        # Create a new results
        else: self.__new(code, name)
        #--

        #
        # Cast values of the parameters
        self.__cast_params_values()

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

        # Magnetic field
        self._B_params['B_0'] = float(self.B_params['B_0'])

        # Initial conditions
        self._ini_conds['T_0'] = float(self.ini_conds['T_0'])
        self._ini_conds['v_0'] = float(self.ini_conds['v_0'])
        self._ini_conds['g_bool'] = int(self.ini_conds['g_bool'])

        # Performance
        self._perform['max_time'] = float(self.perform['max_time'])
        self._perform['wait_time'] = float(self.perform['wait_time'])
        self._perform['dt'] = float(self.perform['dt'])
        self._perform['max_r'] = float(self.perform['max_r'])
        self._perform['max_v'] = float(self.perform['max_v'])
        self._perform['num_sim'] = int(self.perform['num_sim'])
        self._perform['num_bins'] = int(self.perform['num_bins'])
        self._perform['parallel_tasks'] = int(self.perform['parallel_tasks'])

        # Beams (main)
        self._beams['main']['delta'] = float(self.beams['main']['delta'])
        self._beams['main']['s_0'] = float(self.beams['main']['s_0'])
        self._beams['main']['w'] = float(self.beams['main']['w'])

        # Beams (sidebands)
        self._beams['sidebands']['freq'] = float(self.beams['sidebands']['freq'])

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

    # Get attributes
    def __get_attr(self):
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

        # Read parameters
        #--

        # Atom
        path = params_dir + "atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Transition
        path = params_dir + "transition.csv"
        self._transition = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Beams
        #--
        self._beams = {'main': None, 'setup': None, 'sidebands': None}

        # Main
        path = params_dir + "beams/main.csv"
        self._beams['main'] = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Setup
        path = params_dir + "beams/setup.csv"
        self._beams['setup'] = pd.read_csv(path, header=0)
        self._beams['setup'].index += 1

        # Sidebands
        path = params_dir + "beams/sidebands.csv"
        self._beams['sidebands'] = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        #--

        # Initial conditions
        path = params_dir + "initial_conditions.csv"
        self._ini_conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Performance
        path = params_dir + "performance.csv"
        self._perform = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Magnetic field
        path = params_dir + "magnetic_field.csv"
        self._B_params = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Information about the parameters
        path = self._model_dir + "parameters/informations.csv"
        self._info = pd.read_csv(path, header=0)
        self._info.set_index("parameter", inplace=True)
        self._info.fillna("", inplace=True)
        #--

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

        #
        # 3D-Histograms of positions
        #--
        path = self.directory + 'positions.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._pos_3Dhist["freqs"] = np.array(pd.read_csv(path, index_col=0, squeeze=True)).reshape((int(self.perform['num_bins']), int(self.perform['num_bins']), int(self.perform['num_bins'])))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._pos_3Dhist["dens"] = self.pos_3Dhist["freqs"] / np.sum(self.pos_3Dhist["freqs"])
            
            #
            # Bins
            self._pos_3Dhist["bins"] = np.zeros((3, int(self.perform['num_bins']))) - float(self.perform['max_r'])
            
            for i in range(3):
                for j in range(int(self.perform['num_bins'])):
                    delta = 2*float(self.perform['max_r']) / float(int(self.perform['num_bins']))
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
                self._pos_hist[i]["bins"] = - np.ones(int(self.perform['num_bins'])) * float(self.perform['max_r'])
                delta = 2*float(self.perform['max_r']) / float(int(self.perform['num_bins']))

                for j in range(int(self.perform['num_bins'])):
                    self._pos_hist[i]["bins"][j] += j*delta
        #--

        #
        # 3D-Histograms of velocities
        #--
        path = self.directory + 'velocities.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._vel_3Dhist["freqs"] = np.array(pd.read_csv(path, index_col=0, squeeze=True)).reshape((int(self.perform['num_bins']), int(self.perform['num_bins']), int(self.perform['num_bins'])))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._vel_3Dhist["dens"] = self.vel_3Dhist["freqs"] / np.sum(self.vel_3Dhist["freqs"])
            
            #
            # Bins
            self._vel_3Dhist["bins"] = np.zeros((3, int(self.perform['num_bins']))) - float(self.perform['max_v'])
            
            for i in range(3):
                for j in range(int(self.perform['num_bins'])):
                    delta = 2*float(self.perform['max_r']) / float(int(self.perform['num_bins']))
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
                self._vel_hist[i]["bins"] = - np.ones(int(self.perform['num_bins'])) * float(self.perform['max_v'])
                delta = 2*float(self.perform['max_v']) / float(int(self.perform['num_bins']))

                for j in range(int(self.perform['num_bins'])):
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
                self._pos_hist[i]["bins"] = - np.ones(int(int(self.perform['num_bins']))) * float(self.perform['max_r'])
                if check_vel: self._vel_hist[i]["bins"] = - np.ones(int(int(self.perform['num_bins']))) * float(self.perform['max_v'])
                pos_delta = 2*float(self.perform['max_r']) / float(int(self.perform['num_bins']))
                if check_vel: vel_delta = 2*float(self.perform['max_v']) / float(int(self.perform['num_bins']))

                for j in range(int(self.perform['num_bins'])):
                    self._pos_hist[i]["bins"][j] += j*pos_delta
                    if check_vel: self._vel_hist[i]["bins"][j] += j*vel_delta
        #--

    #
    def __get_log(self):
        path = self.directory + 'log.csv'
        if os.path.exists(path):
            log = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
            self._trapped_atoms = log['trapped_atoms']

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
            # Magnetic field
            prohibited_variables = ["B_axial", "B_bias", "B_lin_grad"]

            if self.loop["var"] in prohibited_variables:
                self._loop["var"] = ''

            elif self.loop["var"] in self.B_params.index:
                param = pd.read_csv(obj_dir.path + "/parameters/magnetic_field.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Initial conditions
            prohibited_variables = ["v_0_dir"]
            if self.loop["var"] in self.ini_conds.index:
                param = pd.read_csv(obj_dir.path + "/parameters/initial_conditions.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Performance
            if self.loop["var"] in self.perform.index:
                param = pd.read_csv(obj_dir.path + "/parameters/performance.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Atom
            prohibited_variables = ["symbol"]

            if self.loop["var"] in prohibited_variables:
                self._loop["var"] = ''

            elif self.loop["var"] in self.atom.index:
                param = pd.read_csv(obj_dir.path + "/parameters/atom.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Transition
            if self.loop["var"] in self.transition.index:
                param = pd.read_csv(obj_dir.path + "/parameters/transition.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Beams (main)
            if self.loop["var"] in self.beams['main'].index:
                param = pd.read_csv(obj_dir.path + "/parameters/beams/main.csv", header=0, index_col=0, squeeze=True).astype(object) 
                self._loop["values"].append(param[self.loop["var"]])

            #
            # Beams (sidebands)
            if self.loop["var"] in self.beams['sidebands'].index:
                param = pd.read_csv(obj_dir.path + "/parameters/beams/sidebands.csv", header=0, index_col=0, squeeze=True).astype(object) 
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
        #--
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
            os.mkdir(params_dir + "beams/")

            #
            # Add loop variable 
            if len(self.loop["var"]) > 0:
                #
                # Magnetic field
                #--
                prohibited_variables = ["B_axial", "B_bias"]
                if (self.loop["var"] in self.B_params.index) and not (self.loop["var"] in prohibited_variables):
                    self.B_params[self.loop["var"]] = self.loop["values"][i]
                #--

                #
                # Initial conditions
                #--
                prohibited_variables = ["v_0_dir"]
                if (self.loop["var"] in self.ini_conds.index) and not (self.loop["var"] in prohibited_variables):
                    self.ini_conds[self.loop["var"]] = self.loop["values"][i]

                #
                # Performance
                #--
                prohibited_variables = []
                if (self.loop["var"] in self.perform.index) and not (self.loop["var"] in prohibited_variables):
                    self.perform[self.loop["var"]] = self.loop["values"][i]
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

                #
                # Beams (main)
                #--
                prohibited_variables = []
                if (self.loop["var"] in self.beams['main'].index) and not (self.loop["var"] in prohibited_variables):
                    self.beams['main'][self.loop["var"]] = self.loop["values"][i]
                #--

                #
                # Beams (sidebands)
                #--
                prohibited_variables = []
                if (self.loop["var"] in self.beams['sidebands'].index) and not (self.loop["var"] in prohibited_variables):
                    self.beams['sidebands'][self.loop["var"]] = self.loop["values"][i]
                #--

            self.atom.to_csv(params_dir + "atom.csv", header="atom")
            self.transition.to_csv(params_dir + "transition.csv", header="transition")
            self.ini_conds.to_csv(params_dir + "initial_conditions.csv", header="initial_conditions")
            self.perform.to_csv(params_dir + "performance.csv", header="performance")
            self.B_params.to_csv(params_dir + "magnetic_field.csv", header="magnetic_field")

            self.beams['main'].to_csv(params_dir + "beams/main.csv", header="beams_main")
            self.beams['setup'].to_csv(params_dir + "beams/setup.csv", index=False)
            self.beams['sidebands'].to_csv(params_dir + "beams/sidebands.csv", header="beams_sidebands")

            # Release memory
            del res_dir, params_dir

        # Set results directory
        self._directory = self.directory + "res"

        if len(self.loop["var"]) > 0: 
            self._directory += str(self.loop['active'] + 1) + '_' + self.loop["var"]

        else:
            self._directory += '1'

        self._directory += '/'
        #--

        # Release memory
        del num_res

    # Create attributes
    def __create_attr(self):
        # Parameters directory
        params_dir = self._model_dir + "parameters/"

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
        #--
        self._beams = {'main': None, 'setup': None, 'sidebands': None}

        # Main
        path = params_dir + "/beams/main.csv"
        self._beams['main'] = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)


        # Setup
        path = params_dir + "/beams/setup.csv"
        self._beams['setup'] = pd.read_csv(path, header=0)
        self._beams['setup'].index += 1

        # Sidebands
        path = params_dir + "/beams/sidebands.csv"
        self._beams['sidebands'] = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        #--

        #
        # Initial conditions
        path = params_dir + "initial_conditions.csv"
        self._ini_conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Performance
        path = params_dir + "performance.csv"
        self._perform = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._perform['num_sim'] = int(float(self._perform['num_sim']))
        self._perform['num_bins'] = int(float(self._perform['num_bins']))
        
        #
        # Magnetic field
        path = params_dir + "magnetic_field.csv"
        self._B_params = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)      

        # Check looping values
        self.__set_loop()
        
    # Set looping values
    def __set_loop(self):
        #
        # Magnetic field
        prohibited_variables = ["B_axial", "B_bias", "B_lin_grad"]

        for idx in self.B_params.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.B_params[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._B_params[idx] = self.loop["values"][self.loop["active"]]
        
        #
        # Initial conditions
        prohibited_variables = ["v_0_dir"]

        for idx in self.ini_conds.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.ini_conds[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._ini_conds[idx] = self.loop["values"][self.loop["active"]]

        #
        # Performance
        prohibited_variables = []

        for idx in self.perform.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.perform[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._perform[idx] = self.loop["values"][self.loop["active"]]

        #
        # Atom
        prohibited_variables = ["symbol"]

        for idx in self.atom.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.atom[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._atom[idx] = self.loop["values"][self.loop["active"]]

        #
        # Transition
        prohibited_variables = []

        for idx in self.transition.index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.transition[idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._transition[idx] = self.loop["values"][self.loop["active"]]

        #
        # Beams (main)
        prohibited_variables = []

        for idx in self.beams['main'].index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.beams['main'][idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._beams['main'][idx] = self.loop["values"][self.loop["active"]]

        #
        # Beams (sidebands)
        prohibited_variables = []

        for idx in self.beams['sidebands'].index:
            if not (idx in prohibited_variables):
                values = self.__get_loop_values(str(self.beams['sidebands'][idx]))
                if len(values) > 0:
                    self._loop["var"] = idx
                    self._loop["values"] = values
                    self._beams['sidebands'][idx] = self.loop["values"][self.loop["active"]]

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
            # Convert to numpy array
            x = np.array(x)
            y = np.array(y)

            #
            # Gaussian function
            def gaussian(x, mean, std_dev):
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

    # Get temperatures in uK
    def temperature(self, fixed_loop_idx = False, method=0):
        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            v_av, v_dev = self.average_velocity(axis=[0, 1, 2], fixed_loop_idx=fixed_loop_idx)
            
            if method == 0:
                temp = ((np.sum(v_dev*1e-2)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])
            
            else:
                raise ValueError("Invalid method")

        #--

        #
        # With looping
        #--
        else:
            v_av, v_dev = self.average_velocity(axis=[0, 1, 2])

            if method == 0:
                temp = ((np.sum(v_dev*1e-2, axis=0)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square, axis=0) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])

        return temp   

    # Doppler temperature
    def doppler_temperature(self, power_broadening=False, fixed_loop_idx = False):
        if power_broadening:
            alpha = np.sqrt(1 + self.beams['main']['s_0'])
        else:
            alpha = 0;

        #
        # Check looping
        if self.loop["var"] == "gamma" and fixed_loop_idx:
            temp = np.zeros(len(self.loop["values"]))
            for i, gamma in enumerate(self.loop["values"]):
                temp[i] = 1e9*(self.ctes['hbar'] * gamma * alpha) / (2 * self.ctes['k_B']) # uK

        else:
            temp = 1e9*(self.ctes['h'] * self.transition['gamma'] * alpha) / (2 * self.ctes['k_B']) # uK

        return temp

    # Trapped atoms ratio
    def trapped_atoms_ratio(self, fixed_loop_idx = False):
        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            ratio = (self.trapped_atoms / self.perform['num_sim'])

        #
        # With looping
        #--
        else:
            ratio = np.zeros(len(self.loop["values"]))

            for i in range(len(self.loop["values"])):
                self.loop_idx(i)
                ratio[i] = (self.trapped_atoms / self.perform['num_sim'])
        #--

        return ratio

    # Capture velocity
    def capture_velocity(self):
        v_mean, v_std_dev = (0,0)

        if self.loop["var"] == 'v_0' and len(self.loop["values"]) > 1:
            v_c = self.loop['values']
            ratio = np.zeros(len(self.loop['values']))

            for i, val in enumerate(self.loop["values"]):
                self.loop_idx(i)
                ratio[i] = self.trapped_atoms / (self.perform['num_sim'])

            # General complementary error function
            def general_erfc(t, mean, std_dev):
                return 1 - (erf((t - mean) / np.sqrt(2 * std_dev**2)) - erf((- mean) / np.sqrt(2 * std_dev**2))) / 2

            # Get data
            params, covs = curve_fit(general_erfc, v_c, ratio, bounds=([min(v_c), 0], [max(v_c), (max(v_c) - min(v_c))]))
            v_mean = params[0]
            v_std_dev = params[1]

            '''
            x = np.linspace(v_mean - 5*v_std_dev, v_mean + 5*v_std_dev, 1000)
            y = np.array([general_erfc(xi, v_mean, v_std_dev) for xi in x])

            plt.clf()
            plt.plot(v_c, ratio, marker="o", linestyle="")
            plt.plot(x, y)
            plt.grid(linestyle="--")
            plt.show()
            '''

        return v_mean, v_std_dev

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
        self.__get_log()
        self.__cast_params_values()

    # Add frequencies in the 3D histogram of positions
    def add_positions(self, pos_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array

        indexes = []
        values = []

        for i in range(int(self.perform['num_bins'])):
            for j in range(int(self.perform['num_bins'])):
                for k in range(int(self.perform['num_bins'])):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(pos_freqs_arr[int(self.perform['num_bins'])**2 * i + int(self.perform['num_bins'])*j + k])

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
        self.__get_log()

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

    # Add frequencies in the 3D histogram of velocities
    def add_velocities(self, vel_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array
        indexes = []
        values = []

        for i in range(int(self.perform['num_bins'])):
            for j in range(int(self.perform['num_bins'])):
                for k in range(int(self.perform['num_bins'])):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(vel_freqs_arr[int(self.perform['num_bins'])**2 * i + int(self.perform['num_bins'])*j + k])

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
        self.__get_log()

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

    # Add trapped atoms
    def add_infos(self, data):
        path = self.directory + "log.csv"
        pd.Series(data, name="log").to_csv(path, header="log")