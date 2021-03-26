#
# Libraries and modules
import sys, os, gc
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

            # Get loop
            self.__get_loop()

            # Get parameters
            self.__get_attr()

            # Get distributions
            self.__get_dists()

        # Create a new results
        else: self.__new(code, name)
        #--

    #
    # Get attributes
    def __get_attr(self):
        #
        # Change directory
        #

        self._directory = self.root_dir + str(self.code)

        if self.name:  
            self._directory += '_' + self._name

        self._directory += '/'

        if len(self.loop["var"]) > 0:
            self._directory += "res" + str(self.loop["active"] + 1) \
                + '_' + self.loop["var"] + '/'
        else:
            self._directory += "res1/"

        # Parameters directory
        params_dir = self.directory + "parameters/"

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
        # Environment
        path = params_dir + "environment.csv"
        self._env = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._env["s_0"] = float(self._env["s_0"])
        self._env["w"] = float(self._env["w"])

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
        path = self.directory + 'positions.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._pos_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._pos_3Dhist["dens"] = self.pos_3Dhist["freqs"] / np.sum(self.pos_3Dhist["freqs"])
            
            #
            # Bins
            self._pos_3Dhist["bins"] = np.zeros((3, self.conds["num_bins"])) - float(self.conds['max_r'])
            
            for i in range(3):
                for j in range(self.conds['num_bins']):
                    delta = 2*float(self.conds['max_r']) / float(self.conds['num_bins'])
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
                self._pos_hist[i]["bins"] = - np.ones(self.conds['num_bins']) * float(self.conds['max_r'])
                delta = 2*float(self.conds['max_r']) / float(self.conds['num_bins'])

                for j in range(self.conds['num_bins']):
                    self._pos_hist[i]["bins"][j] += j*delta

        #
        # 3D-Histograms of velocities
        path = self.directory + 'velocities.csv'
        if os.path.exists(path):
            #
            # Read histogram file
            self._vel_3Dhist["freqs"] = pd.read_csv(path, index_col=0, squeeze=True).to_numpy().reshape((self.conds['num_bins'], self.conds['num_bins'], self.conds['num_bins']))

            #
            # Filter frequencies considering the waist size as a threshold

            #
            # Densities
            self._vel_3Dhist["dens"] = self.vel_3Dhist["freqs"] / np.sum(self.vel_3Dhist["freqs"])
            
            #
            # Bins
            self._vel_3Dhist["bins"] = np.zeros((3, self.conds["num_bins"])) - float(self.conds['max_v'])
            
            for i in range(3):
                for j in range(self.conds['num_bins']):
                    delta = 2*float(self.conds['max_r']) / float(self.conds['num_bins'])
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
                self._vel_hist[i]["bins"] = - np.ones(self.conds['num_bins']) * float(self.conds['max_v'])
                delta = 2*float(self.conds['max_v']) / float(self.conds['num_bins'])

                for j in range(self.conds['num_bins']):
                    self._vel_hist[i]["bins"][j] += j*delta

        #
        # Marginal histograms of the positions
        if os.path.exists(self.directory + 'marginals.csv'):
            #
            # Read file
            df = pd.read_csv(self.directory + 'marginals.csv', index_col=0)

            #
            # Frequencies
            self._pos_hist[0]["freqs"] = np.array(df['x'])
            self._pos_hist[1]["freqs"] = np.array(df['y'])
            self._pos_hist[2]["freqs"] = np.array(df['z'])

            self._vel_hist[0]["freqs"] = np.array(df['vx'])
            self._vel_hist[1]["freqs"] = np.array(df['vy'])
            self._vel_hist[2]["freqs"] = np.array(df['vz'])

            #
            # Densities and bins of marginal histograms
            for i in range(3):
                # Densities
                self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"] / np.sum(self._pos_hist[i]["freqs"])
                self._vel_hist[i]["dens"] = self._vel_hist[i]["freqs"] / np.sum(self._vel_hist[i]["freqs"])

                #
                # Bins
                self._pos_hist[i]["bins"] = - np.ones(self.conds['num_bins']) * float(self.conds['max_r'])
                self._vel_hist[i]["bins"] = - np.ones(self.conds['num_bins']) * float(self.conds['max_v'])
                pos_delta = 2*float(self.conds['max_r']) / float(self.conds['num_bins'])
                vel_delta = 2*float(self.conds['max_v']) / float(self.conds['num_bins'])

                for j in range(self.conds['num_bins']):
                    self._pos_hist[i]["bins"][j] += j*pos_delta
                    self._vel_hist[i]["bins"][j] += j*vel_delta
    
    #
    def __get_name(self):
        #
        # Get short name
        obj_scandir = os.scandir(self.root_dir)
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
                self._directory = path.path + "/"
                break

    #
    def __get_loop(self):
        # Variables
        i = 0

        # Scan results directory
        obj_scandir = os.scandir(self.directory)

        for obj_dir in obj_scandir:
            if i == 0: 
                var = obj_dir.name.split("_")

                for j in range(1, len(var)):
                    if j == 1: self._loop["var"] += var[j]
                    else: self._loop["var"] += '_' + var[j]

            if self.loop["var"] == "delta":
                env = pd.read_csv(obj_dir.path + "/parameters/environment.csv", header=0, index_col=0, squeeze=True).astype(object)
                self._loop["values"].append(float(env["delta"]))

            i += 1

        self.loop["values"] = sorted(self.loop["values"], key=(lambda x: abs(x)))

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
            if self.loop["var"] == "delta":
                self.env["delta"] = float(self.loop["values"][i])

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
        self._conds['num_sim'] = int(self._conds['num_sim'])
        self._conds['num_bins'] = int(self._conds['num_bins'])
        
        #
        # Environment
        path = params_dir + "environment.csv"
        self._env = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        #
        # Check delta looping
        values = self.__get_loop_values(str(self.env["delta"]))
        if len(values) > 0:
            self._loop["var"] = "delta"
            self._loop["values"] = values
            self.env["delta"] = self.loop["values"][self.loop["active"]]

    #
    def mass_centre(self, axis=[0,1,2]):
        # Variables
        if len(self.loop["var"]) == 0:
            r_c = [0,0,0]
            std_r_c = [0,0,0]

            for i in axis:
                #
                # Get bins in the waist interval
                first_idx = 0
                last_idx = 0

                x = self.pos_hist[i]["bins"]
                p = self.pos_hist[i]["dens"]
                r_c[i] = np.sum(x*p)

            for i in axis:
                x2 = self.pos_hist[i]["bins"]*self.pos_hist[i]["bins"]
                p = self.pos_hist[i]["dens"]
                std_r_c[i] = np.sqrt(np.sum(x2*p) - r_c[i]**2)

        else:
            if len(axis) == 1:
                r_c = np.zeros(len(self.loop["values"]))
                r_c_square = np.zeros(len(self.loop["values"]))
                std_r_c = np.zeros(len(self.loop["values"]))
            
            else:
                r_c = np.zeros((len(axis),len(self.loop["values"])))
                r_c_square = np.zeros((len(axis),len(self.loop["values"])))
                std_r_c = np.zeros((len(axis),len(self.loop["values"])))

            for i in range(len(self.loop["values"])):
                self.loop_idx(i)

                if len(axis) > 1:
                    for j in axis:
                        x = self.pos_hist[j]["bins"]
                        p = self.pos_hist[j]["dens"]

                        r_c[j][i] = np.sum(x*p)
                        r_c_square[j][i] = np.sum(x*x*p)

                else:
                    x = self.pos_hist[axis[0]]["bins"]
                    p = self.pos_hist[axis[0]]["dens"]

                    r_c[i] = np.sum(x*p)
                    r_c_square[i] = np.sum(x*x*p)

            std_r_c = np.sqrt(r_c_square - r_c*r_c)

        return r_c, std_r_c      

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

    #
    # Add 3D position histogram
    def add_positions(self, pos_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array

        indexes = []
        values = []

        for i in range(self.conds['num_bins']):
            for j in range(self.conds['num_bins']):
                for k in range(self.conds['num_bins']):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(pos_freqs_arr[self.conds['num_bins']**2 * i + self.conds['num_bins']*j + k])

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
    # Add 3D velocities histogram
    def add_velocities(self, vel_freqs_arr):
        #
        # Transform the 3D-array in a 1D-array
        indexes = []
        values = []

        for i in range(self.conds['num_bins']):
            for j in range(self.conds['num_bins']):
                for k in range(self.conds['num_bins']):
                    indexes.append("[%d,%d,%d]" % (i+1, j+1, k+1))
                    values.append(vel_freqs_arr[self.conds['num_bins']**2 * i + self.conds['num_bins']*j + k])

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
    # Add marginal positions
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