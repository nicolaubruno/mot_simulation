#
# Libraries and modules
import os, sys, time
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
import numpy as np

from datetime import datetime as dt
from model import Results
from tqdm import tqdm

#
class View:
    
    ''' Attributes '''

    #
    @property
    def separator(self):
        return self._separator
    
    #
    @property
    def time_update(self):
        return self._time_update
    
    #
    @property
    def header(self):
        return self._header
    

    ''' Methods '''

    #
    def __init__(self, simulation):
        #
        # Model objects
        self.__simulation = simulation

        #
        # Terminal separator line
        self._separator = '\n'
        for i in range(70): self._separator += '-'

        self._header = "\nMonte Carlo simulation of narrow-line magneto-optical traps\n"
        self._header += "Version 2.0, Author: Bruno N. Santos;" + self._separator
        self._header += "\nGlobal options: -1 -> Back | 0 -> Exit" + self.separator + '\n'

    #
    def terminal_menu(self, options, header = None, footer = None, clear_screen = True, show_main_header = True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show headers 
        if show_main_header: print(self.header)
        if header is not None: print(header + self.separator)

        #
        # Show options
        cter = 0
        max_cter = len(options)

        for key, val in options.items():
            cter += 1
            if cter < max_cter: end_delim = '\n'
            else: end_delim = ''

            if val:  print('%d -> %s' % (key, val), end=end_delim)
            else: print('%d' % (key), end=end_delim)

        print(self.separator)

        #
        # Show footer
        if footer is not None: print(footer)

        #
        # Return input
        return input('Option code: ')

    #
    def terminal_input(self, description = None, clear_screen = True, header = True, footer= None):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        if header: print(self.header)

        #
        # Show footer 
        if footer is not None: print(footer)

        #
        # Input
        if description is not None: description += ": "
        else: description = ''

        input_str = input(description)

        return input_str

    #
    # Print the status of the simulation
    def simulation_status(self, clear_screen=True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        print(self.header)

        if not self.__simulation.results.loop["var"]:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.perform['num_sim']))
            print()

        else:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.perform['num_sim']))
            print("Loop " + str(self.__simulation.results.loop["active"]+1) + "/", end='')
            print(str(len(self.__simulation.results.loop['values'])) + " ", end='')
            print("(" + str(self.__simulation.results.loop["var"]) + " ", end='')
            print(str(self.__simulation.loop["values"][self.__simulation.results.loop["active"]]) + ")")
            print()

    #
    # Print general information of a simulation
    def simulation_header(self, header = True, clear_screen = True, sim_opt=None, group='', last_loop=-1):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        if header: print(self.header)

        if sim_opt: str_opt = ".. / " + sim_opt + " / "
        else: str_opt = ".. / "

        str_opt += "Group " + group + " / "
        str_opt += str(self.__simulation.results.code) + " (" + self.__simulation.results.name + ")"
        str_opt += self._separator
        print(str_opt)

        if self.__simulation.results.loop["var"]:
            print(self.__simulation.results.loop["var"] + " = ", end="")
            print(self.__simulation.results.loop['values'])

        if last_loop == -1:
            print("[0/" + str(len(self.__simulation.results.loop['values'])) + "]", end="")

        else:
            print("["+ str(last_loop) +"/" + str(len(self.__simulation.results.loop['values'])) + "] ", end="")
            for i in range(last_loop):
                print(str(self.__simulation.results.loop["values"][i]) + ", ", end="")

        print(self.separator)

    #
    # Print the results
    def results_history(self, num = 5, clear_screen = True):
        #
        # Get results
        res = self.__simulation.available_results(num)

        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        #
        # Show results
        print('Results history' + self.separator)

        for i, val in enumerate(res):
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()
    
    #
    # Print the results groups
    def results_groups(self):
        #
        # Get results
        res = self.__simulation.available_results_groups()

        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        #
        # Show groups
        print('Results groups' + self.separator)

        for i, val in enumerate(res):
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()

    #
    # View position marginal histo    def pos_marg_hist(self, res, axis=0):
        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.mass_centre(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":14,\
                "axes.titlepad":16
            })

        #
        # Set labels
        labels = ['x', 'y', 'z']
        plt.title("Marginal histogram " + labels[axis].upper() + "-position")
        plt.xlabel(labels[axis] + " (cm)")
        plt.ylabel(r"density")

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.pos_hist[axis]["bins"], height=res.pos_hist[axis]["dens"], width=0.1, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.pos_hist[axis]["dens"])
            x = res.pos_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show

        plt.show()

    #
    # View velocity marginal histogram
    def vel_marg_hist(self, res, axis=0):
        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.average_velocity(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":14,\
                "axes.titlepad":16
            })

        #
        # Set labels
        labels = ['x', 'y', 'z']
        plt.title("Marginal histogram of velocities (" + labels[axis].upper() + "-axis)")
        plt.xlabel(labels[axis] + " [cm/s]")
        plt.ylabel(r"density")

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.vel_hist[axis]["bins"], height=res.vel_hist[axis]["dens"], width=0.9, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.vel_hist[axis]["dens"])
            x = res.vel_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        print('Showing graph ...')
        plt.show()

    #
    # Plot mass centre
    def mass_centre(self, res):
        r_c, std_r_c = res.mass_centre()

        #
        # Mass centre as a function of laser detuning
        if res.loop["var"]:
            #
            # Clear stored plots
            plt.clf()

            #
            # Set style
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })


            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #plt.title("Centre of mass as a\nfunction of the " + info['name'].lower())
            if info["unit"]:
                label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
            else:
                label = r"$ " + info['symbol'] + r"$"

            #
            # Set labels
            markers = ['o', '^', 's']
            labels = ['x', 'y', 'z']
            plt.ylabel("centre of mass [mm]")
            plt.xlabel(label)
            delta = np.array(res.loop["values"])*(res.transition["gamma"]*1e-3)

            #
            # Plot simulated date
            for i in range(3):
                plt.plot(x, r_c[i]*10, linestyle="--", label=labels[i], marker=markers[i])

            #
            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)
            plt.close(1)
            plt.tight_layout()

            #
            # Show
            print('Showing graph ...')
            plt.show()

        #
        # Mass centre
        else:
            print()
            print("x = %f +- %f" % (r_c[0], std_r_c[0]))
            print("y = %f +- %f" % (r_c[1], std_r_c[1]))
            print("z = %f +- %f" % (r_c[2], std_r_c[2]))
            print()

    #
    # Plot temperature
    def temperature(self, res, log_scale=1, doppler_temperature=False, method=0):
        #
        # Check looping
        if len(res.loop["var"]) > 0:
            # Data
            temp = res.temperature(method = method)*1e6 # uK

            # Clear stored plots
            plt.clf()

            #
            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            #plt.tight_layout()
            ax = plt.gca()

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #
            # Set label
            #--
            #plt.title("Temperature as a function\nof the " + info['name'].lower())
            
            if info["unit"]:
                label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
            else:
                label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            plt.ylabel(r"T [$\mu K$]")
            #--

            # Plot temperature
            plt.plot(x, temp, label="Simulation", marker='o', linestyle='--')

            # Plot Doppler temperature
            if doppler_temperature:
                plt.plot(x, res.doppler_temperature()*np.ones(len(x)), label="Doppler temperature", linestyle='--', marker='', color='black')

            if log_scale == 0:
                plt.xscale('log')

            elif log_scale == 1:
                plt.yscale('log')

            elif log_scale == 2:
                plt.xscale('log')
                plt.yscale('log')

            # Set plot
            plt.grid(True, linestyle="--", which="both")
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())

            if doppler_temperature:
                plt.legend(frameon=True)

            plt.close(1)
            plt.tight_layout()
            
            #
            # Show
            print('Showing graph ...')
            plt.show()

        #
        # Mass centre
        else:
            temp = res.temperature(fixed_loop_idx=True, method = method)*1e6 # uK
            print()
            print("T [uK] = %f" % (temp))
            print("T_{doppler} [uK] = %f" % (res.doppler_temperature()))
            print()

    #
    # Plot temperature
    def trapped_atoms_ratio(self, res):
        #
        # With looping
        if len(res.loop["var"]) > 0:
            # Data
            ratio = res.trapped_atoms_ratio()

            # Clear stored plots
            plt.clf()

            #
            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            #plt.tight_layout()
            ax = plt.gca()

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #
            # Set label
            #--            
            if info["unit"]:
                label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
            else:
                label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            plt.ylabel(r"$N / N_{total}$")
            #--

            # Plot trapped atoms ratio
            plt.plot(x, ratio, label="Simulation", marker='o', linestyle='--')

            # Set plot
            plt.grid(True, linestyle="--")

            plt.close(1)
            plt.tight_layout()
            
            #
            # Show
            print('Showing graph ...')
            plt.show()

        #
        # Without loop
        else:
            ratio = res.trapped_atoms_ratio()
            print("\nN_trapped / N_total%f\n" % (ratio))

    #
    # Plot r.m.s. cloud size
    def cloud_size(self, res):
        if res.loop["var"]:
            #
            # Get data
            r_c, std_r_c = res.mass_centre()

            #
            # Clear stored plots
            plt.clf()

            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #
            # Set labels
            #--
            markers = ['o', '^', 's']
            labels = [r'$\sigma_x$', r'$\sigma_y$', r'$\sigma_z$']

            #plt.title("R.M.S. cloud size as a function of\n " + info['name'].lower())
            if info["unit"]:
                label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
            else:
                label = r"$ " + info['symbol'] + r"$"

            plt.ylabel("Size [cm]")
            plt.xlabel(label)
            #--

            # Plot simulated date
            for i in range(3):
                plt.plot(x, std_r_c[i], label=labels[i], linestyle="--", marker=markers[i])

            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)
            plt.close(1)
            
            #
            # Show
            print('Showing graph ...')
            plt.show()

        else:
            print('Visualization not implemented')
    
    #
    # Heat map
    def heatmap(self, res, axis, val):
        #
        # Get data
        hist = res.pos_2Dhist(axis, val)*1e3;

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":12,\
                "axes.titlepad":12
            })

        #
        # Set labels
        plane_label = ['yz', 'xz', 'xy']
        axis_label = ['x', 'y', 'z']

        plt.title("Position distribution in " + plane_label[axis] + "-axis and " + axis_label[axis] + " = " + ("%.2f" % float(val)))
        plt.xlabel(plane_label[axis][0] + " [cm]")
        plt.ylabel(plane_label[axis][1] + " [cm]")

        #
        # Plot simulated date
        l = float(res.perform['max_r'])
        plt.imshow(hist, cmap="Blues", vmin=np.min(hist), vmax=np.max(hist), extent=[-l, l, -l, l])

        #
        # Show
        print('Showing graph ...')
        plt.show()

