#
# Libraries and modules
import os, sys, time
import pandas as pd
import matplotlib.pyplot as plt
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
        for i in range(50): self._separator += '-'

        self._header = self.separator + '\n'
        self._header += "Monte Carlo simulation of ultracold atoms in a Magneto-Optical Trap\n"
        self._header += "Version 2.0, Author: Bruno N. Santos;" + self._separator
        self._header += "\n\nGlobal options:" + self.separator + '\n'
        self._header += "-1 -> Back\n"
        self._header += " 0 -> Exit\n"

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
        for key, val in options.items():
            print('%d -> %s;' % (key, val))

        print('')

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

        input_str = input('\n' + description)

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
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.conds['num_sim']))
            print()

        else:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.conds['num_sim']))
            print("Loop " + str(self.__simulation.results.loop["active"]+1) + "/", end='')
            print(str(len(self.__simulation.results.loop['values'])) + " ", end='')
            print("(" + str(self.__simulation.results.loop["var"]) + " ", end='')
            print(str(self.__simulation.loop["values"][self.__simulation.results.loop["active"]]) + ")")
            print()

    #
    # Print general information of a simulation
    def simulation_header(self, header = True, clear_screen = False):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        if header: print(self.header)

        print()
        print("Results " + str(self.__simulation.results.code) + " " + self.__simulation.results.name + self._separator)

        if not self.__simulation.results.loop["var"]:
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.conds['num_sim']))
            print()

        else:
            print()
            print("Number of atoms simulated in each looping: " + str(self.__simulation.results.conds["num_sim"]))
            
            #
            # Looping
            if len(self.__simulation.results.loop["var"]) > 0:
                print("Number of loopings: " + str(len(self.__simulation.results.loop["values"])))
                print(self.__simulation.results.loop["var"] + " = [ ", end='')
                
                for val in self.__simulation.results.loop["values"]:
                    print("%.2f " % val, end='')

                print(']')
            print()

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
    # View position marginal histogram
    def pos_marg_hist(self, res, axis=0):
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

        #
        # Plot
        if axis in [0, 1, 2]:
            style={}
            plt.bar(res.pos_hist[axis]["bins"], height=res.pos_hist[axis]["dens"], width=0.1, **style)

        #
        # Set plot
        plt.grid(linestyle="--")

        #
        # Show
        plt.show()

    #
    # Plot mass centre
    def mass_centre(self, res):
        r_c, std_r_c = res.mass_centre()

        if len(res.loop["var"]) > 0:
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

            style={
                "linestyle":'--'
            }

            markers = ['o', '^', 's']

            #
            # Set labels
            labels = ['x', 'y', 'z']
            plt.title("Centre of mass as a function of the laser detuning")
            plt.xlabel(r"$ \Delta (2\pi \times MHz) $")
            plt.ylabel("position (cm)")
            delta = np.array(res.loop["values"])*(res.transition["gamma"]*1e-3)

            #
            # Plot simulated date
            for i in range(3):
                plt.errorbar(delta, r_c[i], yerr=std_r_c[i], label=labels[i], marker=markers[i], **style)


            #
            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)

            #
            # Show
            plt.show()

        else:
            print()
            print("x = %f +- %f" % (r_c[0], std_r_c[0]))
            print("y = %f +- %f" % (r_c[1], std_r_c[1]))
            print("z = %f +- %f" % (r_c[2], std_r_c[2]))
            print()

    #
    # Plot r.m.s. cloud size
    def cloud_size(self, res):
        #
        # Get data
        r_c, std_r_c = res.mass_centre()

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

        style={
            "linestyle":'--'
        }

        markers = ['o', '^', 's']

        #
        # Set labels
        labels = [r'$\sigma_x$', r'$\sigma_y$', r'$\sigma_z$']
        plt.title("R.M.S. cloud size as a function of the laser detuning")
        plt.xlabel(r"$ \Delta (2\pi \times MHz) $")
        plt.ylabel("Size (mm)")
        delta = np.array(res.loop["values"])*(res.transition["gamma"]*1e-3)

        #
        # Plot simulated date
        for i in range(3):
            plt.plot(delta, std_r_c[i], label=labels[i], marker=markers[i], **style)

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        plt.show()

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
        l = float(res.conds['max_r'])
        plt.imshow(hist, cmap="Blues", vmin=np.min(hist), vmax=np.max(hist), extent=[-l, l, -l, l])

        #
        # Show
        plt.show()

