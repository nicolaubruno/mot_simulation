#
# Libraries and modules
import os, sys, time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime as dt
from model import Result
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
        # Simulation object
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
    def terminal_menu(self, options, header = None, footer = None, clear_screen = True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        print(self.header)
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

        if self.__simulation.loop["var"]:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.conds['num_sim']))
            print()

        else:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.conds['num_sim']))
            print("Loop " + str(self.__simulation.loop["active"]+1) + "/", end='')
            print(str(len(self.__simulation.loop['values'])) + " ", end='')
            print("(" + str(self.__simulation.loop["var"]) + " ", end='')
            print(str(self.__simulation.loop["values"][self.__simulation.loop["active"]]) + ")")
            print()

    #
    # Print the results
    def results_history(self, num = 5, clear_screen = True):
        #
        # Get results
        res = self.__simulation.get_results(num)

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
    # View position marginal histogram
    def pos_marg_hist(self, res, axis=0):
        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        plt.tight_layout()
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
        style={}

        plt.bar(res.pos_hist["margs"][axis]["bins"], height=res.pos_hist["margs"][axis]["dens"], **style)

        #
        # Set plot
        plt.grid(linestyle="--")

        #
        # Show
        plt.show()

    #
    # Print a 3D-vector
    def mass_centre(self, res):
        r_c, std_r_c = res.mass_centre()

        if res.loop["var"]:
            #
            # Clear stored plots
            plt.clf()

            #
            # Set style
            plt.style.use('seaborn-whitegrid')
            plt.tight_layout()
            plt.rcParams.update({
                    "figure.figsize": (7,6),\
                    "font.size":14,\
                    "axes.titlepad":16
                })

            #
            # Set labels
            labels = ['x', 'y', 'z']
            plt.title("Centre of mass as a function of the laser detuning")
            plt.xlabel(r"$ \Delta (2\pi \times MHz) $")
            plt.ylabel("position (cm)")

            #
            # Plot
            style={
                "linestyle":'-'
            }

            markers = ['o', '^', 's']

            delta = np.array(res.loop["values"])*(res.transition["gamma"]*1e-3)

            for i in range(3):
                plt.plot(delta, r_c[i], label=labels[i], marker=markers[i], **style)

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

