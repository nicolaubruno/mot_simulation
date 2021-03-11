#
# Libraries and modules
import os, sys
import pandas as pd
from datetime import datetime as dt
import matplotlib.pyplot as plt
from model import Result
import numpy as np

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
    def __init__(self, model):
        #
        # Model object
        self.__model = model

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
    def terminal_menu(self, options, header = '', footer = '', clear_screen = True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        print(self.header)
        if len(header) > 0: print(header + self.separator)

        #
        # Show options
        for key, val in options.items():
            print('%d -> %s;' % (key, val))

        print("")

        #
        # Show footer
        if len(footer) > 0: print(footer)

        #
        # Return input
        return input('Option code: ')

    #
    def terminal_input(self, description = '', clear_screen = True, header = True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        if header: print(self.header)

        #
        # Input
        if len(description) > 0: description += ": "
        input_str = input('\n' + description)

        return input_str

    #
    # Print the status of the simulation
    def print_simulation_status(self):
        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        print('')
        if self.__model.atoms_simulated == -1:
            print('Starting simulation ...\n')

        elif self.__model.atoms_simulated < self.__model.conds['num_sim']:
            print('Simulating ...\n')
            print('Atoms simulated: %d / %d' % (self.__model.atoms_simulated, self.__model.conds['num_sim']))
            print('Loop counter: %d / %d' % (self.__model.loop_counter, len(self.__model.loop_values)))
            print('')

        else:
            print('Simulation has finished ...\n')
            print('Simulation code: %d' % self.__model.sim_code)
            print("Simulation name: " + self.__model.sim_name)
            print('Atoms simulated: %d / %d' % (self.__model.atoms_simulated, self.__model.conds['num_sim']))
            print('Loop counter: %d / %d' % (self.__model.loop_counter + 1, len(self.__model.loop_values)))
            print('')

    #
    # Print the results
    def print_results(self, num = 10, clear_screen = True):
        #
        # Get results code
        res = self.__model.get_results(num)

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
            print(str(i + 1) + " - (" + str(val[0]) + ") " + str(dt.fromtimestamp(val[0])) + " " + val[1])

    #
    # View position histogram
    def position_marginal_histogram(self, res, axis=0):
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

        plt.bar(res.position_marginal_histogram[axis]["bins"], height=res.position_marginal_histogram[axis]["dens"], **style)

        #
        # Set plot
        plt.grid(linestyle="--")

        #
        # Show
        plt.show()
