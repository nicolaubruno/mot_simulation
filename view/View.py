#
# Libraries and modules
import os, sys
import pandas as pd

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
    

    ''' Methods '''

    #
    def __init__(self, model):
        #
        # Model object
        self.__model = model

        #
        # Terminal separator line
        self._separator = '\n--\n'

    #
    def terminal_menu(self, options, header='', footer=''):
        option_code = False
        call_menu = True

        while call_menu:
            #
            # Clear screen
            if os.name == 'nt':
                os.system('cls')

            else:
                os.system('clear')

            #
            # Show header 
            if len(header) > 0:
                print(header + '\n' + self.separator)

            #
            # Show options
            for key, val in options.items():
                print('%d - %s;' % (key, val))

            print('%d - %s;' % (0, 'Exit'))
            print(self.separator)

            #
            # Show footer
            if len(footer) > 0:         
                print(footer + '\n' + self.separator)

            #
            # Ask option code
            option_code = input('Option code: ')
            print(self.separator)

            #
            # Check option code
            if option_code.isdigit():
                option_code = int(option_code)

                if option_code == 0:
                    print('Exiting ...', end='\n\n')
                    exit()

                elif option_code in options.keys():
                    call_menu = False

                else:
                    call_menu = True
                    print('Invalid option code', end='\n\n')

            else:
                call_menu = True
                print('The option code is not a integer, please enter with an integer number!', end='\n\n')

        return option_code

    #
    # Print the status of the simulation
    def print_simulation_status(self):
        #
        # Clear screen
        if os.name == 'nt':
            os.system('cls')

        else:
            os.system('clear')

        print('')
        if self.__model.atoms_simulated == -1:
            print('Starting simulation ...\n')

        else:
            print('Simulating ...\n')
            print('Atoms simulated: %d / %d' % (self.__model.atoms_simulated, self.__model.conds['num_sim']))
            print('')
