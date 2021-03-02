#
# Libraries and modules
import os

#
class View:
    
    ''' Attributes '''

    #
    @property
    def separator(self):
        return self._separator
    

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
    def initial_menu(self):
        header = '\nMonte Carlo simulation of ultracold atoms in a Magneto-Optical Trap\n'
        header += 'Version 2.0, Authors: Bruno N. Santos, MSc; Ramon G. T. Rosa, PhD;'

        options = {\
            1:"Run Simulation",\
            2:"Parameters",\
            3:"Results"
        }

        return self.terminal_menu(options, header)

    #
    def parameters(self):
        call_menu = True
        header = '\nParameters of the simulation'

        options = {\
            1 : 'Atom',\
            2 : 'Transition',\
            3 : 'Beams setup',\
            4 : 'Magnetic field',\
            5 : 'Initial conditions',\
            6 : 'Settings',\
            7 : 'Initial Menu'
        }

        while call_menu:
            option_code = self.terminal_menu(options, header)

            if option_code == 1:
                header = '\nAtom\n\n' + self.__model.atom.to_string()

            elif option_code == 2:
                header = '\nTransition\n\n' + self.__model.transition.to_string()

            elif option_code == 3:
                header = '\nBeams setup\n\n' + self.__model.beams.to_string()

            elif option_code == 4:
                header = '\nMagnetic field\n\n' + self.__model.magnetic_field.to_string()

            elif option_code == 5:
                header = '\nInitial conditions\n\n'
                header += 'Initial Temperature: ' + str(self.__model.T_0) + '\n'
                header += 'Gravity status: ' + str(self.__model.g_bool)

            elif option_code == 6:
                header = '\nSettings\n\n' + self.__model.settings.to_string()

            elif option_code == 7:
                call_menu = False

        return option_code
