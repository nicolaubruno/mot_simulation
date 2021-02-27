#
# Libraries and modules
import os
from model import Model

#
# Controller class
#
class Controller:
    
    ''' Attributes '''


    ''' Methods '''

    #
    def __init__(self):
        self.__model = Model()

        # Terminal separator line
        self._separator = ''
        for i in range(50):
            self._separator += '-'

    #
    # This method shows a menu structure on the terminal
    # options -> Python dictionary with the ID code of the option and its description
    def __menu(self, options, head=''):
        option_code = False
        call_menu = True

        while call_menu:
            #
            # Clear screen
            if os.name == 'nt':
                os.system('cls')

            else:
                os.system('clear')

            print(self._separator)
            print(head)
            print(self._separator)

            for key, val in options.items():
                print('%d - %s;' % (key, val))

            print('%d - %s;' % (0, 'Exit'))
            print(self._separator)

            option_code = input('Option code: ')
            print(self._separator)

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
        options = {\
            1:"Run simulation",\
            2:"Parameters",\
            3:"Results"
        }

        head = 'Magneto-optical trap simulation'
        option_code = self.__menu(options, head)

        if option_code == 1:
            pass

        elif option_code == 2:
            self.show_parameters()

        elif option_code == 3:
            pass

    #
    def show_parameters(self):
        call_menu = True

        options = {\
            1 : 'Atom',\
            2 : 'Transition',\
            3 : 'Beams setup',\
            4 : 'Magnetic field',\
            5 : 'Initial conditions',\
            6 : 'Settings',\
            7 : 'Initial Menu',\
        }

        head = 'Parameters'

        while call_menu:
            option_code = self.__menu(options, head)

            if option_code == 1:
                head = 'Atom\n' + self._separator  + '\n'
                head += self.__model.atom.to_string() + '\n' + self._separator + '\nParameters'

            elif option_code == 2:
                head = 'Transition\n' + self._separator  + '\n'
                head += self.__model.transition.to_string() + '\n' + self._separator + '\nParameters'

            elif option_code == 3:
                head = 'Beams setup\n' + self._separator  + '\n'
                head += self.__model.beams.to_string() + '\n' + self._separator + '\nParameters'

            elif option_code == 4:
                head = 'Magnetic field\n' + self._separator  + '\n'
                head += self.__model.magnetic_field.to_string() + '\n' + self._separator + '\nParameters'

            elif option_code == 5:
                head = 'Initial conditions\n' + self._separator  + '\n'
                head += 'Initial Temperature: ' + str(self.__model.T_0) + '\n'
                head += 'Gravity status: ' + str(self.__model.g_bool) + '\n' + self._separator + '\nParameters'

            elif option_code == 6:
                head = 'Settings\n' + self._separator  + '\n'
                head += self.__model.settings.to_string() + '\n' + self._separator + '\nParameters'

            elif option_code == 7:
                self.initial_menu()
