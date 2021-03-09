#
# Libraries and modules
from model import Model
from view import View
import gc
from datetime import datetime as dt

#
# Controller class
#
class Controller:
    
    ''' Attributes '''


    ''' Methods '''

    #
    def __init__(self):
        self.__model = Model()
        self.__view = View(self.__model)

    #
    def main(self):
        #
        # Header
        header = '\nMonte Carlo simulation of ultracold atoms in a Magneto-Optical Trap\n'
        header += 'Version 2.0, Authors: Bruno N. Santos, MSc; Ramon G. T. Rosa, PhD;'

        #
        # Footer
        footer = ''

        #
        # Options
        options = {\
            1:"Run Simulation",\
            2:"Parameters",\
            3:"Results"
        }

        #
        # Call menu
        call_menu = True
        while call_menu:
            #
            # Get option
            option_code = self.__view.terminal_menu(options, header, footer)
            call_menu = False

            #
            # Check option code
            if type(option_code) == int:
                option_code = int(option_code)

                #
                # Exit
                if option_code == 0:
                    print('Exiting ...', end='\n\n')
                    exit()

                #
                # Valid options
                elif option_code in options.keys(): 
                    #
                    # Run Simulation
                    if option_code == 1:
                        atoms_simulated = self.run_simulation()
                        footer = 'Simulation finished\n'
                        footer += "Atoms simulated: " + str(atoms_simulated) + '\n'
                        call_menu = True

                    #
                    # Show parameters of the simulation
                    elif option_code == 2:
                        call_menu = self.view_parameters()

                    #
                    # Show results
                    elif option_code == 3:
                        pass

                #
                # Invalid options
                else: call_menu = True
            else: call_menu = True

    #
    def run_simulation(self):
        #
        # Time update
        time = dt.now().timestamp()

        #
        # Start simulation
        self.__model.start_simulation()
        self.__view.print_simulation_status()

        #
        # Simulate atoms
        for i in range(self.__model.conds.num_sim):
            self.__model.simulate_atom()

            if (dt.now().timestamp() - time) > 0.5:
                self.__view.print_simulation_status()
                time = dt.now().timestamp()

        #
        # Finish simulation
        self.__model.finish_simulation()
        self.__view.print_simulation_status()

        # Clear memory
        del self.__model
        del self.__view
        self.__model = Model()
        self.__view = View(self.__model)

        # Return simulated atoms
        return i+1

    #
    def view_parameters(self):
        #
        # Header
        header = '\nParameters of the simulation'

        #
        # Options
        options = {\
            1 : 'Atom',\
            2 : 'Transition',\
            3 : 'Beams setup',\
            4 : 'Conditions',\
            5 : 'Initial Menu'
        }

        #
        # Call menu
        call_menu = True
        call_initial_menu = False
        while call_menu:
            option_code = self.__view.terminal_menu(options, header)
            call_menu = False

            #
            # Check option code
            if type(option_code) == int:
                option_code = int(option_code)

                #
                # Exit
                if option_code == 0:
                    print('Exiting ...', end='\n\n')
                    exit()

                #
                # Valid options
                elif option_code in options.keys(): 

                    if option_code == 1:
                        header = '\nAtom\n\n' + self.__model.atom.to_string()

                    elif option_code == 2:
                        header = '\nTransition\n\n' + self.__model.transition.to_string()

                    elif option_code == 3:
                        header = '\nBeams setup\n\n' + self.__model.beams.to_string()

                    elif option_code == 4:
                        header = '\nConditions\n\n' + self.__model.conds.to_string()

                    elif option_code == 5:
                        call_initial_menu = True

                #
                # Invalid options
                else: call_menu = True
            else: call_menu = True

            if (not 
                call_menu) and (not call_initial_menu):
                call_menu = True

        return call_initial_menu
