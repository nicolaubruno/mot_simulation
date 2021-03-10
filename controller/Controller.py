#
# Libraries and modules
from model import Model, Result
from view import View
from datetime import datetime as dt
import gc

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
        header = "Initial menu"

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
            opt = self.__call_menu(options, header)
            call_menu = False

            #
            # Call menu again
            if opt == -1:
                call_menu = True

            #
            # Run Simulation
            if opt == 1:
                call_menu = self.run_simulation()

            #
            # Show parameters of the simulation
            elif opt == 2:
                call_menu = self.view_parameters()

            #
            # Show results
            elif opt == 3:
                call_menu = self.view_results()

    #
    def run_simulation(self):
        #
        # Call initial menu
        call_initial_menu = False

        #
        # Ask for shortname
        shortname = self.__call_input("Insert a short name for the simulation")

        #
        # Check shortname
        if shortname == '-1':
            call_initial_menu = True

        else:
            #
            # Time update
            time = dt.now().timestamp()

            #
            # Start simulation
            self.__model.start_simulation(shortname)
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

        return call_initial_menu

    #
    def view_parameters(self):
        #
        # Header
        header = 'Parameters of the simulation'

        #
        # Options
        options = {\
            1 : 'Atom',\
            2 : 'Transition',\
            3 : 'Beams setup',\
            4 : 'Conditions'
        }

        #
        # Get option
        call_menu = True
        call_initial_menu = False

        while call_menu:
            opt = self.__call_menu(options, header)

            if opt == 1:
                header = 'Atom\n\n' + self.__model.atom.to_string() + '\n'

            elif opt == 2:
                header = 'Transition\n\n' + self.__model.transition.to_string() + '\n'

            elif opt == 3:
                header = 'Beams setup\n\n' + self.__model.beams.to_string() + '\n'

            elif opt == 4:
                header = 'Conditions\n\n' + self.__model.conds.to_string() + '\n'

            elif opt == -1:
                call_initial_menu = True
                call_menu = False

        return call_initial_menu

    #
    def view_results(self):
        #
        # Call terminal
        msg = ''
        call_term = True
        back = False

        while call_term:
            call_term = False

            #
            # Show the first results
            self.__view.print_results()

            #
            # Get code
            code = self.__call_input(msg + "Insert simulation code", header=False, clear_screen=False)

            #
            # Check code
            if code.isdigit() or (code[0] == "-" and code[1:].isdigit()): 
                code = int(code)

                if code == -1:
                    back = True
                    call_term = False

                #
                # Call results menu
                elif self.__model.check_sim_code(code):
                    #
                    # Result
                    res = Result(code)

                    #
                    # Options
                    options = {
                        1 : "Position histogram",\
                        2 : "Heat map"
                    }

                    header = "Simulation " + str(res.sim_code) + " " + res.sim_name
                    
                    opt = self.__call_menu(options, header)

                    if opt == -1:
                        call_term = False

                    elif opt == 1 or opt == 2:
                        call_axis = True
                        while call_axis:
                            opts = {
                                1 : 'x-axis',\
                                2 : 'y-axis',\
                                3 : 'z-axis'
                            }

                            opt = self.__call_menu(opts, "Choose the axis")
                            
                            if opt == -1:
                                call_term = True
                                call_axis = False

                            else:
                                self.__view.position_histogram(res, opt-1)

                else:
                    msg = "Code %d invalid!\n" % code
                    call_term = True

            else:
                msg = "Invalid code!\n" 
                call_term = True

        return back

    #
    def __call_menu(self, options, header = ''):
        #
        # Variables
        call = True
        msg = ''

        #
        # Call menu
        while call:
            #
            # Get option
            opt = self.__view.terminal_menu(options, header=header, footer=msg)
            msg = "Invalid code! "

            #
            # Check option code
            if len(opt) > 0 and (opt.isdigit() or (opt[0] == "-" and opt[1:].isdigit())):
                opt = int(opt)

                #
                # Exit
                if opt == 0:
                    print('\nExiting ...', end='\n\n')
                    exit()

                elif opt in options.keys():
                    call = False

                elif opt == -1:
                    call = False

        return opt

    #
    def __call_input(self, description, header = True, clear_screen=True):
        opt = self.__view.terminal_input(description, header=header, clear_screen=clear_screen)

        #
        # Check return
        if opt == '0':
            print('\nExiting ...\n')
            exit(0)

        return opt
