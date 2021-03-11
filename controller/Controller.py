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

    #
    @property
    def menu_request(self):
        return self._menu_request


    ''' Methods '''

    #
    def __init__(self):
        self.__model = Model()
        self.__view = View(self.__model)
        self._menu_request = -1

    #
    def main_menu(self):
        #
        # Header
        header = "Main menu"

        #
        # Options
        options = {\
            1:"Run Simulation",\
            2:"Parameters",\
            3:"Results"
        }

        #
        # Call menu
        while self.menu_request == -1:
            opt = self.__call_menu(options, header)

            #
            # Run Simulation
            if opt == 1:
                self.run_simulation()

            #
            # Show parameters of the simulation
            elif opt == 2:
                self.view_parameters()

            #
            # Show results
            elif opt == 3:
                self.view_results()

    #
    def run_simulation(self):
        #
        # Get shortname
        shortname = self.__call_input("Insert a short name for the simulation")

        #
        # Chose a parameter
        if self.menu_request == 1:
            #
            # Change menu request
            self._menu_request = -1

            #
            # Time update
            time = dt.now().timestamp()

            #
            # Print simulation status
            self.__model.start_simulation(shortname)
            self.__view.print_simulation_status()
            
            #
            # Simulate atoms
            while not self.__model.sim_end:
                self.__model.simulate()

                if (dt.now().timestamp() - time) > 0.5:
                    self.__view.print_simulation_status()
                    time = dt.now().timestamp()

            #
            # Finish simulation
            self.__view.print_simulation_status()

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
                    back = True
                    while back:
                        #
                        # Result
                        res = Result(code)

                        #
                        # Options
                        options = {
                            1 : "Position histogram"
                        }

                        header = "Simulation " + str(res.sim_code) + " " + res.sim_name
                        
                        opt = self.__call_menu(options, header)

                        if opt == -1:
                            back = False

                        #
                        # Position histogram
                        elif opt == 1:
                            call_loop = True
                            while call_loop:
                                call_axis = True
                                call_loop = False

                                #
                                # Check loop
                                if res.loop["var"] != "default":
                                    header = "Choose an option"

                                    idx = [i+1 for i in range(len(res.loop["values"]))]
                                    loop_idx = [res.loop["var"] + " = " + str(res.loop["values"][i]) for i in range(len(res.loop["values"]))]
                                    opts = dict(zip(idx, loop_idx))
                                    print(opts)
                                    
                                    opt = self.__call_menu(opts, header)

                                    if opt == -1:
                                        call_axis = False
                                        call_loop = False

                                    else:
                                        res.loop_idx(opt)

                                #
                                # Choose axis
                                while call_axis:
                                    opts = {
                                        1 : 'x-axis',\
                                        2 : 'y-axis',\
                                        3 : 'z-axis'
                                    }

                                    opt = self.__call_menu(opts, "Choose the axis")
                                    
                                    if opt == -1:
                                        call_axis = False
                                        if res.loop["var"] != "default":
                                            call_loop = True

                                    else:
                                        self.__view.pos_marg_hist(res, opt-1)

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
                    self._menu_request = 1
                    call = False

                elif opt == -1:
                    self._menu_request = -1
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

        elif opt == '-1':
            self._menu_request = -1

        else:
            self._menu_request = 1

        return opt