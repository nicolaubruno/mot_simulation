#
# Libraries and modules
import numpy as np
import time, gc

from model import Simulation, Results
from view import View

from datetime import datetime as dt
#from threading import Thread
from multiprocessing import Process
from pathos.multiprocessing import ProcessPool as Pool
from tqdm import tqdm

#
# Controller class
#
class Controller:
    
    ''' Attributes '''

    __slots__ = ["_menu_level", "__simulation", "__view"]

    #
    @property
    def menu_level(self):
        return self._menu_level
    

    ''' Methods '''

    #
    def __init__(self):
        self.__simulation = Simulation()
        self.__view = View(self.__simulation)
        self._menu_level = 0

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
        # Call Main menu
        while self.menu_level == 0:
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
        # Menu level
        self._menu_level = 1

        #
        # Get shortname
        shortname = self.__call_input("Insert a short name for the simulation")
        if shortname == -1: self._menu_level -= 1

        #
        # Chose a parameter
        if self.menu_level == 1:
            opts = {
                1: "Marginal distributions",\
                2: "Complete 3D distribution"
            }

            opt = self.__call_menu(opts, "Choose a simulation option:")
            if opt ==  -1: self._menu_level -= 1
            elif opt == 2: opt = 0
            else:
                # Create a new simulation
                self.__simulation.new(shortname, opt)

                # Check looping
                loop_num = len(self.__simulation.results.loop["values"])
                if loop_num == 0: loop_num = 1

                # Update time
                check_time = dt.now().timestamp()

                # Print simulation status
                self.__view.simulation_header(header=False)

                pbars = []
                for i in range(loop_num):
                    desc = "Atoms simulated" if loop_num == 1 else self.__simulation.results.loop["var"] + " = " + ("%.2f" % self.__simulation.results.loop["values"][i])
                    pbars.append(tqdm(total=self.__simulation.results.conds["num_sim"], desc=desc, position=i))

                # Run simulation
                for i in range(loop_num):
                    #
                    # Open new simulation for each looping value
                    if i > 0: self.__simulation.open(self.__simulation.results.code, i, opt)

                    #
                    # Loop all atoms
                    while self.__simulation.atoms_simulated < self.__simulation.results.conds["num_sim"]:
                        if opt == 0: times = int(512 / self.__simulation.results.conds["num_bins"])
                        else: times = int(1024 / self.__simulation.results.conds["num_bins"])
                        if times < 1: times = 1

                        #
                        # Simulate atoms
                        times = self.__simulation.run(times)
                        
                        pbars[i].update(times)

                    # Save simulation
                    self.__simulation.save()

                for i in range(loop_num): pbars[i].close()

                #
                # Release memory
                gc.collect()

                # Information of the simulation
                self.__view.simulation_header(clear_screen=True)

        #
        # Set menu level
        self._menu_level = 0

    #
    def view_parameters(self):
        #
        # Set menu level
        self._menu_level = 1

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
        while self.menu_level == 1:
            opt = self.__call_menu(options, header)

            if opt == 1:
                header = 'Atom\n\n' + self.__simulation.atom.to_string() + '\n'

            elif opt == 2:
                header = 'Transition\n\n' + self.__simulation.transition.to_string() + '\n'

            elif opt == 3:
                header = 'Beams setup\n\n' + self.__simulation.beams.to_string() + '\n'

            elif opt == 4:
                header = 'Conditions\n\n' + self.__simulation.conds.to_string() + '\n'

            else:
                self._menu_level = 0

    #
    def view_results(self):
        #
        # Set menu level
        self._menu_level = 1

        #
        # Menu level 1
        while self.menu_level == 1:
            #
            # Show the first results
            self.__view.results_history(10)

            #
            # Get code
            def val_fun(code): 
                check = code.isdigit()
                check = check or (code[0] == "-" and code[1:].isdigit())

                if check:
                    check = self.__simulation.check_results_code(code)

                return check
            
            code = self.__call_input("Insert simulation code", header=False, clear_screen=False, validation=val_fun)
            code = int(code)

            if code == -1:
                self._menu_level -= 1

            #
            # Set menu level
            self._menu_level += 1

            #
            # Menu level 2
            while self._menu_level == 2:
                #
                # Result
                res = Results(code)

                #
                # Header
                header = "Simulation " + str(res.code) + " " + res.name

                #
                # Options
                options = {
                    1: "Position histogram",\
                    2: "Centre of mass"
                }

                #
                # Add options
                if res.loop["var"] == "delta":
                    options[3] = "R.M.S Clouds sizes"

                if res.pos_3Dhist["dens"] is not None:
                    options[4] = "Heat map"
                
                opt = self.__call_menu(options, header)

                #
                # Position histogram
                if opt == 1:
                    #
                    # Set menu level
                    self._menu_level += 1
                    while self.menu_level == 3:
                        #
                        # Check loop
                        if bool(res.loop["var"]):
                            header = "Choose an option"

                            idx = [i+1 for i in range(len(res.loop["values"]))]
                            loop_idx = [res.loop["var"] + " = " + ("%.2f" % res.loop["values"][i]) for i in range(len(res.loop["values"]))]
                            opts = dict(zip(idx, loop_idx))

                            opt = self.__call_menu(opts, header)
                            opt = int(opt)

                            if opt != -1: 
                                res.loop_idx(opt-1)
                                self._menu_level += 1

                        else: self._menu_level += 1

                        #
                        # Set menu level
                        while self.menu_level == 4:
                            opts = {
                                1 : 'x-axis',\
                                2 : 'y-axis',\
                                3 : 'z-axis'
                            }
                            
                            if len(res.loop["var"]) > 0:
                                header = "(" +res.loop["var"] + " = "
                                header += ("%.2f" % res.loop["values"][res.loop["active"]])
                                header += ") Choose the axis"

                            else:
                                header = "Choose the axis"

                            opt = self.__call_menu(opts, header)

                            if (opt == "-1") and not (res.loop["var"]):
                                self._menu_level -= 1

                            elif opt != "-1": 
                                self.__view.pos_marg_hist(res, int(opt)-1)

                #
                # Centre of mass
                elif opt == 2:
                    #
                    # Set menu level
                    self._menu_level += 1
                    while self.menu_level == 3:
                        if res.loop["var"] == "delta":
                            self.__view.mass_centre(res)
                            self._menu_level -= 1

                        elif not res.loop["var"]:
                            self.__view.mass_centre(res)
                            opt = self.__call_input("Enter with any key to continue", header = False,clear_screen=False)
                            if opt != "-1": self._menu_level -= 1

                        else:
                            self._menu_level -= 1

                #
                # R.M.S. Clouds size
                elif opt == 3:
                    #
                    # Set menu level
                    self._menu_level += 1

                    #
                    # Menu level 3
                    while self.menu_level == 3:
                        if res.loop["var"] == "delta":
                            self.__view.cloud_size(res)
                            self._menu_level -= 1

                        elif not res.loop["var"]:
                            self.__view.cloud_size(res)
                            opt = self.__call_input("Enter with any key to continue", header = False, clear_screen=False)
                            if opt != "-1": self._menu_level -= 1

                        else:
                            self._menu_level -= 1

                #
                # Heat map
                elif opt == 4:
                    #
                    # Set menu level
                    self._menu_level += 1

                    #
                    # Menu level 3
                    while self.menu_level == 3:
                        #
                        # Axis
                        opts = {
                            1:"xy-axis",\
                            2:"xz-axis",\
                            3:"yz-axis"
                        }

                        header = "Choose axis"

                        opt = self.__call_menu(opts, header)

                        if opt > 0:
                            axis = 3 - opt
                            self.__view.heatmap(res, axis, 0)

                        else: self._menu_level -= 1

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
            call = False

            #
            # Check option code
            #
            # Exit
            if opt == '0':
                print('\nExiting ...', end='\n\n')
                exit()

            # Back
            elif opt == "-1":
                self._menu_level -= 1

            # Check value
            elif len(opt) > 0 and (opt.isdigit() or (opt[0] == "-" and opt[1:].isdigit())):
                opt = int(opt)

                if not (opt in options.keys()):
                    call = True

            else: call = True

        return opt

    #
    def __call_input(self, description, header = True, clear_screen=True, validation=None):
        #
        # Get value
        opt = self.__view.terminal_input(description, header=header, clear_screen=clear_screen)

        #
        # Check value
        if opt == '0':
            print('\nExiting ...\n')
            exit(0)

        elif opt == "-1":
            self._menu_level -= 1

        elif (validation is not None) and (opt != "-1"):
            while (not validation(opt)) and (opt != "-1"):
                msg = "Invalid value!"
                opt = self.__view.terminal_input(description, header=header, clear_screen=clear_screen, footer=msg)

                #
                # Check value
                if opt == '0':
                    print('\nExiting ...\n')
                    exit(0)

                elif opt == -1:
                    self._menu_level -= 1

        return opt
