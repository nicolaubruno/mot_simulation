#
# Libraries and modules
import numpy as np
import time, gc

from model import Simulation, Results
from view import View

from datetime import datetime as dt
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
                self.run_simulation(last_header = "Main Menu")

            #
            # Show parameters of the simulation
            elif opt == 2:
                self.view_parameters(last_header = "Main Menu")

            #
            # Show results
            elif opt == 3:
                self.view_results(last_header = "Main Menu")

    #
    def run_simulation(self, last_header = None):
        #
        # Menu level
        self._menu_level = 1

        #
        # Select a result group
        while self.menu_level == 1:
            # Get result group
            if last_header: header = " .. / Run Simulation / Select a results group:"
            else: header = "Run Simulation / Select a results group:"

            available_groups = self.__simulation.available_results_groups()
            results_group = self.__call_menu(available_groups, header)

            # Check back option
            if self.menu_level < 1:
                break

            # Set menu level
            self._menu_level += 1

            #
            # Get result code
            while self.menu_level == 2:
                #
                # Get shortname
                shortname = self.__call_input(".. / Group " + available_groups[results_group] + " / Short name")
                if shortname == -1: self._menu_level -= 1

                #
                # Back option
                if self.menu_level < 2:
                    break

                #
                # Chose a parameter
                opts = {
                    1: "Marginal distributions",\
                    2: "Complete 3D distribution",\
                    3: "Analysis of trapped atoms"
                }

                if shortname: show_shortname = " (" + shortname + ") "
                else: show_shortname = ""

                opt = self.__call_menu(opts, ".. / Group " + available_groups[results_group] + show_shortname + " / Choose a simulation option:")

                #
                # Back option
                if self.menu_level < 2:
                    break

                # Create a new simulation
                self.__simulation.new(shortname, opt, results_group)

                # Check looping
                loop_num = len(self.__simulation.results.loop["values"])
                if loop_num == 0: loop_num = 1

                # Update time
                check_time = dt.now().timestamp()

                # Print simulation status
                self.__view.simulation_header(group=available_groups[results_group], sim_opt=opts[opt], clear_screen=True)

                # Set progress bars
                pbars = []
                for i in range(loop_num):
                    desc = "Atoms simulated" if loop_num == 1 else self.__simulation.results.loop["var"] + " = " + ("%.2f" % self.__simulation.results.loop["values"][i])
                    pbars.append(tqdm(total=self.__simulation.results.perform["num_sim"], desc=desc, position=i))

                #
                # Run simulation
                #--
                for i in range(loop_num):
                    # Open new simulation for each looping value
                    if i > 0: self.__simulation.open(self.__simulation.results.code, i, opt)
                    
                    #
                    # Simulate atoms
                    #--
                    while self.__simulation.atoms_simulated < self.__simulation.results.perform["num_sim"]:
                        # Simulate atoms
                        times = self.__simulation.run()
                        
                        pbars[i].update(times)
                    #--

                    # Save simulation
                    self.__simulation.save()
                #--

                # Close progress bars
                for i in range(loop_num): pbars[i].close()

                #
                # Release memory
                gc.collect()

                # Information about the simulation
                #self.__view.simulation_header(group=available_groups[results_group], sim_opt=opts[opt], clear_screen=True)
                input("Enter any key to continue ... ")

                #
                # Set menu level
                self._menu_level = 0

    #
    def view_parameters(self, last_header = None):
        # Set menu level
        self._menu_level = 1

        # Header
        header = 'Parameters of the simulation'

        #
        # Options
        options = {\
            1 : 'Atom',\
            2 : 'Transition',\
            3 : 'Environment',\
            4 : 'Beams setup',\
            5 : 'Conditions'
        }

        #
        # Get option
        while self.menu_level == 1:
            opt = self.__call_menu(options, header)

            if opt == 1:
                header = 'Atom\n\n' + self.__simulation.results.atom.to_string() + '\n'

            elif opt == 2:
                header = 'Transition\n\n' + self.__simulation.results.transition.to_string() + '\n'
            
            elif opt == 3:
                header = 'Environment\n\n' + self.__simulation.results.env.to_string() + '\n'

            elif opt == 4:
                header = 'Beams setup\n\n' + self.__simulation.results.beams.to_string() + '\n'

            elif opt == 5:
                header = 'Conditions\n\n' + self.__simulation.results.perform.to_string() + '\n'

            else:
                self._menu_level = 0

    #
    def view_results(self, last_header = None):
        # Available groups
        available_groups = self.__simulation.available_results_groups()

        # Set menu level
        self._menu_level = 1

        # Looping menu
        while True:
            # Get result group
            #--
            if self.menu_level == 1:
                # Header
                if last_header: header = last_header + " / Results / Select a results group:"
                else: header = "Results / Select a results group:"
                
                # Get group
                results_group = self.__call_menu(available_groups, header, add_menu_level=True)
            #--

            # Collective or individual results set
            #--
            elif self.menu_level == 2:
                # Header
                header = ".. / Group " + available_groups[results_group] + " / Collective or Individual Results Sets"

                # Options
                opts = {
                    1: "Individual",\
                    2: "Collective"
                }

                # Check option
                opt = self.__call_menu(opts, header, add_menu_level=True)
                if opt == 1: collective = False
                elif opt == 2: collective = True
            #--

            # Visualization option of collective results and
            # get simulation code for individual results
            #--
            elif self.menu_level == 3:
                # Collective results set
                #--
                if collective:
                    # Header
                    header = ".. / (Collective) Group " + available_groups[results_group] + " / Visualization"

                    # Options
                    opts = {
                        1: "Trap depth vs detuning"
                    }

                    # Get visualization option
                    opt = self.__call_menu(opts, header, add_menu_level=True)
                #--

                # Individual results set
                #--
                elif not collective:
                    # Get code
                    header = ".. / (Individual) Group " + available_groups[results_group] + " / Simulation Code"
                    code = self.__call_menu(self.__simulation.available_results(results_group), header, order_command=True, add_menu_level=True)
                #--
            #--

            # Visualization of collective results set and
            # Visualization option of individual results
            #--
            elif self.menu_level == 4:
                # Collective results set
                #--
                if collective:
                    # Trap depth vs detuning
                    if opt == 1:
                        results = self.__simulation.available_results(results_group)
                        self.__view.trap_depth_vs_detuning(results)
                        self._menu_level -= 1
                #--

                # Individual results set
                if not collective:
                    #
                    # Result
                    res = Results(code)

                    # Header
                    header = ".. / Group " + available_groups[results_group] + " / " + str(res.code)
                    if res.name: header += " (" + str(res.name) + ")"
                    header += " / Visualizations"

                    # Options
                    options = {}

                    # Add options
                    #--
                    if len(res.pos_hist[0]["freqs"]) > 0:
                        options[1] = "Histogram of positions"
                        options[2] = "Centre of mass"
                        options[3] = "R.M.S Clouds sizes"

                    if len(res.vel_hist[0]["freqs"]) > 0:
                        options[4] = "Histogram of velocities"
                        options[5] = "Temperature"

                    if res.pos_3Dhist["dens"] is not None:
                        options[6] = "Heat map"

                    if res.trapped_atoms >= 0:
                        options[7] = "Trapped atoms ratio"
                    #--

                    # Get visualization option
                    opt = self.__call_menu(options, header, add_menu_level=True)
                #--
            #--

            # Visualization of individual results
            #--
            elif self.menu_level == 5:
                #
                # Individual results
                if not collective:
                    # Trapped atoms ratio
                    if opt == 7:
                        self.__view.trapped_atoms_ratio(res)

                        if not res.loop["var"]:
                            opt = self.__call_input("Enter with any key to continue", header = False, clear_screen=False)

                        self._menu_level -= 1
            
                # Collective results
                elif collective:
                    break
            #--

            # Leave menu
            else: break

        # Check results group
        #--
        while self.menu_level == 1:
            # Get result group
            if last_header: header = " .. / Results / Select a results group:"
            else: header = "Results / Select a results group:"
            available_groups = self.__simulation.available_results_groups()
            results_group = self.__call_menu(available_groups, header)

            # Check back option
            if self.menu_level < 1: break

            # Set menu level
            self._menu_level += 1

            # Collective or individual results set
            #--
            while self.menu_level == 2:
                # Get code
                header = ".. / Group " + available_groups[results_group] + " / Collective or Individual visualization"
                code = self.__call_menu(self.__simulation.available_results(results_group), header, order_command=True)

                # Check back option
                if self.menu_level < 2: break

                # Set menu level
                self._menu_level += 1

                # View option
                #--
                while self._menu_level == 3:
                    #
                    # Result
                    res = Results(code)

                    #
                    # Header
                    header = ".. / Group " + available_groups[results_group] + " / " + str(res.code)
                    if res.name: header += " (" + str(res.name) + ")"
                    header += " / Visualizations"

                    #
                    # Options
                    options = {}

                    #
                    # Add options
                    if len(res.pos_hist[0]["freqs"]) > 0:
                        options[1] = "Histogram of positions"
                        options[2] = "Centre of mass"
                        options[3] = "R.M.S Clouds sizes"

                    if len(res.vel_hist[0]["freqs"]) > 0:
                        options[4] = "Histogram of velocities"
                        options[5] = "Temperature"

                    if res.pos_3Dhist["dens"] is not None:
                        options[6] = "Heat map"

                    if res.trapped_atoms >= 0:
                        options[7] = "Trapped atoms ratio"

                    opt = self.__call_menu(options, header)

                    #
                    # Histogram of positions
                    if opt == 1:
                        #
                        # Set menu level
                        self._menu_level += 1
                        while self.menu_level == 4:
                            #
                            # Check loop
                            if len(res.loop["var"]) > 0:
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
                            while self.menu_level == 5:
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
                        while self.menu_level == 4:
                            if res.loop["var"]:
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
                        while self.menu_level == 4:
                            if res.loop["var"]:
                                self.__view.cloud_size(res)
                                self._menu_level -= 1

                            elif not res.loop["var"]:
                                self.__view.cloud_size(res)
                                opt = self.__call_input("Enter with any key to continue", header = False, clear_screen=False)
                                if opt != "-1": self._menu_level -= 1

                            else:
                                self._menu_level -= 1

                    #
                    # Histogram of velocities
                    elif opt == 4:
                        #
                        # Set menu level
                        self._menu_level += 1
                        while self.menu_level == 4:
                            #
                            # Check loop
                            if bool(res.loop["var"]):
                                header = "Choose an option"

                                idx = [i+1 for i in range(len(res.loop["values"]))]
                                loop_idx = [res.loop["var"] + " = " + ("%.2f" % float(res.loop["values"][i])) for i in range(len(res.loop["values"]))]
                                opts = dict(zip(idx, loop_idx))

                                opt = self.__call_menu(opts, header)
                                opt = int(opt)

                                if opt != -1: 
                                    res.loop_idx(opt-1)
                                    self._menu_level += 1

                            else: self._menu_level += 1

                            #
                            # Set menu level
                            while self.menu_level == 5:
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
                                    self.__view.vel_marg_hist(res, int(opt)-1)

                    #
                    # Temperature
                    elif opt == 5:
                        #
                        # Set menu level
                        self._menu_level += 1

                        #
                        # Menu level 3
                        while self.menu_level == 4:
                            if len(res.loop["var"]) > 0:
                                self.__view.temperature(res)
                                self._menu_level -= 1

                            elif not res.loop["var"]:
                                self.__view.temperature(res)
                                opt = self.__call_input("Enter with any key to continue", header = False, clear_screen=False)
                                if opt != "-1": self._menu_level -= 1

                            else:
                                self._menu_level -= 1

                    #
                    # Heat map
                    elif opt == 6:
                        #
                        # Set menu level
                        self._menu_level += 1

                        #
                        # Menu level 3
                        while self.menu_level == 4:
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
                    # Trapped atoms ratio
                    elif opt == 7:
                        #
                        # Set menu level
                        self._menu_level += 1

                        #
                        # Menu level 3
                        while self.menu_level == 4:
                            if len(res.loop["var"]) > 0:
                                self.__view.trapped_atoms_ratio(res)
                                self._menu_level -= 1

                            elif not res.loop["var"]:
                                self.__view.trapped_atoms_ratio(res)
                                opt = self.__call_input("Enter with any key to continue", header = False, clear_screen=False)
                                if opt != "-1": self._menu_level -= 1

                            else:
                                self._menu_level -= 1
                #--
            #--
        #--

    #
    def __call_menu(self, options, header = '', clear_screen = True, order_command = False, add_menu_level = False):
        #
        # Variables
        call = True
        msg = ''
        old_menu_level = self.menu_level

        # Call menu
        #--
        while call:
            # Get option
            opt = self.__view.terminal_menu(options, header=header, footer=msg, clear_screen=clear_screen)
            msg = "Invalid code! "
            call = False

            # Exit
            if opt == '0':
                print('\nExiting ...', end='\n\n')
                exit()

            # Back
            elif opt == "-1":
                self._menu_level -= 1

            elif opt[:2] == "_o" and len(opt) > 2:
                opts_keys = list(options.keys())
                if ((int(opt[2:]) - 1) < len(opts_keys)) and int(opt[2:]) > 0:
                    opt = opts_keys[int(opt[2:]) - 1]
                else: call = True

            # Check value
            elif len(opt) > 0 and (opt.isdigit() or (opt[0] == "-" and opt[1:].isdigit())):
                opt = int(opt)

                if not (opt in options.keys()):
                    call = True

            else: call = True
        #--

        if old_menu_level == self.menu_level and add_menu_level:
            self._menu_level += 1

        return opt

    #
    def __call_input(self, description, header = True, clear_screen=True, validation=None, back_opt=True):
        #
        # Get value
        opt = self.__view.terminal_input(description, header=header, clear_screen=clear_screen)

        #
        # Check value
        if opt == '0':
            print('\nExiting ...\n')
            exit(0)

        elif opt == "-1" and back_opt:
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
