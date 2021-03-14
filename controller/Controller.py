#
# Libraries and modules
import numpy as np
import time, gc

from model import Simulation, Result
from view import View

from datetime import datetime as dt
from multiprocessing import cpu_count, Process
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

        #
        # Chose a parameter
        if self.menu_level == 1:
            opts = {
                1: "Simulate just marginal histograms",\
                2: "Complete 3D histogram"
            }

            only_marginals = self.__call_menu(opts, "Choose a simulation option:")

            if only_marginals == 2:
                only_marginals = 0

            if only_marginals != -1:
                #
                # Start simulation
                self.__simulation.start(shortname, only_marginals)
                sim_code = self.__simulation.code

                #
                # Release memory
                del self.__simulation
                del self.__view

                self.__simulation = Simulation(sim_code, 0, only_marginals)
                self.__view = View(self.__simulation)

                #
                # Number of loopings
                loop_num = len(self.__simulation.loop["values"])
                if loop_num == 0: loop_num = 1

                #
                # Number of the parallel processes
                num_proc = cpu_count() - 1
                if num_proc == 0: num_proc = 1

                for idx in range(loop_num):
                    #
                    # Update loop counter
                    self.__simulation.loop["active"] = idx
                    self.__simulation.atoms_simulated = 0

                    if self.__simulation.loop["var"]: 
                        desc = 'Looping ' + str(idx) + \
                            " (" + self.__simulation.loop["var"] +\
                            " = " + str("%.2f" % self.__simulation.loop["values"][idx]) + ")"

                    else: desc = "Atoms simulated"

                    #
                    # Progress bar
                    all_freqs = None
                    with tqdm(total=self.__simulation.conds["num_sim"], desc=desc) as pbar:
                        with Pool(num_proc) as pool:
                            only_marg = [only_marginals for i in range(self.__simulation.conds["num_sim"])]
                            for freqs in pool.map(self.__simulation.simulate, only_marg):
                                if all_freqs is None: all_freqs = freqs
                                else: all_freqs += freqs

                                pbar.update()

                                del freqs
                                gc.collect()

                        self.__simulation.update_pos_freqs(all_freqs)
                        self.__simulation.save()
                        gc.collect()

                    #
                    # Release memory
                    del self.__simulation
                    del self.__view

                    self.__simulation = Simulation(sim_code, idx, only_marginals)
                    self.__view = View(self.__simulation)

                gc.collect()
                exit(0)

        #
        # Set menu level
        self._menu_level = 0

    #
    def bk_run_simulation(self):
        #
        # Menu level
        self._menu_level = 1

        #
        # Get shortname
        shortname = self.__call_input("Insert a short name for the simulation")

        #
        # Chose a parameter
        if self.menu_level == 1:
            opts = {
                1: "Simulate just marginal histograms",\
                2: "Complete 3D histogram"
            }

            only_marginals = self.__call_menu(opts, "Choose a simulation option:")

            if only_marginals == 2:
                only_marginals = 0

            if only_marginals != -1:
                #
                # Time update
                time = dt.now().timestamp()

                #
                # Start simulation
                self.__simulation.start(shortname, only_marginals)
                self.__view.simulation_status()

                def simulate(i):
                    self.__simulation.simulate()

                pool = Pool(process=cpu_count())
                num_sim = self.__simulation.conds['num_sim']

                while not self.__simulation.end:
                    pool.map(simulate, range(int(0.5*num_sim)))
                    #self.__view.simulation_status()
                    print(self.__simulation.atoms_simulated)
                    exit(0);

                    #
                    # Print simulation status
                    if (dt.now().timestamp() - time) > 0.5:
                        self.__view.simulation_status()
                        time = dt.now().timestamp()

                    if not self.__simulation.status:
                        self.__simulation.save()
                        self.__simulation.status = True

                #
                # Finish simulation
                self.__view.simulation_status()
                exit(0)

                '''


                #
                # Time update
                time = dt.now().timestamp()

                #
                # Process
                simulate = Process(target=self.__simulation.simulate)

                #
                # Simulate atoms
                while not self.__simulation.end:
                    simulate.start()

                    if (dt.now().timestamp() - time) > 0.5:
                        self.__view.simulation_status()
                        time = dt.now().timestamp()

                    simulate.join()
                '''
                #
                # Finish simulation
                self.__view.simulation_status()

                #
                # Release memory
                gc.collect()

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
        while self.menu_level == 1:
            #
            # Show the first results
            self.__view.results_history()

            #
            # Get code
            def val_fun(code): 
                check = code.isdigit()
                check = check or (code[0] == "-" and code[1:].isdigit())

                if check:
                    check = self.__simulation.code(code)

                return check
            
            code = self.__call_input("Insert simulation code", header=False, clear_screen=False, validation=val_fun)
            code = int(code)

            #
            # Set menu level
            if self.menu_level == 1: 
                self._menu_level += 1

            while self._menu_level == 2:
                #
                # Result
                res = Result(code)

                #
                # Header
                header = "Simulation " + str(res.sim_code) + " " + res.sim_name

                #
                # Options
                options = {
                    1: "Position histogram",\
                    2: "Centre of mass",\
                    3: "Heat map"
                }
                
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
                            loop_idx = [res.loop["var"] + " = " + str(res.loop["values"][i]) for i in range(len(res.loop["values"]))]
                            opts = dict(zip(idx, loop_idx))

                            opt = self.__call_menu(opts, header)
                            opt = int(opt)

                            if opt != -1: 
                                res.loop_idx(opt)
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

                            opt = self.__call_menu(opts, "Choose the axis")

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
                            opt = self.__call_input("Enter with any key to continue:", header = False,clear_screen=False)
                            if opt != "-1": self._menu_level -= 1

                        else:
                            self._menu_level -= 1

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