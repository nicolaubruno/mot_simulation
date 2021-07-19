#
# Libraries and modules
import os, sys, time
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
import numpy as np

from datetime import datetime as dt
from model import Results
from tqdm import tqdm

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
    def __init__(self, simulation):
        # Model objects
        self.__simulation = simulation

        # Terminal separator line
        #--
        self._separator = '\n'
        for i in range(85): self._separator += '-'

        self._header = "\nMonte Carlo simulation of narrow-line magneto-optical traps\n"
        self._header += "Version 2.0, Author: Bruno N. Santos;" + self._separator
        self._header += "\nGlobal options: -1 -> Back | 0 -> Exit" + self.separator + '\n'
        #--

    #
    def terminal_menu(self, options, header = None, footer = None, clear_screen = True, show_main_header = True, enumerated_list = False):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show headers 
        if show_main_header: print(self.header)
        if header is not None: print(header + self.separator)

        #
        # Show options
        i = 0
        for key, val in options.items():
            i += 1
            if enumerated_list: print("[%d] " % i, end='')

            if val:
                print('%d -> %s;' % (key, val))

            else:
                print('%d;' % (key))

        print('')

        #
        # Show footer
        if footer is not None: print(footer)

        #
        # Return input
        return input('Option code: ')

    #
    def terminal_input(self, description = None, clear_screen = True, header = True, footer= None):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        if header: print(self.header)

        #
        # Show footer 
        if footer is not None: print(footer)

        #
        # Input
        if description is not None: description += ": "
        else: description = ''

        input_str = input('\n' + description)

        return input_str

    # Print the status of the simulation
    def simulation_status(self, clear_screen=True):
        #
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        #
        # Show header 
        print(self.header)

        if not self.__simulation.results.loop["var"]:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.perform['num_sim']))
            print()

        else:
            print()
            print("Atoms simulated = %d / %d" % (self.__simulation.atoms_simulated, self.__simulation.results.perform['num_sim']))
            print("Loop " + str(self.__simulation.results.loop["active"]+1) + "/", end='')
            print(str(len(self.__simulation.results.loop['values'])) + " ", end='')
            print("(" + str(self.__simulation.results.loop["var"]) + " ", end='')
            print(str(self.__simulation.loop["values"][self.__simulation.results.loop["active"]]) + ")")
            print()

    # Print general information of a simulation
    def simulation_header(self, group, header = True, clear_screen = False, sim_opt = None):
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        # Show header 
        if header: print(self.header)

        # Simulation option
        if sim_opt: str_opt = ".. / " + sim_opt + " / "
        else: str_opt = ".. / "

        str_opt += "Group " + group + " / " + str(self.__simulation.results.code)
        if self.__simulation.results.name: str_opt += " (" + self.__simulation.results.name + ")"
        str_opt += self._separator
        print(str_opt)

        if self.__simulation.results.loop["var"]:
            print(self.__simulation.results.loop["var"] + " = ", end="")
            print(self.__simulation.results.loop['values'])

        print()
        print("Results " + str(self.__simulation.results.code) + " " + self.__simulation.results.name + self._separator)

    # Print the results
    def results_history(self, num = 5, clear_screen = True):
        #
        # Get results
        res = self.__simulation.available_results(num)

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
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()
    
    # Print the results groups
    def results_groups(self):
        #
        # Get results
        res = self.__simulation.available_results_groups()

        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        #
        # Show groups
        print('Results groups' + self.separator)

        for i, val in enumerate(res):
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()

    # View position marginal histogram
    def pos_marg_hist(self, res, axis=0):
        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.mass_centre(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
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

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.pos_hist[axis]["bins"], height=res.pos_hist[axis]["dens"], width=0.1, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.pos_hist[axis]["dens"])
            x = res.pos_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        plt.tight_layout()
        plt.show()

    # View velocity marginal histogram
    def vel_marg_hist(self, res, axis=0):
        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.average_velocity(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":14,\
                "axes.titlepad":16
            })

        #
        # Set labels
        labels = ['x', 'y', 'z']
        plt.title("Marginal histogram of velocities (" + labels[axis].upper() + "-axis)")
        plt.xlabel(labels[axis] + " [cm/s]")
        plt.ylabel(r"density")

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.vel_hist[axis]["bins"], height=res.vel_hist[axis]["dens"], width=0.9, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.vel_hist[axis]["dens"])
            x = res.vel_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        plt.tight_layout()
        plt.show()

    # Plot mass centre
    def mass_centre(self, res):
        r_c, std_r_c = res.mass_centre()

        #
        # Mass centre as a function of laser detuning
        if res.loop["var"]:
            #
            # Clear stored plots
            plt.clf()

            #
            # Set style
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })


            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1

                if info["unit"]:
                    label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else:
                    label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            #
            # Set labels
            markers = ['o', '^', 's']
            labels = ['x', 'y', 'z']
            plt.ylabel("centre of mass [mm]")
            delta = np.array(res.loop["values"])*(res.transition["gamma"]*1e-3)

            #
            # Plot simulated date
            for i in range(3):
                plt.plot(x_scale_factor * x, r_c[i]*10, linestyle="--", label=labels[i], marker=markers[i])

            #
            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)

            #
            # Show
            plt.close(1)
            plt.tight_layout()
            plt.show()

        #
        # Mass centre
        else:
            print()
            print("x = %f +- %f" % (r_c[0], std_r_c[0]))
            print("y = %f +- %f" % (r_c[1], std_r_c[1]))
            print("z = %f +- %f" % (r_c[2], std_r_c[2]))
            print()

    # Plot temperature
    def temperature(self, res, log_scale=1, doppler_temperature=False, method=0):
        #
        # Check looping
        if len(res.loop["var"]) > 0:
            # Data
            temp = res.temperature(method = method)*1e6 # uK

            # Clear stored plots
            plt.clf()

            #
            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            #plt.tight_layout()
            ax = plt.gca()

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            # Set label
            #--
            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1

                if info["unit"]:
                    label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else:
                    label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            # y label
            y_scale_factor, label = self.__temperature_axis(np.max(temp))
            plt.ylabel(label)
            #--

            # Plot temperature
            plt.plot(x_scale_factor*x, y_scale_factor*temp, label="Simulation", marker='o', linestyle='--')

            # Plot Doppler temperature
            if doppler_temperature:
                plt.plot(x, res.doppler_temperature()*np.ones(len(x)), label="Doppler temperature", linestyle='--', marker='', color='black')

            if log_scale == 0:
                plt.xscale('log')

            elif log_scale == 1:
                plt.yscale('log')

            elif log_scale == 2:
                plt.xscale('log')
                plt.yscale('log')

            # Set plot
            plt.grid(True, linestyle="--", which="both")
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())

            if doppler_temperature:
                plt.legend(frameon=True)

            plt.close(1)
            plt.tight_layout()
            plt.show()

        #
        # Mass centre
        else:
            temp = res.temperature(fixed_loop_idx=True, method = method)*1e6 # uK
            print()
            print("T [uK] = %f" % (temp))
            print("T_{doppler} [uK] = %f" % (res.doppler_temperature()))
            print()

    # Plot trapped atoms ratio
    def trapped_atoms_ratio(self, res, fitting = True):
        # With looping
        #--
        if len(res.loop["var"]) > 0:
            # Data
            ratio = res.trapped_atoms_ratio()

            # Set figure
            plt.clf() # Clear plots
            plt.figure(figsize=(6,5)) # Set figure size
            plt.style.use('seaborn-whitegrid') # Theme
            #plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            ax = plt.gca()

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            # Set labels
            #--          
            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, x_label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1
                x_label = r"$ " + info['symbol']
                if info["unit"]: x_label += r"\ [" + info['unit'] + r"] $"

            plt.xlabel(x_label)
            #--

            # y label
            plt.ylabel(r"$N / N_total$")
            #--

            # Plot trapped atoms ratio
            plt.plot(x_scale_factor * x, ratio, label="Simulated points", marker='o', linestyle="")

            # Fitting
            #--
            # Initial Velocity
            if fitting and res.loop["var"] in ["v_0", "T_0"]:
                if res.loop["var"] == "v_0":
                    vel_c, c = res.capture_velocity()

                elif res.loop["var"] == "T_0":
                    vel_c, c = res.capture_temperature()

                
                # Polynomial function
                deg = len(c) - 1
                def f(x):
                    y = 0.0

                    for i in range(deg+1):
                        y += c[deg - i]*x**i

                    return y

                x_fit = np.linspace(np.min(x), np.max(x), 1000)
                y_fit = np.array(list(map(f, x_fit)))

                plt.plot(vel_c, 0.5, marker="o", linestyle="", label=(r"$v_{cap} = %.2f\ cm/s $" % vel_c))
                plt.plot(x_scale_factor * x_fit, y_fit, label="Fitting", marker="", linestyle="--", color="Black")
            #--

            # Set plot
            plt.grid(True, linestyle="--")
            plt.legend(frameon = True)

            plt.close(1)
            plt.tight_layout()
            
            # Show
            print('Showing graph ...')
            plt.show()
        #--

        # Without loop
        #--
        else:
            ratio = res.trapped_atoms_ratio()
            print("\nN_trapped / N_total = %f\n" % (ratio))
        #--

    # Plot trap depth vs detuning
    def trap_depth_vs_detuning(self, results):
        # Variables
        vels = np.zeros(len(results))
        delta = []
        opt = None

        # Get capture values
        i = 0
        for code, name in results.items():
            res = Results(code)

            if len(res.loop["values"]) > 0:
                opt = res.loop["var"]
                delta.append(float(res.beams['main'].delta))

                if opt == "T_0":
                    vel_c, c = res.capture_temperature()

                elif opt == "v_0":
                    vel_c, c = res.capture_velocity()

                else:
                    raise ValueError("This visualization option requires looping in T_0 or v_0 variables")

                vels[i] = vel_c
                i += 1

            else:
                raise ValueError("This visualization option requires looping in each results set")

        # Capture velocity
        if opt == "v_0":
            y_label = r"$ v_c\ [cm/s] $"
            y_scale_factor = 1

        # Initial temperature
        if opt == "T_0":          
            y_scale_factor, y_label = self.__temperature_axis(np.max(mean)) 

        # Clear stored plots
        plt.clf()

        #
        # Set figure
        plt.figure(figsize=(5,4))
        plt.style.use('seaborn-whitegrid')
        plt.subplots_adjust(top=0.80, bottom=0.15)
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })
        #plt.tight_layout()
        ax = plt.gca()

        # Looping info
        info = res.info.loc["delta"]

        # Set label
        #--            
        if info["unit"]:
            x_label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
        else:
            x_label = r"$ " + info['symbol'] + r"$"

        # x label
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        #--

        # Plot trapped atoms ratio
        plt.plot(delta, y_scale_factor * vels, label="Simulated Points", marker='o', linestyle='')

        # Set plot
        plt.grid(True, linestyle="--")

        plt.close(1)
        plt.tight_layout()
        
        #
        # Show
        print('Showing graph ...')
        plt.show()

    # Plot escape flux of atoms
    def escape_flux_atoms(self, res):
        X = np.array(res.loop["values"])
        Y = np.array(res.escape_flux_atoms())

        #
        # Set figure
        plt.clf() # Clear previously plots
        plt.figure(figsize=(5,4))
        plt.style.use('seaborn-whitegrid')
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })

        # Looping info
        info = res.info.loc[res.loop["var"]]

        # Set labels
        #--
        # x label        
        if info["symbol"] == "T_0":
            x_scale_factor, x_label = self.__temperature_axis(np.max(X))
        else:
            x_scale_factor = 1
            x_label = r"$ " + info['symbol'] + r"$"           
            if info["unit"]: x_label += r"$\ [" + info['unit'] + r"] $"

        plt.xlabel(x_label)

        # y label
        plt.ylabel(r"$ \Phi\ [atoms\ /\ s]$")
        #--

        # Plot trapped atoms ratio
        plt.errorbar(x_scale_factor*X, Y, label="MOT On", marker='o', linestyle='')

        # Set plot
        plt.grid(True, linestyle="--")
        plt.close(1)
        plt.tight_layout()
        
        # Show
        print('Showing graph ...')
        plt.show()

    # Plot r.m.s. cloud size
    def cloud_size(self, res):
        if res.loop["var"]:
            #
            # Get data
            r_c, std_r_c = res.mass_centre()

            #
            # Clear stored plots
            plt.clf()

            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #
            # Set labels
            #--
            markers = ['o', '^', 's']
            labels = [r'$\sigma_x$', r'$\sigma_y$', r'$\sigma_z$']

            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1

                if info["unit"]:
                    label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else:
                    label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            plt.ylabel("Size [mm]")
            #--

            # Plot simulated date
            for i in range(3):
                plt.plot(x_scale_factor * x, std_r_c[i]*10, label=labels[i], linestyle="--", marker=markers[i])

            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)

            # Show
            plt.close(1)
            plt.tight_layout()
            plt.show()

        else:
            print('Visualization not implemented')
    
    # Heat map
    def heatmap(self, res, axis, val):
        #
        # Get data
        hist = res.pos_2Dhist(axis, val)*1e3;

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":12,\
                "axes.titlepad":12
            })

        #
        # Set labels
        plane_label = ['yz', 'xz', 'xy']
        axis_label = ['x', 'y', 'z']

        plt.title("Position distribution in " + plane_label[axis] + "-axis and " + axis_label[axis] + " = " + ("%.2f" % float(val)))
        plt.xlabel(plane_label[axis][0] + " [cm]")
        plt.ylabel(plane_label[axis][1] + " [cm]")

        #
        # Plot simulated date
        l = float(res.perform['max_r'])
        plt.imshow(hist, cmap="Blues", vmin=np.min(hist), vmax=np.max(hist), extent=[-l, l, -l, l])

        #
        # Show
        plt.tight_layout()
        plt.show()

    # Temperature axis
    def __temperature_axis(self, T_max, symbol="T"):
        if T_max >= 1e6:
            scale_factor = 1e-6
            label = r"$" + symbol + r"\ [K]$"

        elif T_max >= 1e3:
            scale_factor = 1e-3
            label = r"$" + symbol + r"\ [m K]$"

        else:
            scale_factor = 1
            label = r"$" + symbol + r"\ [\mu K]$"

        return scale_factor, label