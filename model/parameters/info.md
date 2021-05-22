# Parameters of the simulation

## Loop function:

    If some parameter has the value started with the word "loop", the 
    simulation will loop over this parameter. Its value will be replaced
    by the values described after the word "loop". We have two options:

        - [start_value end_value step]: the value start at "start_value" 
          and goes until "end_value" adding "step" in each iteration;

        - {val_1, ..., val_n}: the value start at "val_1" and goes 
          until "val_n".

    It is just allow one loop variable in the parameters. The allowed 
    variables to use the loop function are marked with "(l)".

# Atom (atom.csv)

- symbol	->	(char[2])	Atom symbol
- Z			->	(int) 		Atomic number
- mass		->	(float) 	Mass [Da or u]

# Transition (transition.csv)
 
- gamma     ->  (float)     Linewidth [kHz / (2 pi)]
- lambda    ->  (float)     Resonant wavelength [nm]
- J_gnd     ->  (int)       Total angular momentum of the ground state
- J_exc     ->  (int)       Total angular momentum of the excited state
- g_gnd     ->  (float)     Landè factor of the ground state
- g_exc     ->  (float)     Landè factor of excited state

# Beams  (beams.csv)

	We are considering two bases to handle with the beams. The first one is 
    the basis C = {c1, c2, c3}, which c1, c2, and c3 are complex vectors 
    related to the polarizations sigma+, sigma-, and pi on the beams frame.
    The second one is the basis D = {d1, d2, d3}, which d1, d2, and d3 are
    also complex vectors related to the polarizations sigma+, sigma-, and 
    pi on the magnetic field frame (d3 is parallel to the magnetic field 
    direction). We are not interested in both real and imaginary components
    of the vectors b and c, we only need  its module, therefore we consider
    the non-unit vector eps = (e1, e2, e3) to define the polarization of 
    each beam in the simulation. This vector is defined on the basis C and
    its components only can be 0 or 1.

- k_dic		-> (float[3])   non-unit vector on the basis A parallel to the wave vector
- eps		-> (int[3])     Polarizations on the basis C

# Conditions (conditions.csv)

- T_0       -> (float)      Initial temperature [uK]
- i_max     -> (int)		Maximum number of iteration (simulation of individual atoms)
- r_max     -> (float)		Maximum distance (threshold) [cm]
- num_sim   -> (int)		Number of simulations (individual atoms) stopped by reach the maximum distance (r_max)
- num_bins  -> (int)		Number of bins in each histogram axis
- ini_iters -> (int)        Initial iterations necessary for the simulation to reach the equilibrium

# Environment (conditions.csv)

- B_0       -> (float)      Axial magnetic field gradient [G / cm]
- delta     -> (float)(l)   Laser detuning in units of the transition rate (gamma)
- s_0       -> (float)      Peak of the saturation parameter (I_peak / I_sat)
- w         -> (float)      Waist Radius [cm]
- g_bool    -> (int)        1- use gravity, 0 - do not use gravity