# Parameters of the simulation

# Atom

- symbol	->	(char[2])	Atom symbol
- Z			->	(int) 		Atomic number
- mass		->	(float) 	Mass [Da or u]

# Beams

	We going to consider three bases to handle with the beams. The first one is the basis A = {a1, a2, a3}, which is a real constant basis define on the lab frame. The second one is the basis B = {b1, b2, b3}, which b1, b2, and b3 are complex vectors related to the polarizations sigma+, sigma-, and pi on the beams frame. The last one is the basis C = {c1, c2, c3}, which c1, c2, and c3 are also complex vectors related to the polarizations sigma+, sigma-, and pi on the magnetic field frame (c3 is parellel to the magnetic field direction). We are not interested in both real and imaginary components of the vectors b and c, we only need its module, therefore we going to consider the non-unit vector eps = (e1, e2, e3) to define the polarization of each beam in the simulation. This vector is defined on the basis B and its components only can be 0 or 1.

- delta		-> (float) 		Laser detuning in units of the transition rate (gamma)
- k_dic		-> (float[3]) 	non-unit vector on the basis A parallel to the wave vector
- eps		-> (int[3])		Polarizations on the basis B
- s_0  		-> (float)		Peak of the saturation parameter (I_peak / I_sat)
- w    		-> (float)		Waist Radius [mm]

# Conditions

- T_0       -> (float) 		Initial temperature [uK]
- B_0       -> (float) 		Magnetic field gradient [G / cm]
- g_bool    -> (int)		1- use gravity, 0 - do not use gravity
- i_max     -> (int)		Maximum number of iteration (simulation of individual atoms)
- r_max     -> (float)		Maximum distance (threshold) [cm]
- num_sim   -> (int)		Number of simulations (individual atoms) stopped by reach the maximum distance (r_max)
- num_bins  -> (int)		Number of bins in each histogram

# Constants

- h			-> (float)		Planck constant [10^{-34} J s]
- e			-> (float)		Elementary charge [10^{-19} C]
- c			-> (float)		Speed of light [10^{8} m / s]
- k_B		-> (float)		Boltzmann constant [10^{-23} J / K]
- mu_B		-> (float)		Bohr magneton [10^{-24} J / T]
- u			-> (float) 		Atomic mass unit [10^{-27} kg]

# Transition
 
- gamma		->	(float)		Transition rate [kHz / (2 pi)]
- lambda	->	(float)		Resonant wave length [nm]
- J_gnd		->	(int) 		Total angular momentum of the ground state
- J_exc		->	(int) 		Total angular momentum of the excited state
- g_gnd		->	(float)		Landè factor of the ground state
- g_exc		->	(float)		Landè factor of excited state