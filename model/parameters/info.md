# Parameters of the simulation

# Atom

- symbol	->	(char[2])	Atom symbol
- Z			->	(int) 		Atomic number
- mass		->	(float) 	Mass [Da or u]

# Beams

	We going to use two bases to define the beams. The first one is the basis B1 = {e_x, e_y, e_z} related to the Cartesian vectors in lab frame. The second one is the basis B2 = {sigma+, pi, sigma-} related to the polarizations being pi a vector parallel to the wave vector (k_dic). The coordinates of the basis B1 are real numbers, whereas the coordinate of the basis B2 are 0 or 1 (there is or there is not polarization).

- delta		-> (float) 		Laser detuning in units of the transition rate (gamma)
- k_dic		-> (float[3]) 	Proportional vector to the wave vector on the basis B1
- eps		-> (int[3])		Polarization on the basis B2
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