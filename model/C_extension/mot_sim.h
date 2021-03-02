
//
// Structures
//

// Results
typedef struct{
    float trapping_time;
    int num_iters; /* Number of iterations */
    int flag; /* Status of the simulation */
} results_t;

// Atom
typedef struct{
    char symbol[3]; /* Atom symbol */
    int Z; /* Atomic number */
    float mass; /* [float (Da or u)] Atomic mass */
} atom_t;

// Transition
typedef struct{
    float gamma; /* [float (kHz / (2*pi))] Transition rate */
    float lambda; /* [float (nm)] Resonant wave length */
    int J_gnd; /* [integer] Total angular momentum of the ground state */
    int J_exc; /* [integer] Total angular momentum of the excited state */
    float g_gnd; /* [float] Landè factor of the ground state */
    float g_exc; /* [float] Landè factor of the excited state */
} transition_t;

// Settings
typedef struct{
    int i_max; /* [positive integer] Maximum number of iteration (simulation of individual atoms) */
    float r_max; /* [float (cm)] Maximum distance (threshold) */
    float num_bins; /* [positive integer] Number of bins in the histogram */
} settings_t;

// Beam
typedef struct {
    float delta; /* [gamma] Laser detuning in units of the transition rate (gamma) related to the used transition */
    float k_dic[2]; /* [2D-array] Azimuthal and polar angles of the wave vector on the lab frame */
    float eps[3]; /* [3D-array] (sigma-, pi, sigma+) Polarization vector on the frame with pi parallel to the wave vector */
    float s_0; /* [float] Peak of the saturation parameter (I_peak / I_sat) */
    float w; /* [millimetre] "Waist Radius" */
} beam_t;

// Parameters
typedef struct{
    atom_t atom; /* Atom */
    transition_t transition; /* Transition */
    settings_t settings; /* Settings */
    float B; /* Magnetic Field */
    float g_bool; /* [0 or 1] 1 - use gravity, 0 - do not use gravity */
} parameters_t;

//
// Variables
//

// [J * s] Planck constant
float h;

// [C] Elementary charge
float e;

// [m / s] Speed of light
float c;

// [J / K] Boltzmann constant
float k_B;

// [J / T] Bohr magneton
float mu_B;

//
// Functions
//

// Read the parameters from CSV files and generate the results
int C_run();

// Read simulation parameters from CSV files
parameters_t C_read_parameters();

// Read physical constant from CSV files
int C_read_constants(float *h, float *e, float *c, float *k_B, float *mu_B);
