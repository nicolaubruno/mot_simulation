//
// Libraries
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// 
// Constants
//

#define STRING_BUFFER_SIZE 1024
#define NUM_ACTIVE_PARAMS 4
#define DELIM ","
#define MAX_BEAMS 128

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
    char *symbol; /* Atom symbol */
    int Z; /* Atomic number */
    float mass; /* [float (Da or u)] Atomic mass */
} atom_t;

// Transition
typedef struct{
    atom_t atom; /* Atom related to transition */
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
    int num_bins; /* [positive integer] Number of bins in the histogram */
} settings_t;

// Initial Conditions
typedef struct {
    float T_0; /* [float (uK)] Initial temperature  */
    int g_bool; /* [0 or 1] 1 - use gravity, 0 - do not use gravity */
} initial_conditions_t;

// Beam
typedef struct {
    float delta; /* [gamma] Laser detuning in units of the transition rate (gamma) related to the used transition */
    float *k_dic; /* [2D-array] Azimuthal and polar angles of the wave vector on the lab frame */
    float *eps; /* [3D-array] (sigma-, pi, sigma+) Polarization vector on the frame with pi parallel to the wave vector */
    float s_0; /* [float] Peak of the saturation parameter (I_peak / I_sat) */
    float w; /* [millimetre] "Waist Radius" */
} beam_t;

// Beams setup
typedef struct {
    int num; /* Number of beams */
    beam_t *beams;
} beams_setup_t;

// Environment
typedef struct{
    beams_setup_t beams_setup; /* Beams setup */
    float B; /* Magnetic Field */
    float g_bool; /* [0 or 1] 1 - use gravity, 0 - do not use gravity */
} environment_t;

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
// Main functions
//

// Execute simulation for a single atom
int C_run();

// Get active parameters IDs
int *C_get_active_params();

// Get environment
environment_t C_get_environment(int id);

// Get transition
transition_t C_get_transition(int id);

// Get atom
atom_t C_get_atom(int id);

// Get settings
settings_t C_get_settings(int id);

// Get initial conditions
initial_conditions_t C_get_initial_conditions(int id);

// Get beams setup
beams_setup_t C_get_beams_setup(int num_beams, int *beams_id);

// Get physical constant from CSV files
int C_get_constants(float *h, float *e, float *c, float *k_B, float *mu_B);

//
// Utility functions
//

// Get int array from string in the format [i1 i2 ... in]
int C_get_int_array(char *str, int max_size, int *arr);

// Get float array from string in the format [f1 f2 ... fn]
int C_get_float_array(char *str, int max_size, float *arr);