//
// Libraries
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

// 
// Constants
//

#define STRING_BUFFER_SIZE 1024
#define DELIM ","
#define MAX_BEAMS 16
//#define ROOT_PATH "../"
#define ROOT_PATH "model/"

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
    float gamma; /* [float (kHz / (2*pi))] Transition rate */
    float lambda; /* [float (nm)] Resonant wave length */
    int J_gnd; /* [integer] Total angular momentum of the ground state */
    int J_exc; /* [integer] Total angular momentum of the excited state */
    float g_gnd; /* [float] Landè factor of the ground state */
    float g_exc; /* [float] Landè factor of the excited state */
} transition_t;

// Conditions
typedef struct{
    float T_0; /* [float (uK)] Initial temperature  */
    float B_0; /* Magnetic Field gradient */
    int g_bool; /* [0 or 1] 1 - use gravity, 0 - do not use gravity */
    int i_max; /* [positive integer] Maximum number of iteration (simulation of individual atoms) */
    float r_max; /* [float (cm)] Maximum distance (threshold) */
    int num_bins; /* [positive integer] Number of bins in the histogram */
} conditions_t;

// Constants
typedef struct{
    float h; /* Planck constant */
    float e; /* Elementary charge */
    float c; /* Speed of light */
    float k_B; /* Boltzmann constant */
    float mu_B; /* Bohr magneton */
} constants_t;

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

// Get atom
atom_t C_get_atom();

// Get transition
transition_t C_get_transition();

// Get conditions
conditions_t C_get_conditions();

// Get beams setup
beams_setup_t C_get_beams();

// Get physical constant from CSV files
constants_t C_get_constants();

//
// Utility functions
//

// Get int array from string in the format [i1 i2 ... in]
int *C_get_int_array(char *str, int *size);

// Get float array from string in the format [f1 f2 ... fn]
float *C_get_float_array(char *str, int *size);

// Concatenate ROOT_PATH to a filename
char *concatenate_ROOT_PATH(char *filename);