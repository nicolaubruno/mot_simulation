//
// Libraries
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "vectors.h"

// 
// Constants
//

#define PI 3.14159265
#define STRING_BUFFER_SIZE 1024
#define DELIM ","
#define MAX_BEAMS 16
#define ROOT_PATH "../"
//#define ROOT_PATH "model/"

//
// Structures
//

// Transition
typedef struct{
    float gamma;    /* Transition rate */
    float lambda;   /* Resonant wave length */
    int J_gnd;      /* Total angular momentum of the ground state */
    int J_exc;      /* Total angular momentum of the excited state */
    float g_gnd;    /* Landè factor of the ground state */
    float g_exc;    /* Landè factor of the excited state */
} transition_t;

// Atom
typedef struct{
    char *symbol;               /* Atom symbol */
    int Z;                      /* Atomic number */
    float mass;                 /* Atomic mass */
    float *pos;                 /* Position */
    float *vel;                 /* Velocity */
    transition_t transition;    /* Optical transition */
} atom_t;

// Conditions
typedef struct{
    float T_0;      /* Initial temperature  */
    float B_0;      /* Magnetic Field gradient */
    int g_bool;     /* Use gravity */
    int i_max;      /* Maximum number of iteration */
    float r_max;    /* Maximum distance (threshold) */
    int num_bins;   /* Number of bins in each histogram */
} conditions_t;

// Beam
typedef struct {
    float delta;    /* Laser detuning */
    float *k_dic;   /* Wave vector direction */
    float *eps;       /* Polarization vector */
    float s_0;      /* Peak of the saturation parameter */
    float w;        /* Waist radius */
} beam_t;

// Beams setup
typedef struct {
    int num;        /* Number of beams */
    beam_t *beams;  /* All beams */
} beams_setup_t;

// Histogram
typedef struct {
    float num_bins; /* Number of bins */
    float bin_size; /* Bin size */
    float coord0;   /* Initial value */
    int *freqs;     /* Frequencies */
} histogram_t;

// Results
typedef struct{
    histogram_t *pos_hist;  /* Histogram of the position */
    int num_iters;          /* Number of iterations */
    int flag;               /* Status of the simulation */
} results_t;

/*
    Flags:
    0 - Simulations has stopped due an error
    1 - Simulations has finished
    2 - Simulation is running
*/

//
// Constants
//

// Planck constant [10^{-34} J s]
float h; 

// Elementary charge [10^{-19} C]
float e;

// Speed of light [10^{8} m / s]
float c;

// Boltzmann constant [10^{-23} J / K]
float k_B;

// Bohr magneton [10^{-24} J / T]
float mu_B;

// Atomic mass unit [10^{-27} kg]
float u;

//
// Main functions
//

// Run simulation for a single atom
int simulate_atom();

// Get atom
atom_t get_atom();

// Get transition
transition_t get_transition();

// Get conditions
conditions_t get_conditions();

// Get beams setup
beams_setup_t get_beams();

// Get physical constant from CSV files
int get_constants();

// Apply the photonic recoil in the atom returning the time interval of the process
float compute_photonic_recoil(atom_t atom, beams_setup_t beams, conditions_t conds);

// Compute the momentum due to the gravitational force
int compute_gravitational_force(atom_t atom);

// Compute the momentum due the magnetic force
int compute_magnetic_force(atom_t atom, float B_0);

// Get magnetic field vector in the lab frame
float *get_magnetic_field(float B_0, float *r);

//
// Utility functions
//

// Print all parameters of the simulation
int show_all_parameters(atom_t atom, conditions_t conditions, beams_setup_t beams_setup);

// Get int array from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size);

// Get float array from string in the format [f1 f2 ... fn]
float *get_float_array(char *str, int *size);

// Concatenate ROOT_PATH to a filename
char *concatenate_ROOT_PATH(char *filename);

// Generate a float random number following a Gaussian distribution given a mean and a standard deviation
float norm(float mean, float std_dev);

// Update histogram
int update_hist(histogram_t *hist, float val);

// Get coordinates of a polarization vector on the basis with pi transition parallel to the magnetic field B
float *update_polarization_vector(beam_t beam, float *B);

// Generate a orthonormal basis given a vector
float **orthonormal_basis(float *v3);