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

#define PI 3.14159265358979323
#define STRING_BUFFER_SIZE 1024
#define DELIM ","
#define MAX_BEAMS 16
#define ROOT_PATH "../"
//#define ROOT_PATH "model/"

#define h 6.62607004        // Planck constant [10^{-34} J s]
#define e 1.60217662        // Elementary charge [10^{-19} C]s
#define c 2.99792458        // Speed of light [10^{8} m / s]
#define k_B 1.38064852      // Boltzmann constant [10^{-23} J / K]
#define mu_B 9.274009994    // Bohr magneton [10^{-24} J / T]
#define u 1.660539040       // Atomic mass unit [10^{-27} kg]

//
// Structures
//

// Transition
typedef struct{
    double gamma;    /* Transition rate [kHz] */
    double lambda;   /* Resonant wave length */
    int J_gnd;      /* Total angular momentum of the ground state */
    int J_exc;      /* Total angular momentum of the excited state */
    double g_gnd;    /* Landè factor of the ground state */
    double g_exc;    /* Landè factor of the excited state */
} transition_t;

// Atom
typedef struct{
    char *symbol;               /* Atom symbol */
    int Z;                      /* Atomic number */
    double mass;                 /* Atomic mass */
    double *pos;                 /* Position [cm] */
    double *vel;                 /* Velocity [cm / s] */
    transition_t transition;    /* Optical transition */
} atom_t;

// Conditions
typedef struct{
    double T_0;      /* Initial temperature  */
    double B_0;      /* Magnetic Field gradient */
    int g_bool;     /* Use gravity */
    int i_max;      /* Maximum number of iteration */
    double r_max;    /* Maximum distance (threshold) */
    int num_bins;   /* Number of bins in each histogram */
} conditions_t;

// Beam
typedef struct {
    double delta;    /* Laser detuning */
    double *k_dic;   /* Wave vector direction */
    double *eps;       /* Polarization vector */
    double s_0;      /* Peak of the saturation parameter */
    double w;        /* Waist radius */
} beam_t;

// Beams setup
typedef struct {
    int num;        /* Number of beams */
    beam_t *beams;  /* All beams */
} beams_setup_t;

// Histogram
typedef struct {
    double num_bins; /* Number of bins */
    double bin_size; /* Bin size */
    double coord0;   /* Initial value */
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

// Apply the photonic recoil in the atom returning the time interval of the process
double compute_photonic_recoil(atom_t atom, beams_setup_t beams, conditions_t conds);

// Compute the momentum due to the gravitational force
int compute_gravitational_force(atom_t atom);

// Compute the momentum due the magnetic force
int compute_magnetic_force(atom_t atom, double B_0);

// Get magnetic field vector in the lab frame
double *get_magnetic_field(double B_0, double *r);

// Get scattering rate of a given transition
double get_scattering_rate(atom_t atom, beam_t beam, double *B, double transition, double **C, double **D);

//
// Utility functions
//

// Print all parameters of the simulation
int show_all_parameters(atom_t atom, conditions_t conditions, beams_setup_t beams_setup);

// Get int array from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size);

// Get double array from string in the format [f1 f2 ... fn]
double *get_double_array(char *str, int *size);

// Concatenate ROOT_PATH to a filename
char *concatenate_ROOT_PATH(char *filename);

// Generate a double random number following a Gaussian distribution given a mean and a standard deviation
double norm(double mean, double std_dev);

// Update histogram
int update_hist(histogram_t *hist, double val);

// Convert the components (module) of a polarization vector on the basis C to a basis D
double *update_polarization_vector(double *eps, double **r3_C, double **r3_D);

// Generate a orthonormal basis given a vector
double **orthonormal_basis(double *v3);