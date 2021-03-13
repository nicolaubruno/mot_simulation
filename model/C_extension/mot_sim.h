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

#define STRING_BUFFER_SIZE 1024
#define DELIM ","
#define MAX_BEAMS 16
#define MAX_dt 0.01
//#define ROOT_PATH "../"
#define ROOT_PATH "model/"

#define h 6.62607004        // Planck constant [10^{-34} J s]
#define e 1.60217662        // Elementary charge [10^{-19} C]s
#define c 2.99792458        // Speed of light [10^{8} m / s]
#define k_B 1.38064852      // Boltzmann constant [10^{-23} J / K]
#define mu_B 9.274009994    // Bohr magneton [10^{-24} J / T]
#define u 1.660539040       // Atomic mass unit [10^{-27} kg]
#define g 980.665           // Gravitational acceleration [cm / s^2]
#define PI 3.14159265358

//
// Structures
//

// Transition
typedef struct{
    double gamma;   /* Transition rate [kHz] */
    double lambda;  /* Resonant wave length */
    int J_gnd;      /* Total angular momentum of the ground state */
    int J_exc;      /* Total angular momentum of the excited state */
    double g_gnd;   /* Landè factor of the ground state */
    double g_exc;   /* Landè factor of the excited state */
} transition_t;

// Atom
typedef struct{
    char *symbol;               /* Atom symbol */
    int Z;                      /* Atomic number */
    double mass;                /* Atomic mass */
    double *pos;                /* Position [cm] */
    double *vel;                /* Velocity [cm / s] */
    transition_t transition;    /* Optical transition */
} atom_t;

// Conditions
typedef struct{
    double T_0;         /* Initial temperature  */
    double B_0;         /* Magnetic Field gradient */
    double *B_axial;    /* Magnetic Field axial direction */
    double delta;       /* Laser detuning */
    int g_bool;         /* Use gravity */
    int i_max;          /* Maximum number of iteration */
    double r_max;       /* Maximum distance (threshold) */
    int num_bins;       /* Number of bins in each histogram */
} conditions_t;

// Beam
typedef struct {
    double *k_dic;   /* Wave vector direction */
    double *eps;     /* Polarization vector */
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
    double num_bins;    /* Number of bins */
    double bin_size;    /* Bin size */
    double coord0;      /* Initial value */
    int *freqs;         /* Frequencies */
} histogram_t;

// 3D-dimensional Histogram
typedef struct {
    int *num_bins;    /* Number of bins in each axis */
    double *bins_size;    /* Bin size in each axis */
    double *coord0;      /* Initial value in each axis */
    int ***freqs;        /* Frequencies */
} histogram_3d_t;

// Photon-atom scattering
typedef struct {
    double *vel;    /* Velocity gain */
    double R;       /* Scattering rate */
    double dt;      /* Process time */
} scattering_t;

// Results
typedef struct{
    histogram_3d_t pos_3Dhist;      /* Histogram of the position */
    histogram_t *pos_hist;          /* Marginal histograms */
    int num_iters;                  /* Number of iterations */
    double time;                    /* Total time [s] */
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
results_t simulate_atom(char *params_path, int marginals, long time);

// Get atom
atom_t get_atom(conditions_t conds, char *params_path);

// Get transition
transition_t get_transition(char *params_path);

// Get conditions
conditions_t get_conditions(char *params_path);

// Get beams setup
beams_setup_t get_beams(char *params_path);

// Compute scattering variables (scattering_t) in a photon-atom scattering event
scattering_t photonic_recoil(atom_t atom, beams_setup_t beams_setup, conditions_t conds, double **B);

// Get magnetic acceleration
double *magnetic_acceleration(atom_t atom, double B_0, double **B_basis);

// Get magnetic field vector in the lab frame
double *get_magnetic_field(double B_0, double **B_basis, double *r);

// Get scattering rate
double scattering_rate(atom_t atom, beam_t beam, conditions_t conds, double *B, int pol);

// Get components of a vector v on the basis B given the components on basis A
double *change_basis(double *v, double **A, double **B_basis);

// Get polarization probabilities given the magnetic field direction
double *polarization_probs(beam_t beam, double *eB);

// Print simulation status
int print_status(atom_t atom, results_t res);

// Print results
int print_results(results_t res);

// Write results in a CSV file
int write_results(char *dir_code, results_t res);

//
// Utility functions
//

// Get int array from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size);

// Get double array from string in the format [f1 f2 ... fn]
double *get_double_array(char *str, int *size);

// Concatenate ROOT_PATH to a filename
char *concatenate_ROOT_PATH(char *filename);

// Concatenate strings
char *str_concatenate(char *str1, char *str2);

// Generate a double random number following a Gaussian distribution given a mean and a standard deviation
double random_norm(double mean, double std_dev);

// Generate a double random number following a Exponential distribution given a mean
double random_exp(double mean);

// Update histogram
int update_hist(histogram_t *hist, double val);

// Update multidimensional histogram
int update_hist_3d(histogram_3d_t *hist, double *vals);

// Convert the components (module) of a polarization vector on the basis C to a basis D
double *update_polarization_vector(double *eps, double **r3_C, double **r3_D);

// Generate a orthonormal basis given a vector
double **orthonormal_basis(double *v3);

// Pick randomly a element of an integer array given an array of probabilities
int random_pick(int *arr, double *probs, int size);

// Replace a character in a string
char* str_replace(char *orig, char *rep, char *with);