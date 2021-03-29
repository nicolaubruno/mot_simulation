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
#define MAX_LINES 124
#define DELIM ","
#define MAX_BEAMS 16
#define Py_MODULE 1

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
    int J;                      /* Total angular momentum */
    int mJ;                     /* Magnetic quantum number */
    transition_t transition;    /* Optical transition */
} atom_t;

// Conditions
typedef struct{
    double T_0;         /* Initial temperature [uK] */
    double max_time;    /* Maximum time of simulation [1/gamma] */
    double max_r;       /* Maximum distance (threshold) [cm] */
    double max_v;       /* Maximum speed [cm/s] */
    int num_bins;       /* Number of bins in each histogram */
    double wait_time;   /* Time to reach the equilibrium [1/gamma] */
    double dt;          /* Time interval [1/gamma] */
} conditions_t;

// Environment
typedef struct{
    double B_0;         /* Magnetic Field gradient */
    double **B_basis;   /* Basis with z-axis parallel to the axial direction of the magnetic field */
    double local_B;     /* Local magnetic field gradient on the z direction */
    double delta;       /* Laser detuning */
    double s_0;         /* Peak of the saturation parameter */
    double w;           /* Waist radius */
    int g_bool;         /* Use gravity */
} environment_t;

// Beam
typedef struct {
    double s_0;      /* Resonant saturation parameter */
    double *k_dic;   /* Wave vector direction */
    double *eps;     /* Polarization vector */
} beam_t;

// Polarized Beam
typedef struct {
    double s_0;      /* Resonant saturation parameter */
    double *k_dic;   /* Wave vector direction */
    int eps;         /* Polarization (+1, -1, 0) */
} polarized_beam_t;

// Beams setup
typedef struct {
    int num;        /* Number of beams */
    beam_t *beams;  /* All beams */
} beams_setup_t;

// Photon
typedef struct{
    int m;           /* Magnetic quantum number */
    double *eK;         /* Wave vector direction */
    double lambda;      /* Wavelength [nm] */
} photon_t;

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

// Results
typedef struct{
    histogram_3d_t pos_3Dhist;      /* 3D-Histogram of position */
    histogram_t *pos_hist;          /* Marginal positions histograms */
    histogram_3d_t vel_3Dhist;      /* 3D-Histogram of velocity */
    histogram_t *vel_hist;          /* Marginal velocities histograms */
    histogram_t speed_hist;        /* Speed (velocity module) histogram */
    double time;                    /* Total time [s] */
    int *transitions;               /* Counter of occurred transitions */
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
results_t simulate_atom(char *params_path, int marginals, long seed_time);

// Get atom
atom_t get_atom(conditions_t conds, environment_t env, char *params_path);

// Get transition
transition_t get_transition(char *params_path);

// Get conditions
conditions_t get_conditions(char *params_path);

// Get environment
environment_t get_environment(char *params_path);

// Get beams setup
beams_setup_t get_beams(char *params_path);

// Set histograms
int set_hist(int only_marginals, results_t *res, conditions_t conds);

// Get polarizations amplitudes
double *polarizations_amplitudes(beam_t beam, environment_t env, double *eB);

// Apply movement on the atom due to photonic recoil, gravitational force and magnetic force
double move(atom_t *atom, beams_setup_t beams_setup, conditions_t conds, environment_t env);

// Photon absorption event
double photon_absorption(atom_t *atom, beams_setup_t beams_setup, conditions_t conds, environment_t env, double dt);

// Photon emission event
double photon_emission(atom_t *atom, environment_t env);

// Get magnetic acceleration
double *magnetic_acceleration(atom_t atom, environment_t env);

// Get magnetic field vector in the lab frame
double *magnetic_field(environment_t env, double *r);

// Get scattering rate
double scattering_rate(atom_t atom, polarized_beam_t beam, conditions_t conds, environment_t env, double *B);

// Get components of a vector v on the basis B given the components on basis A
double *change_basis(double *v, double **A, double **B_basis);

//
// Utility functions
//

// (Debug) print parameters
int print_params(atom_t atom, conditions_t conds, beams_setup_t beams, environment_t env);

// (Debug) Print simulation status
int print_status(atom_t atom, results_t res);

// (Debug) Print results
int print_results(results_t res, atom_t atom, int only_marginals);

// Get int array from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size);

// Get double array from string in the format [f1 f2 ... fn]
double *get_double_array(char *str, int *size);

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
int random_pick(double *probs, int size);

// Replace a character in a string
char* str_replace(char *orig, char *rep, char *with);

// Read lines from a file
char **read_lines(char *path);