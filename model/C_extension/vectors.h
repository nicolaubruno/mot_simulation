//
// Libraries
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//
// Structures
//

// Complex number
typedef struct {
    double re;    /* Real part */
    double im;    /* Imaginary part */
} complex_t;

//
// Complex Space
//

// Sum
complex_t c_sum(complex_t z1, complex_t z2);

// Difference
complex_t c_diff(complex_t z1, complex_t z2);

// Product
complex_t c_product(complex_t z1, complex_t z2);

// Conjugate
complex_t c_conj(complex_t z);

// Module
double c_mod(complex_t z);

// Inner product
complex_t c_inner_product(complex_t z1, complex_t z2);

// Sum between C3 vectors
complex_t *c3_sum(complex_t *z1, complex_t *z2);

// Difference between C3 vectors
complex_t *c3_diff(complex_t *z1,complex_t *z2);

// Product between a complex number and a real vector
complex_t *c3_scalar_product(complex_t a, complex_t *v);

// Inner product between C3 vectors
complex_t c3_inner_product(complex_t *z1, complex_t *z2);

// Module
double c3_mod(complex_t *z);

// View complex vector
int c3_print(complex_t *z, char *name);

// Apply complex operator
complex_t *c3_apply_operator(complex_t **A, complex_t *v);

// Print C3 operator
int c3_operator_print(complex_t **A, char *name);

// Convert a R3 vector to a C3 vector
complex_t *r3_to_c3(double *v);

//
// Real space
//

// Module
double r3_mod(double *z);

// Inner product between R3 vectors
double r3_inner_product(double *r1, double *r2);

// Cross product between R3 vectors
double *r3_cross_product(double *r1, double *r2);

// Product between a real number and a real vector
double *r3_scalar_product(double a, double *v);

// Sum between R3 vectors
double *r3_sum(double *r1, double *r2);

// Difference between R3 vectors
double *r3_diff(double *r1, double *r2);

// Normalize a R3 vector
double *r3_normalize(double *r);

// Apply real operator
double *r3_apply_operator(double **A, double *v);

// View real vector
int r3_print(double *z, char *name);

// Print R3 operator
int r3_operator_print(double **A, char *name);
