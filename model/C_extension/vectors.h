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

//
// C3 Space
//

// Sum
complex_t *c3_sum(complex_t *z1, complex_t *z2);

// Difference
complex_t *c3_diff(complex_t *z1,complex_t *z2);

// Product between a complex number and a C3 vector
complex_t *c3_scalar_product(complex_t a, complex_t *v);

// Inner product
complex_t c3_inner_product(complex_t *z1, complex_t *z2);

// Module
double c3_mod(complex_t *z);

// Print C3 vector
int c3_print(complex_t *z, char *name);

// Apply complex operator to a c3 vector
complex_t *c3_apply_operator(complex_t **A, complex_t *v);

// Print C3 operator
int c3_operator_print(complex_t **A, char *name);

// Convert a R3 vector to a C3 vector
complex_t *r3_to_c3(double *v);

// Convert a R3 operator to a C3 operator
complex_t **r3_oper_to_c3_oper(double **A);

// Empty C3 operator
complex_t **c3_operator_zeros();

// Complex conjugate operator
complex_t **c3_operator_dagger(complex_t **A);

//
// R3 space
//

// Module
double r3_mod(double *z);

// Inner product
double r3_inner_product(double *r1, double *r2);

// Cross product
double *r3_cross_product(double *r1, double *r2);

// Product between a real number and a R3 vector
double *r3_scalar_product(double a, double *v);

// Sum
double *r3_sum(double *r1, double *r2);

// Difference
double *r3_diff(double *r1, double *r2);

// Normalize
double *r3_normalize(double *r);

// Apply real operator to a R3 vector
double *r3_apply_operator(double **A, double *v);

// View R3 vector
int r3_print(double *z, char *name);

// Print R3 operator
int r3_operator_print(double **A, char *name);
