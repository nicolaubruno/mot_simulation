//
// Libraries
//

#include "vectors.h"

//
// Complex space
//

// Sum
complex_t c_sum(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re + z2.re;
    z.im = z1.im + z2.im;

    return z;
}

// Difference
complex_t c_diff(complex_t z1, complex_t z2){
    complex_t z;

    z.re = -z2.re;
    z.im = -z2.im;

    return c_sum(z1, z);
}

// Product
complex_t c_product(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re * z2.re - z1.im * z2.im;
    z.im = z1.re * z2.im + z1.im * z2.re;

    return z; 
}

// Conjugate
complex_t c_conj(complex_t z){
    z.im = -z.im;
    return z;
}

// Inner product
complex_t c_inner_product(complex_t z1, complex_t z2){
    return c_product(z1, c_conj(z2));
}

// Module
double c_mod(complex_t z){
    return sqrt(c_inner_product(z, z).re);
}

// Sum between C3 vectors
complex_t *c3_sum(complex_t *z1, complex_t *z2){
    //
    // Variables
    //

    int i;
    complex_t *res;

    res = (complex_t*) calloc(3, sizeof(complex_t));

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = c_sum(z1[i], z2[i]);

    return res;
}

// Difference between C3 vectors
complex_t *c3_diff(complex_t *z1, complex_t *z2){
    //
    // Variables
    //

    int i;
    complex_t *res;

    res = (complex_t*) calloc(3, sizeof(complex_t));

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = c_diff(z1[i], z2[i]);

    return res;
}

// Product between a complex number and a real vector
complex_t *c3_scalar_product(complex_t a, complex_t *v){
    // Variables
    int i;
    complex_t *w;

    w = (complex_t*) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++)
        w[i] = c_product(a, v[i]);

    return w;
}

// Inner product between C3 vectors
complex_t c3_inner_product(complex_t *z1, complex_t *z2){
    int i;
    complex_t z;

    z.re = 0;
    z.im = 0;

    for(i = 0; i < 3; i++){
        z = c_sum(z, c_inner_product(z1[i], z2[i]));
    }

    return z;
}

// Module
double c3_mod(complex_t *z){
    return sqrt(c3_inner_product(z, z).re);
}

// Print complex vector
int c3_print(complex_t *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("(%f, %f), ", z[i].re, z[i].im);
    printf("(%f, %f)]\n", z[i].re, z[i].im);

    return 1;
}

// Apply complex operator
complex_t *c3_apply_operator(complex_t **A, complex_t *v){
    int i,  j;
    complex_t *a;

    a = (complex_t*) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++){
        a[i].re = 0;
        a[i].im = 0;

        for(j = 0; j < 3; j++){
            a[i] = c_sum(a[i], c_product(A[i][j], v[j]));
        }
    }

    return a;
}

// Print C3 operator
int c3_operator_print(complex_t **A, char *name){
    int i;

    printf("%s = [\n", name);
    for(i = 0; i < 3; i++) {
        printf("\t(%f, %f) ", A[i][0].re, A[i][0].im);
        printf("(%f, %f) ", A[i][1].re, A[i][1].im);
        printf("(%f, %f)\n", A[i][2].re, A[i][2].im);
    }
    printf("]\n");

    return 1;
}

// Convert a R3 vector to a C3 vector
complex_t *r3_to_c3(double *v){
    int i;
    complex_t *res;

    res = (complex_t *) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++){
        res[i].re = v[i];
        res[i].im = 0;
    }

    return res;
}

//
// Real space
//

// Module
double r3_mod(double *z){
    return sqrt(r3_inner_product(z, z));
}

// Inner product between R3 vectors
double r3_inner_product(double *r1, double *r2){
    //
    // Variables
    //

    int i;
    double res = 0;

    // Compute dot product
    for(i = 0; i < 3; i++) res += r1[i] * r2[i];

    return res;
}

// Cross product between R3 vectors
double *r3_cross_product(double *r1, double *r2){
    //
    // Variables
    //

    double *res;

    res = (double*) calloc(3, sizeof(double));

    // Compute dot product
    res[0] = r1[1]*r2[2] - r1[2]*r2[1];
    res[1] = r1[2]*r2[0] - r1[0]*r2[2];
    res[2] = r1[0]*r2[1] - r1[1]*r2[0];

    return res;
}

// Product between a real number and a real vector
double *r3_scalar_product(double a, double *v){
    // Variables
    int i;
    double *w;

    w = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++)
        w[i] = a * v[i];

    return w;
}

// Sum between R3 vectors
double *r3_sum(double *r1, double *r2){
    //
    // Variables
    //

    int i;
    double *res;

    res = (double*) calloc(3, sizeof(double));

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = r1[i] + r2[i];

    return res;
}

// Subtraction between R3 vectors
double *r3_diff(double *r1, double *r2){
    //
    // Variables
    //

    int i;
    double *res;

    res = (double*) calloc(3, sizeof(double));

    // Compute dot product
    for(i = 0; i < 3; i++) res[i] = r1[i] - r2[i];

    return res;
}

// Normalize a R3 vector
double *r3_normalize(double *r){
    //
    // Variables
    //

    double *new_r;
    double mod_r;
    int i;

    new_r = (double*) calloc(3, sizeof(double));

    //
    // Normalization
    //

    mod_r = sqrt(r3_inner_product(r, r));
    if(mod_r == 0.0){
        printf("Division by zero in normalize_vector() function\n");
        exit(0);
    }
    for(i = 0; i < 3; i++) new_r[i] = (r[i] / mod_r);

    return new_r;
}

// Apply real operator
double *r3_apply_operator(double **A, double *v){
    int i,  j;
    double *a;

    a = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            a[i] += A[i][j]*v[j];
        }
    }

    return a;
}

// Print R3 vector
int r3_print(double *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("%f, ", z[i]);
    printf("%f]\n", z[i]);

    return 1;
}

// Print R3 operator
int r3_operator_print(double **A, char *name){
    int i, j;

    printf("%s \n--\n", name);
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 2; j++) printf("%f ", A[i][j]);
        printf("%f\n", A[i][2]);
    }
    printf("\n");

    return 1;
}