//
// Libraries
//

#include "vectors.h"

//
// Complex space
//--

complex_t c_sum(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re + z2.re;
    z.im = z1.im + z2.im;

    return z;
}

complex_t c_diff(complex_t z1, complex_t z2){
    complex_t z;

    z.re = -z2.re;
    z.im = -z2.im;

    return c_sum(z1, z);
}

complex_t c_product(complex_t z1, complex_t z2){
    complex_t z;

    z.re = z1.re * z2.re - z1.im * z2.im;
    z.im = z1.re * z2.im + z1.im * z2.re;

    return z; 
}

complex_t c_conj(complex_t z){
    z.im = -z.im;
    return z;
}

complex_t c_inner_product(complex_t z1, complex_t z2){
    return c_product(z1, c_conj(z2));
}

double c_mod(complex_t z){
    return sqrt(c_inner_product(z, z).re);
}

//--

//
// C3 space
//--

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

complex_t *c3_scalar_product(complex_t a, complex_t *v){
    // Variables
    int i;
    complex_t *w;

    w = (complex_t*) calloc(3, sizeof(complex_t));

    for(i = 0; i < 3; i++)
        w[i] = c_product(a, v[i]);

    return w;
}

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

double c3_mod(complex_t *z){
    return sqrt(c3_inner_product(z, z).re);
}

int c3_print(complex_t *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("(%f, %f), ", z[i].re, z[i].im);
    printf("(%f, %f)]\n", z[i].re, z[i].im);

    return 1;
}

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

complex_t **r3_oper_to_c3_oper(double **A){
    int i, j;
    complex_t **B;

    B = (complex_t**) calloc(3, sizeof(complex_t*));

    for(i = 0; i < 3; i++){
        B[i] = (complex_t*) calloc(3, sizeof(complex_t));

        for(j = 0; j < 3; j++){
            B[i][j].re = A[i][j];
            B[i][j].im = 0;
        }
    }

    return B;
}

complex_t **c3_operator_zeros(){
    // Variables
    int i, j;
    complex_t **A;

    A = (complex_t**) calloc(3, sizeof(complex_t*));

    for(i = 0; i < 3; i++){
        A[i] = (complex_t*) calloc(3, sizeof(complex_t));

        for(j = 0; j < 3; j++){
            A[i][j].re = 0;
            A[i][j].im = 0;
        } 
    }

    return A;
}

// Conjugate transpose
complex_t **c3_operator_dagger(complex_t **A){
    // Variables
    int i, j;
    complex_t **A_dagger;

    A_dagger = (complex_t**) calloc(3, sizeof(complex_t*));

    for(i = 0; i < 3; i++){
        A_dagger[i] = (complex_t*) calloc(3, sizeof(complex_t));

        for(j = 0; j < 3; j++){
            // Transpose
            A_dagger[i][j] = A[j][i];

            // Conjugate
            A_dagger[i][j].im = - A_dagger[i][j].im;
        } 
    }

    return A_dagger;
}

//--

//
// R3 space
//--

double r3_mod(double *z){
    return sqrt(r3_inner_product(z, z));
}

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

double *r3_scalar_product(double a, double *v){
    // Variables
    int i;
    double *w;

    w = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++)
        w[i] = a * v[i];

    return w;
}

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

    mod_r = r3_mod(r);
    if(mod_r == 0.0){
        printf("Division by zero in r3_normalize() function\n");
        exit(0);
    }
    for(i = 0; i < 3; i++) new_r[i] = (r[i] / mod_r);

    return new_r;
}

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

int r3_print(double *z, char *name){
    int i;

    printf("%s = [", name);
    for(i = 0; i < 2; i++) printf("%f, ", z[i]);
    printf("%f]\n", z[i]);

    return 1;
}

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

//--