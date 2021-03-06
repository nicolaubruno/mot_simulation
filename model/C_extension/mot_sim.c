
//  Header

#include "mot_sim.h"

// --

// Execute simulation for a single atom
int simulate_atom(){
    // Seed used by the rand() function
    srand(time(NULL));

    //
    // Variables
    //

    int i;              // Iterations
    float r, v, dt;     // Dynamics
    get_constants();    // Get constants variables from the CSV files

    // Parameters of the simulation
    conditions_t conds = get_conditions();
    beams_setup_t beams_setup = get_beams();
    atom_t atom = get_atom(conds);

    // Show all parameters
    //show_all_parameters(atom, conds, beams_setup);

    // Position histogram
    histogram_t pos_hist[3];
    for(i = 0; i < 3; i++){
        // Histogram
        pos_hist[i].num_bins = conds.num_bins;
        pos_hist[i].bin_size = 2 * conds.r_max / pos_hist[i].num_bins;
        pos_hist[i].coord0 = - conds.r_max;
        pos_hist[i].freqs = (int *) calloc(pos_hist[i].num_bins, sizeof(int));

        // Insert initial position
        update_hist(&pos_hist[i], atom.pos[i]);
    }

    // Results
    results_t res;
    res.pos_hist = pos_hist;
    res.num_iters = 0;
    res.flag = 2;

    // Compute distance and speed of the atom
    r = sqrt(r3_inner_product(atom.pos, atom.pos)); // Distance from centre
    v = sqrt(r3_inner_product(atom.vel, atom.vel)); // Total speed

    //  Iterations
    res.num_iters = conds.i_max - 1; 
    while((res.num_iters < conds.i_max) && (r < conds.r_max)){
        // Compute the photonic recoil
        dt = compute_photonic_recoil(atom, beams_setup, conds);

        // Compute gravitational force
        if(conds.g_bool) compute_gravitational_force(atom);

        // Compute magnetic force
        compute_magnetic_force(atom, conds.B_0);

        res.num_iters++;
    }

    return 0;
}

// Get atom
atom_t get_atom(conditions_t conds){
    //
    // Variables
    //

    atom_t atom;

    int row_cter = 0;
    int i;

    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/atom.csv");
    char *token, *rest;
    
    float std_dev;

    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    fgets(row, STRING_BUFFER_SIZE, fp);

    // Symbol
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"symbol\" in the file \"%s\"\n", path);
            exit(0);
        } else {
            if(token[2] == '\n') token[2] = '\0';

            atom.symbol = (char *) malloc(strlen(token) * sizeof(char));
            strcpy(atom.symbol, token);
        };
    }

    // Atomic number
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"Z\" in the file \"%s\"\n", path);
            exit(0);
        } else atom.Z = atoi(token);
    }

    // Mass
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"mass\" in the file \"%s\"\n", path);
            exit(0);
        } else atom.mass = atoi(token);
    }

    //
    // Initial position and velocity
    //

    atom.pos = (float *) calloc(3, sizeof(float));
    atom.vel = (float *) calloc(3, sizeof(float));

    for(i = 0; i < 3; i++){
        // Position
        atom.pos[i] = 0; // cm

        // Velocity
        std_dev = sqrt(k_B * conds.T_0 / (atom.mass * u)) * 10; // cm / s
        atom.vel[i] = norm(0, std_dev); // cm / s
    }

    // Optical transition
    atom.transition = get_transition();

    fclose(fp);
    return atom;
}

// Get transition
transition_t get_transition(){
    //
    // Variables
    //

    transition_t transition;
    int row_cter = 0;
    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/transition.csv");
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    fgets(row, STRING_BUFFER_SIZE, fp);

    // Transition rate
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"gamma\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.gamma = atof(token);
    }

    // Resonant wave length
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"lambda\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.lambda = atof(token);
    }

    // Total angular momentum of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"J_gnd\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.J_gnd = atoi(token);
    }

    // Total angular momentum of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"J_exc\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.J_exc = atoi(token);
    }

    // Landè factor of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_gnd\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.g_gnd = atof(token);
    }

    // Landè factor of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_exc\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.g_exc = atof(token);
    }

    //
    // Check values
    //

    if((transition.J_exc - transition.J_gnd) < 0){
        printf("J_exc must be grater than J_gnd.\n");
        exit(0);
    } else if(transition.J_exc < 0 || transition.J_gnd < 0){        
        printf("J_exc and J_gnd must be positive values.\n");
        exit(0);
    }

    fclose(fp);
    return transition;
}

// Get conditions
conditions_t get_conditions(){
    //
    // Variables
    //

    conditions_t conditions;
    int row_cter = 0;
    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/conditions.csv");
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    fgets(row, STRING_BUFFER_SIZE, fp);

    // Initial temperature 
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"T_0\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.T_0 = atof(token);
    }

    // Magnetic field gradient
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"B_0\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.B_0 = atof(token);
    }

    // Gravity
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_bool\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.g_bool = atoi(token);
    }

    // Maximum number of iteration
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"i_max\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.i_max = atoi(token);
    }

    // Maximum distance
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"r_max\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.r_max = atof(token);
    }

    // Number of simulations
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"num_sim\" in the file \"%s\"\n", path);
            exit(0);
        }
    }

    // Number of bins
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"num_bins\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.num_bins = atoi(token);
    }

    fclose(fp);
    return conditions;
}

// Get beams setup
beams_setup_t get_beams(){
    //
    // Variables
    //

    beams_setup_t beams_setup;
    beam_t *beams = (beam_t*) malloc(MAX_BEAMS * sizeof(beam_t));
    beam_t *c_beams;

    int row_cter = 0, num_beams = 0, n;
    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/beams.csv");
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    fgets(row, STRING_BUFFER_SIZE, fp);

    // Get beams
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        //
        // Laser detuning
        //

        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"delta\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].delta = atof(token);

        //
        // Wave vector direction
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"k_dic\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].k_dic = r3_normalize(get_float_array(token, &n));

        //
        // Polarization vector
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"eps\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].eps = r3_normalize(get_float_array(token, &n));

        //
        // Peak of the saturation parameter
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"s_0\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].s_0 = atof(token);

        //
        // Waist Radius
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"w\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].w = atof(token);

        num_beams++;
    }

    c_beams = (beam_t *) malloc(num_beams * sizeof(beam_t));
    for(n = 0; n < num_beams; n++) c_beams[n] = beams[n];

    beams_setup.num = num_beams;
    beams_setup.beams = c_beams;

    fclose(fp);
    return beams_setup;
}

// Get physical constant from CSV files
int get_constants(){
    //
    // Variables
    //

    int row_cter = 0;
    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/constants.csv");
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    fgets(row, STRING_BUFFER_SIZE, fp);

    // Planck constant
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"h\" in the file \"%s\"\n", path);
            exit(0);
        } else h = atof(token);
    }

    // Elementary charge
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"e\" in the file \"%s\"\n", path);
            exit(0);
        } else e = atof(token);
    }

    // Speed of light
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"c\" in the file \"%s\"\n", path);
            exit(0);
        } else c = atof(token);
    }

    // Boltzmann constant
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"k_B\" in the file \"%s\"\n", path);
            exit(0);
        } else k_B = atof(token);
    }

    // Bohr magneton
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"mu_B\" in the file \"%s\"\n", path);
            exit(0);
        } else mu_B = atof(token);
    }

    // Atomic mass unit
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"u\" in the file \"%s\"\n", path);
            exit(0);
        } else u = atof(token);
    }

    fclose(fp);
    return 1;
}

// Apply the photonic recoil in the atom returning the time interval of the process
float compute_photonic_recoil(atom_t atom, beams_setup_t beams_setup, conditions_t conds){
    //
    // Variables
    //

    int i, j, k;                // Iterations variables
    float aux_dt, dt = 0;       // Times variables
    float *B, *eps;             // Electromagnetic fields

    // Get Magnetic field
    B = get_magnetic_field(conds.B_0, atom.pos);

    //
    // Check transitions
    //

    // Loop each beam
    beams_setup.num = 1;
    for(i = 0; i < beams_setup.num; i++){
        B[0] = 0;
        B[1] = 0;
        B[2] = 1;
        
        // Get a R3 vector whose components are the module of the polarization vector components on the B frame
        eps = update_polarization_vector(beams_setup.beams[i], B);
    }


    return dt;
}

// Compute the momentum due to the gravitational force
int compute_gravitational_force(atom_t atom){
    return 1;
}

// Compute the momentum due the magnetic force
int compute_magnetic_force(atom_t atom, float B_0){
    return 1;
}

// Get magnetic field vector on the lab frame in the position r = (x, y, z)
float *get_magnetic_field(float B_0, float *r){
    //
    // Variables
    //
    static float B[3];      // Magnetic field vector

    B[0] = B_0 * r[0];
    B[1] = B_0 * r[1];
    B[2] = - 2 * B_0 * r[2];

    return B;
}

// Get coordinates of a polarization vector on the basis with pi transition parallel to the magnetic field B
/*
float *get_polarization_vector(beam_t beam, float *B){
    //
    // Variables
    //

    static float eps[3] = {0, 0, 0};    // Desired polarization vector  
    complex_t B[3][3];                  // Polarization basis of the beam frame on the lab frame
    complex_t C[3][3];                  // Polarization basis of the magnetic field frame on the lab frame


    complex_t B1_eps[3];                // Polarization vector on B1 (modules)
    complex_t B2_eps[3];                // Polarization vector on B2 (modules)
    float **B1, **B2;                   // Basis
    complex_t z;
    int i, j;

    // Check magnetic field
    if(B[0] == 0.0 && B[1] == 0.0 && B[2] == 0.0) B[2] = 1;
    else B = r3_normalize(B);

    //
    // Get orthonormal basis
    //    
    
    B1 = get_orthonormal_basis(beam.k_dic);
    B2 = get_orthonormal_basis(B);

    printf("B1 = {\n");
    for(i = 0; i < 3; i++) printf("\t[%f, %f, %f]\n", B1[i][0], B1[i][1], B1[i][2]);
    printf("}\n");

    // Define polarization on basis B2
    z.re = (beam.eps[0] + beam.eps[2]) / sqrt(2);
    z.im = 0;

    B1_eps[0] = z;

    z.im = (beam.eps[0] - beam.eps[2]) / sqrt(2);
    z.re = 0;

    B1_eps[1] = z;

    z.re = beam.eps[1];
    z.im = 0;

    B1_eps[2] = z;

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            z.re = r3_inner_product(B1[j], B2[i]);
            z.im = 0;

            B2_eps[i] = c_sum(B2_eps[i], c_inner_product(B1_eps[j], z));
        }
    }

    z.re = 0;
    z.im = 1;

    eps[0] = c_mod(c_sum(B2_eps[0], B2_eps[2])) / sqrt(2);
    eps[1] = c_mod(B2_eps[1]);
    eps[2] = c_mod(c_product(z, c_diff(B2_eps[0], B2_eps[2]))) / sqrt(2);

    return eps;
}
*/

// Get coordinates of a polarization vector on the basis with pi transition parallel to the magnetic field B
/*
float *update_polarization_vector_v1(beam_t beam, float *B){
    //
    // Variables
    //

    int i, j;
    static float *eps;                  // Desired polarization vector  
    float **r3_C, **r3_D;               // R3 Bases
    complex_t **c3_A, **c3_C, **c3_D;   // C3 Bases
    complex_t *C_eps, *A_eps, *D_eps;   // Complex polarization vector

    //
    // Define basis 
    //

    // Real orthonormal bases
    r3_C = orthonormal_basis(beam.k_dic);
    r3_D = orthonormal_basis(B);

    // Polarization basis of the magnetic field frame
    c3_A = (complex_t **) malloc(3 * sizeof(complex_t *));
    c3_C = (complex_t **) malloc(3 * sizeof(complex_t *));
    c3_D = (complex_t **) malloc(3 * sizeof(complex_t *));

    // Allocate memory
    for(i = 0; i < 3; i++){
        c3_A[i] = (complex_t *) calloc(3, sizeof(complex_t));
        c3_C[i] = (complex_t *) calloc(3, sizeof(complex_t));
        c3_D[i] = (complex_t *) calloc(3, sizeof(complex_t));
    }

    // Components on the lab frame
    for(i = 0; i < 3; i++){
        // sigma+
        c3_D[0][i].re = r3_D[0][i] / sqrt(2);
        c3_D[0][i].im = r3_D[1][i] / sqrt(2);

        c3_C[0][i].re = r3_C[0][i] / sqrt(2);
        c3_C[0][i].im = r3_C[1][i] / sqrt(2);

        c3_A[0][i].re = 0;
        c3_A[0][i].im = 0;

        // sigma-
        c3_D[1][i].re = r3_D[0][i] / sqrt(2);
        c3_D[1][i].im = -r3_D[1][i] / sqrt(2);

        c3_C[1][i].re = r3_C[0][i] / sqrt(2);
        c3_C[1][i].im = -r3_C[1][i] / sqrt(2);

        c3_A[1][i].re = 0;
        c3_A[1][i].im = 0;

        // pi
        c3_D[2][i].re = r3_D[2][i];
        c3_D[2][i].im = 0;

        c3_C[2][i].re = r3_C[2][i];
        c3_C[2][i].im = 0;

        c3_A[2][i].re = 0;
        c3_A[2][i].im = 0;
    }

    // Lab frame
    c3_A[0][0].re = 1;
    c3_A[1][1].re = 1;
    c3_A[2][2].re = 1;

    //
    // Compute polarization on different bases
    //

    C_eps = (complex_t *) calloc(3, sizeof(complex_t));
    for(i = 0; i < 3; i++){
        C_eps[i].re = beam.eps[i];
        C_eps[i].im = 0;
    }

    c3_view(C_eps, "C_eps");

    D_eps = (complex_t *) calloc(3, sizeof(complex_t));
    for(i = 0; i < 3; i++) {
        D_eps[i].re = 0;
        D_eps[i].im = 0;

        for(j = 0; j < 3; j++){
            D_eps[i] = c_sum(D_eps[i], c_product(C_eps[j], c3_inner_product(c3_D[i], c3_C[j])));
        }
    }

    c3_view(D_eps, "D_eps");

    eps = (float*) calloc(3, sizeof(float));
    for(i = 0; i < 3; i++) eps[i] = c_mod(D_eps[i]); 
    eps = r3_normalize(eps);

    // Return
    return eps;
}
*/

// Get a R3 vector whose components are the module of the polarization vector components on the B frame
float *update_polarization_vector(beam_t beam, float *B){
    //
    // Variables
    //

    int i, j;
    float *eps;                         // Desired vector
    complex_t *C_eps, *D_eps;           // Polarization vector on different Cartesian bases    
    complex_t *Dp_eps, *Cp_eps, *aux;   // Polarization vector on different polarization bases
    complex_t **A1, **A1_i, **A2;       // Change-of-basis matrices
    float **r3_C, **r3_D;               // Real Bases
    complex_t **c3_C, **c3_D;           // Complex Bases
    complex_t im_p;

    // Imaginary particle
    im_p.re = 0;
    im_p.im = 1;

    //
    // Change-of-basis matrix from the polarization basis to the Cartesian basis
    //

    A1 = (complex_t **) calloc(3, sizeof(complex_t*));
    for(i = 0; i < 3; i++){ 
        A1[i] = (complex_t *) calloc(3, sizeof(complex_t)); 

        for(j = 0; j < 3; j++){
            A1[i][j].re = 0;
            A1[i][j].im = 0;
        }
    }

    A1[0][0].re = 1 / sqrt(2);
    A1[0][1].re = 1 / sqrt(2);
    A1[1][0].im = 1 / sqrt(2);
    A1[1][1].im = - 1 / sqrt(2);
    A1[2][2].re = 1;
    
    
    //
    // Change-of-basis matrix from the Cartesian basis to the polarization basis
    //

    A1_i = (complex_t **) malloc(3 * sizeof(complex_t*));
    for(i = 0; i < 3; i++){ 
        A1_i[i] = (complex_t *) calloc(3, sizeof(complex_t)); 

        for(j = 0; j < 3; j++){
            A1_i[i][j].re = A1[j][i].re;
            A1_i[i][j].im = -A1[j][i].im;
        }
    }

    //
    // Change-of-basis matrix from the Cartesian beam frame to the Cartesian B frame
    //

    // Real orthonormal bases
    r3_C = orthonormal_basis(beam.k_dic);
    r3_D = orthonormal_basis(B);

    c3_C = (complex_t**) calloc(3, sizeof(complex_t *));
    c3_D = (complex_t**) calloc(3, sizeof(complex_t *));
    
    for(i = 0; i < 3; i++){
        c3_C[i] = r3_to_c3(r3_C[i]);
        c3_D[i] = r3_to_c3(r3_D[i]);
    }

    A2 = (complex_t **) malloc(3 * sizeof(complex_t*));
    for(i = 0; i < 3; i++){ 
        A2[i] = (complex_t *) calloc(3, sizeof(complex_t)); 

        for(j = 0; j < 3; j++){
            A2[i][j] = c3_inner_product(c3_D[i], c3_C[j]);
        }
    }

    //
    // Polarization vector on the polarization beam frame
    //

    Cp_eps = r3_to_c3(beam.eps);
    C_eps = c3_apply_operator(A1, Cp_eps);

    // Polarization vector on the Cartesian B frame
    D_eps = c3_apply_operator(A2, C_eps);

    // Polarization vector on the polarization B frame
    Dp_eps = c3_apply_operator(A1_i, D_eps);

    // Compute desired vector
    eps = (float*) calloc(3, sizeof(float));
    for(i = 0; i < 3; i++) eps[i] = c_mod(Dp_eps[i]);

    return eps;
}

//
// Utility functions
//

// Get int array (max 124 numbers) from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size){
    //
    // Variables
    //
    int i, j, max_size = 124;
    char *token;
    int aux_arr[124];
    static int *arr;

    str = str + 1;
    str[strlen(str)-1] = '\0';
    token = strtok_r(str, " ", &str);    

    //
    // Parse string
    //

    i = 0;
    while(token && (i < max_size)){
        if(token[0] == '-') aux_arr[i] = -atoi((token + 1)); 
        else aux_arr[i] = atoi(token);

        token = strtok_r(str, " ", &str);
        i++;
    }

    *size = i;

    arr = (int*) malloc(i * sizeof(int));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    return arr;
}

// Get float array from string in the format [f1 f2 ... fn]
float *get_float_array(char *str, int *size){
    //
    // Variables
    //
    int i, j, max_size = 124;
    char *token;
    float aux_arr[max_size];
    static float *arr;

    str = str + 1;
    str[strlen(str)-1] = '\0';
    token = strtok_r(str, " ", &str);    

    //
    // Parse string
    //

    i = 0;
    while(token && (i < max_size)){
        if(token[0] == '-') aux_arr[i] = -atof((token + 1)); 
        else aux_arr[i] = atof(token);

        token = strtok_r(str, " ", &str);
        i++;
    }

    arr = (float *) malloc(i * sizeof(float));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    return arr;
}

// Concatenate ROOT_PATH to a filename
char *concatenate_ROOT_PATH(char *filename){
    int i;
    char aux_path[124];
    static char *path;

    aux_path[0] = '\0';

    strcat(aux_path, ROOT_PATH);
    strcat(aux_path, filename);

    i = 1;
    while((aux_path[i-1] != '\0') && (i < 124)) i++;

    path = (char*) malloc(i*sizeof(char));
    strcpy(path, aux_path);

    return path;
}

// Print all parameters of the simulation
int show_all_parameters(atom_t atom, conditions_t conditions, beams_setup_t beams_setup){
    //
    // Variables
    //

    int i;

    //
    // Prints
    //

    printf("\nAtom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass = %f\n", atom.mass);

    printf("\nTransition\n--\n");
    printf("gamma = %f\n", atom.transition.gamma);
    printf("J_gnd = %d\n", atom.transition.J_gnd);
    printf("J_exc = %d\n", atom.transition.J_exc);
    printf("g_gnd = %f\n", atom.transition.g_gnd);
    printf("g_exc = %f\n", atom.transition.g_exc);

    printf("\nConditions\n--\n");
    printf("T_0 = %f\n", conditions.T_0);
    printf("B_0 = %f\n", conditions.B_0);
    printf("g_bool = %d\n", conditions.g_bool);
    printf("i_max = %d\n", conditions.i_max);
    printf("r_max = %f\n", conditions.r_max);
    printf("num_bins = %d\n", conditions.num_bins);

    printf("\nBeams setup\n--\n");
    printf("num_beams = %d\n", beams_setup.num);

    for(i = 0; i < beams_setup.num; i++){
        printf("\nBeam %d\n--\n", (i+1));
        printf("delta = %f\n", beams_setup.beams[i].delta);
        printf("k_dic = [%f, %f, %f]\n", beams_setup.beams[i].k_dic[0], beams_setup.beams[i].k_dic[1], beams_setup.beams[i].k_dic[2]);
        printf("eps = [%f, %f, %f]\n", beams_setup.beams[i].eps[0], beams_setup.beams[i].eps[1], beams_setup.beams[i].eps[2]);
        printf("s_0 = %f\n", beams_setup.beams[i].s_0);
        printf("w = %f\n", beams_setup.beams[i].w);
    }
}

// Generate a float random number following a Gaussian distribution given a mean and a standard deviation
float norm(float mean, float std_dev){
    //
    // Variables
    //

    int i;
    float u[2], r, theta;   // Variables for Box-Muller method
    float std_norm;         // Normal(0, 1)
    float norm;             // Adjusted normal

    //
    // Box-Muller transform
    //

    // Generate uniform random numbers
    for(i = 0; i < 2; i++) u[i] = ((float) rand()) / ((float) RAND_MAX);

    // Compute r
    r = sqrt(-2 * log(u[0]));

    // Generate theta
    theta = 0.0;
    while(theta == 0.0) theta = 2.0 * PI * u[1];

    // Generate std_norm value
    std_norm = r * cos(theta);

    // Adjust std_norm
    norm = (std_norm * std_dev) + mean;

    return norm;
}

// Update histogram
int update_hist(histogram_t *hist, float val){
    //
    // Variables
    //

    int bin;
    float lower_lim, upper_lim;

    // Add frequency
    for(bin = 0; bin < (*hist).num_bins; bin++){
        lower_lim = (*hist).coord0 + bin * (*hist).bin_size;
        upper_lim = lower_lim + (*hist).bin_size;

        if((val >= lower_lim) && (val <= upper_lim)){
            (*hist).freqs[bin]++;
            break;
        }
    }

    return 1;
}

// Generate a orthonormal basis given a vector
float **orthonormal_basis(float *v3){
    //
    // Variables
    //

    int i;
    static float **B;   // Desired basis
    float *v1, *v2;     // Auxiliary vectors

    // Normalize vector v
    v3 = r3_normalize(v3);

    // Generate a random vector  
    v1 = (float*) calloc(3, sizeof(float));
    for(i = 0; i < 3; i++) v1[i] = ((float) rand()) / ((float) RAND_MAX);

    // Define a orthonormal vector
    v2 = r3_scalar_product(r3_inner_product(v1, v3), v3);
    v1 = r3_normalize(r3_diff(v1, v2));

    // Define the last vector of the basis
    v2 = r3_cross_product(v3, v1);

    // Define basis
    B = (float **) malloc(3 * sizeof(float *));
    
    B[0] = v1;
    B[1] = v2;
    B[2] = v3;

    return B;
}