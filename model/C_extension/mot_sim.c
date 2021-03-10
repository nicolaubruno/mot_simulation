/**
 * Monte Carlos simulation of a single atom in a Magneto-Optical Trap (MOT)
 * Bruno N. Santos <nicolau.bruno@gmail.com>
 * Version 2.0
 */

//  Header

#include "mot_sim.h"

// --

results_t simulate_atom(){
    // Seed used by the rand() function
    srand(time(NULL));

    //
    // Variables
    //

    int i;                  // Iterations
    double r, dt;           // Dynamics
    results_t res;          // Results

    // Parameters of the simulation
    conditions_t conds = get_conditions();
    beams_setup_t beams_setup = get_beams();
    atom_t atom = get_atom(conds);

    // Scattering process
    scattering_t scatt;

    // Position histogram
    res.pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

    for(i = 0; i < 3; i++){
        // Histogram
        res.pos_hist[i].num_bins = conds.num_bins;
        res.pos_hist[i].bin_size = 2 * conds.r_max / res.pos_hist[i].num_bins;
        res.pos_hist[i].coord0 = - conds.r_max;
        res.pos_hist[i].freqs = (int *) calloc(res.pos_hist[i].num_bins, sizeof(int));

        // Insert initial position
        update_hist(&res.pos_hist[i], atom.pos[i]);
    }

    // Compute distance and speed of the atom
    r = sqrt(r3_inner_product(atom.pos, atom.pos)); // Distance from centre

    //
    //  Iterations
    //

    res.num_iters = 0; 
    dt = 0;
    res.time = 0;

    while((res.num_iters < conds.i_max) && (r < conds.r_max)){
        // Compute the photonic recoil
        scatt = photonic_recoil(atom, beams_setup, conds);

        // Compute gravitational force
        if(conds.g_bool) compute_gravitational_force(atom);

        // Compute magnetic force
        compute_magnetic_force(atom, conds.B_0);

        // Update time
        if(scatt.dt > 0 && scatt.dt < MAX_dt) dt = scatt.dt;
        else dt = 1 / atom.transition.gamma;
        res.time += dt;

        // Update position
        for(i = 0; i < 3; i++) atom.pos[i] += atom.vel[i] * dt;

        // Update velocity
        for(i = 0; i < 3; i++){
            if(scatt.dt > 0 && scatt.dt < MAX_dt) atom.vel[i] += scatt.vel[i];
        }

        // Update results
        for(i = 0; i < 3; i++) update_hist(&res.pos_hist[i], atom.pos[i]);

        r = r3_mod(atom.pos);
        res.num_iters++;

        // Print status
        //if(res.num_iters % 500 == 0) print_status(atom, res);
    }

    // Print status
    //print_status(atom, res);

    return res;
}

atom_t get_atom(conditions_t conds){
    //
    // Variables
    //

    atom_t atom;

    int i;

    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/atom.csv");
    char *token, *rest;
    
    double std_dev;

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

    atom.pos = (double *) calloc(3, sizeof(double));
    atom.vel = (double *) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++){
        // Position
        atom.pos[i] = 0; // cm

        // Velocity
        std_dev = sqrt(k_B * conds.T_0 / (atom.mass * u)) * 10; // cm / s
        atom.vel[i] = random_norm(0, std_dev); // cm / s
    }

    // Optical transition
    atom.transition = get_transition();

    // Close file
    fclose(fp);

    return atom;
}

transition_t get_transition(){
    //
    // Variables
    //

    transition_t transition;
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
        } else transition.gamma = atof(token) * 2 * PI;
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

    // Close file
    fclose(fp);

    // Release memory

    // Return
    return transition;
}

conditions_t get_conditions(){
    //
    // Variables
    //

    conditions_t conditions;
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
        } else conditions.i_max = (int) atof(token);
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

beams_setup_t get_beams(){
    //
    // Variables
    //

    beams_setup_t beams_setup;
    beam_t *beams = (beam_t*) malloc(MAX_BEAMS * sizeof(beam_t));
    beam_t *c_beams;

    int num_beams = 0, n;
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
        } else beams[num_beams].k_dic = r3_normalize(get_double_array(token, &n));

        //
        // Polarization vector
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"eps\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].eps = r3_normalize(get_double_array(token, &n));

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

    // Release memory
    free(beams);

    // Close file
    fclose(fp);

    // Return
    return beams_setup;
}

scattering_t photonic_recoil(atom_t atom, beams_setup_t beams_setup, conditions_t conds){
    //
    // Variables
    //

    int i, j;                                   // Iterations variables
    double *B, *eps_probs, *eB, *eK, lambda;    // Electromagnetic fields
    scattering_t *scatt_opt;                    // Scattering options
    double vel_mod, *rd_v;                      // Auxiliary variables
    int *pol_opt, pol;                          // Polarization

    // Get Magnetic field
    B = get_magnetic_field(conds.B_0, atom.pos);

    //
    // Get scattering length
    //

    scatt_opt = (scattering_t*) calloc(beams_setup.num, sizeof(scattering_t));

    // Loop each beam
    for(i = 0; i < beams_setup.num; i++){
        // Get wave length
        lambda = atom.transition.lambda; // nm

        // Directions
        if(r3_mod(B) == 0){
            eB = (double*) calloc(3, sizeof(double));
            eB[2] = 1.0;
        } else eB = r3_normalize(B);
        eK = r3_normalize(beams_setup.beams[i].k_dic);

        //
        // Polarization
        //

        pol_opt = (int*) calloc(3, sizeof(int));

        pol_opt[0] = +1;
        pol_opt[1] = -1;
        pol_opt[2] = 0;

        // Get polarization probabilities
        eps_probs = polarization_probs(beams_setup.beams[i], eB);
        pol = random_pick(pol_opt, eps_probs, 3);

        // Compute scattering
        scatt_opt[i].R = scattering_rate(atom, beams_setup.beams[i], B, pol); // Hz
        if(scatt_opt[i].R > 0) scatt_opt[i].dt = (1 / scatt_opt[i].R); // s
        else scatt_opt[i].dt = 0;

        // Absorption event
        vel_mod = 1e4 * h / (lambda * atom.mass * u); // cm / s
        scatt_opt[i].vel = r3_scalar_product(vel_mod, eK); // cm / s

        // Emission event
        rd_v = (double*) calloc(3, sizeof(double));  // Random vector
        for(j = 0; j < 3; j++) {
            rd_v[j] = ((double) rand()) / ((double) RAND_MAX);
            if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[j] = -rd_v[j];
        }
        rd_v = r3_normalize(rd_v); // Normalization
        rd_v = r3_scalar_product(vel_mod, rd_v); // cm / s
        scatt_opt[i].vel = r3_sum(scatt_opt[i].vel, rd_v); // cm / s
    }

    //
    // Pick scattering option
    //

    j = 0; 

    for(i = 0; i < beams_setup.num; i++){
        if(scatt_opt[i].dt > 0 && scatt_opt[i].dt < scatt_opt[j].dt)
            j = i; 
    }

    // Release memory
    free(eB);
    free(pol_opt);
    free(rd_v);

    // Return
    return scatt_opt[j];
}

int compute_gravitational_force(atom_t atom){
    return 1;
}

int compute_magnetic_force(atom_t atom, double B_0){
    return 1;
}

double *get_magnetic_field(double B_0, double *r){
    //
    // Variables
    //
    double *B;      // Magnetic field vector

    B = (double*) calloc(3, sizeof(double));

    B[0] = B_0 * r[0];
    B[1] = B_0 * r[1];
    B[2] = - 2 * B_0 * r[2];

    // Return
    return B;
}

double *polarization_probs(beam_t beam, double *eB){
    //
    // Variables
    //

    int i, j;
    double *eps_probs, *eK;             // Polarization probabilities and wave vector direction
    double **r3_C,  **r3_D;             // Real bases
    complex_t **c3_C, **c3_D;           // Complex Bases                   
    complex_t *C_eps, *D_eps;           // Polarization vector on different Cartesian bases    
    complex_t *Dp_eps, *Cp_eps;         // Polarization vector on different polarization bases
    complex_t **A1, **A1_i, **A2;       // Change-of-basis matrices

    //
    // Bases
    //

    eB = r3_normalize(eB);
    eK = r3_normalize(beam.k_dic);
    r3_C = orthonormal_basis(eK);
    r3_D = orthonormal_basis(eB);

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
    // Change-of-basis matrix from the Cartesian beam frame to the Cartesian magnetic field frame
    //

    // Complex orthonormal bases
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
            A2[i][j] = c3_inner_product(c3_C[i], c3_D[j]);
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
    eps_probs = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) eps_probs[i] = pow(c_mod(Dp_eps[i]), 2);

    // Release memory
    for(i = 0; i < 3; i++){
        free(A1[i]);
        free(A1_i[i]);
        free(A2[i]);
        free(c3_C[i]);
        free(c3_D[i]);
    }

    free(A1);
    free(A1_i);
    free(A2);
    free(c3_C);
    free(c3_D);

    // Return
    return eps_probs;
}

double scattering_rate(atom_t atom, beam_t beam, double *B, int pol){
    //
    // Variables
    //

    double R;                                           // Scattering rate
    double r, s;                                        // Saturation parameter variables
    double delta, gamma, doppler_shift, zeeman_shift;   // Detuning and Transition rate
    double lambda;                                      // Resonant wave length
    int mj_gnd, mj_exc;                                 // Zeeman shift
    double g_gnd, g_exc;                                // Zeeman shift
    double **C;                                         // Basis of the beam frame

    // Basis of the beam frame
    C = orthonormal_basis(r3_normalize(beam.k_dic));

    //
    // Saturation parameter
    //

    // Distance from the propagation axis
    r = pow(r3_inner_product(C[0], atom.pos), 2);
    r += pow(r3_inner_product(C[1], atom.pos), 2);

    s = beam.s_0;
    s = s * exp(-2 * pow((r / beam.w), 2));

    //
    // Detuning (delta / gamma)
    //

    delta = 0;
    gamma = atom.transition.gamma;      // kHz
    lambda = atom.transition.lambda;    // nm

    // Laser detuning
    delta += beam.delta;

    // Doppler shift
    doppler_shift = -1e4 * r3_inner_product(atom.vel, C[2]) / (lambda * gamma); // Hz
    delta += doppler_shift;

    // Zeeman shift
    mj_gnd = -atom.transition.J_gnd;        // Ground state
    mj_exc = mj_gnd + pol;                  // Excited state
    g_gnd = atom.transition.g_gnd;          // Landè factor of the ground state
    g_exc = atom.transition.g_exc;          // Landè factor of the excited state

    zeeman_shift = 1e4 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc) / gamma;  // Hz
    delta += zeeman_shift;

    // Scattering rate
    R = ((1e3 * gamma)/2) * s / (1 + s + 4*delta*delta); // Hz

    //printf("R = %f\n", R);
    //printf("mj_gnd = %d\n", mj_gnd);
    //printf("mj_exc = %d\n", mj_exc);
    //printf("s = %f\n", s);
    //printf("doppler_shift = %f\n", doppler_shift);
    //printf("zeeman_shift = %f\n", zeeman_shift);
    //printf("|B| (G / cm) = %f\n", r3_mod(B));

    // Release memory    
    free(C);

    // Return
    return R;
}

int print_status(atom_t atom, results_t res){
    printf("Simulation status\n--\n");
    printf("number of iterations = %d\n", res.num_iters);
    printf("total time (s) = %f\n", res.time);
    r3_print(atom.pos, "atom position (cm)");
    r3_print(atom.vel, "atom velocity (cm / s)");
    printf("distance from origin (cm) = %f\n", r3_mod(atom.pos));
    printf("atom speed (m / s) = %f\n", 1e-2 * r3_mod(atom.vel));
    printf("\n");

    return 1;
}

int write_results(char *dir_code, results_t res){
    //
    // Variables
    //

    int i, size;                // Auxiliary variables
    char *path;                 // File stream
    FILE *fp;                   // File stream

    // Open file
    path = str_concatenate(str_concatenate(concatenate_ROOT_PATH("results/"), dir_code), "/position.csv");

    fp = fopen(path, "w");
    
    if(fp == NULL){
        printf("Error to open file \"%s\"", path);
        exit(0);
    }

    // Write file
    fprintf(fp, "id,x,y,z\n");

    size = res.pos_hist[0].num_bins;
    for(i = 0; i < size; i++){
        fprintf(fp, "%d,", i+1);
        fprintf(fp, "%d,", res.pos_hist[0].freqs[i]);
        fprintf(fp, "%d,", res.pos_hist[1].freqs[i]);
        fprintf(fp, "%d\n", res.pos_hist[2].freqs[i]);
    }

    // Close file
    fclose(fp);

    return 1;
}

//
// Utility functions
//

int *get_int_array(char *str, int *size){
    //
    // Variables
    //
    int i, j, max_size = 124;
    char *token;
    int aux_arr[124];
    int *arr;

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

double *get_double_array(char *str, int *size){
    //
    // Variables
    //
    int i, j, max_size = 124;
    char *token;
    double aux_arr[max_size];
    double *arr;

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

    arr = (double *) malloc(i * sizeof(double));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    return arr;
}

// Concatenate strings
char *str_concatenate(char *str1, char *str2){
    // Variables    
    int i, j, size;
    char *str;

    //
    // Concatenation
    //

    size = (int) (strlen(str1) + strlen(str2) - 1);
    str = (char*) malloc(size * sizeof(char));
    
    for(i = 0; i < ((int) strlen(str1)); i++)
        str[i] = str1[i];

    for(j = 0; j < ((int) strlen(str2)); j++)
        str[i + j] = str2[j];

    str[size+1] = '\0';

    return str;
}

char *concatenate_ROOT_PATH(char *filename){
    int i;
    char aux_path[124];
    char *path;

    aux_path[0] = '\0';

    strcat(aux_path, ROOT_PATH);
    strcat(aux_path, filename);

    i = 1;
    while((aux_path[i-1] != '\0') && (i < 124)) i++;

    path = (char*) malloc(i*sizeof(char));
    strcpy(path, aux_path);

    return path;
}

double random_norm(double mean, double std_dev){
    //
    // Variables
    //

    int i;
    double v[2], r, theta;   // Variables for Box-Muller method
    double std_norm;         // Normal(0, 1)
    double norm;             // Adjusted normal

    //
    // Box-Muller transform
    //

    // Generate uniform random numbers
    for(i = 0; i < 2; i++) v[i] = ((double) rand()) / ((double) RAND_MAX);

    // Compute r
    r = sqrt(-2 * log(v[0]));

    // Generate theta
    theta = 0.0;
    while(theta == 0.0) theta = 2.0 * PI * v[1];

    // Generate std_norm value
    std_norm = r * cos(theta);

    // Adjust std_norm
    norm = (std_norm * std_dev) + mean;

    return norm;
}

double random_exp(double mean){
    return - mean * log(1 - ((double) rand()) / ((double) RAND_MAX));
}

int update_hist(histogram_t *hist, double val){
    //
    // Variables
    //

    int bin;
    double lower_lim, upper_lim;

    // Add frequency
    for(bin = 0; bin < (*hist).num_bins; bin++){
        lower_lim = (*hist).coord0 + bin * (*hist).bin_size;
        upper_lim = lower_lim + (*hist).bin_size;

        if((val >= lower_lim) && (val < upper_lim)){
            (*hist).freqs[bin] += 1;
            break;
        }
    }

    return 1;
}

double **orthonormal_basis(double *v){
    //
    // Variables
    //

    int i;
    double **B;       // Desired basis
    double *v1, *v2, *v3;    // Auxiliary vectors

    // Normalize vector v
    v3 = r3_normalize(v);

    // Generate a random vector  
    v1 = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) v1[i] = ((double) rand()) / ((double) RAND_MAX);

    // Define a orthonormal vector
    v2 = r3_scalar_product(r3_inner_product(v1, v3), v3);
    v1 = r3_normalize(r3_diff(v1, v2));

    // Define the last vector of the basis
    v2 = r3_cross_product(v3, v1);

    // Define basis
    B = (double **) malloc(3 * sizeof(double *));
    
    B[0] = v1;
    B[1] = v2;
    B[2] = v3;

    // Release memory
    free(v);
    free(v1);
    free(v2);
    free(v3);

    // Return
    return B;
}

int random_pick(int *arr, double *probs, int size){
    //
    // Variables
    //

    int i, picked = 0;
    double module;
    double rd_n;          // Random number
    double *cum_probs;    // Cumulative probabilities

    // Normalize probabilities
    module = 0;
    for(i = 0; i < size; i++) module += probs[i];
    for(i = 0; i < size; i++) probs[i] = probs[i] / module;

    // Generate a random number  
    rd_n = ((double) rand()) / ((double) RAND_MAX);

    // Cumulative probabilities
    cum_probs = (double*) calloc(size, sizeof(double));
    cum_probs[0] = probs[0];
    for(i = 1; i < size; i++){
        cum_probs[i] += cum_probs[i-1] + probs[i];
    }

    // Pick
    for(i = 0; i < size; i++){
        if(rd_n < cum_probs[i]){
            picked = arr[i];
            break;
        }
    }

    // Release memory
    free(cum_probs);

    // Return
    return picked;
}