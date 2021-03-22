/**
 * Monte Carlos simulation of a single atom in a Magneto-Optical Trap (MOT)
 * Bruno N. Santos <nicolau.bruno@gmail.com>
 * Version 2.0
 */

//  Header

#include "mot_sim.h"

// --

results_t simulate_atom(char *params_path, int only_marginals, long time){
    //
    // Variables
    int i, j, k;                // Iterations
    double r, dt, dt_c;         // Dynamics
    results_t res;              // Results
    double *a_B;                // Magnetic acceleration
    double **B_basis;           // Basis
    int get_values;             // Simulation dynamics
    scattering_t scatt;

    // Seeds random variable
    srand(time);

    // Parameters of the simulation
    conditions_t conds = get_conditions(params_path);
    environment_t env = get_environment(params_path);
    beams_setup_t beams_setup = get_beams(params_path);
    atom_t atom = get_atom(conds, params_path);

    //print_params(atom, conds, beams_setup, env);
    //exit(0);

    // Basis with z-axis parallel to the axial direction of the magnetic field
    B_basis = orthonormal_basis(env.B_axial);

    //
    // Position histogram

    // Only marginals
    if(only_marginals){
        res.pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

        for(i = 0; i < 3; i++){
            res.pos_hist[i].num_bins = conds.num_bins;
            res.pos_hist[i].bin_size = 2 * conds.r_max / res.pos_hist[i].num_bins;
            res.pos_hist[i].coord0 = - conds.r_max;
            res.pos_hist[i].freqs = (int*) calloc(res.pos_hist[i].num_bins, sizeof(int));
            //update_hist(&res.pos_hist[i], atom.pos[i]);
        }

    // Complete histograms
    } else{
        res.pos_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res.pos_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res.pos_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        for(i = 0; i < 3; i++){
            res.pos_3Dhist.num_bins[i] = conds.num_bins;
            res.pos_3Dhist.bins_size[i] = 2 * conds.r_max / res.pos_3Dhist.num_bins[i];
            res.pos_3Dhist.coord0[i] = - conds.r_max;
        }

        res.pos_3Dhist.freqs = (int***) calloc(res.pos_3Dhist.num_bins[0], sizeof(int**));

        for(i = 0; i < res.pos_3Dhist.num_bins[0]; i++){
            res.pos_3Dhist.freqs[i] = (int**) calloc(res.pos_3Dhist.num_bins[1], sizeof(int*));
            for(j = 0; j < res.pos_3Dhist.num_bins[1]; j++){
                res.pos_3Dhist.freqs[i][j] = (int*) calloc(res.pos_3Dhist.num_bins[2], sizeof(int));
                for(k = 0; k < res.pos_3Dhist.num_bins[2]; k++){
                   res.pos_3Dhist.freqs[i][j][k] = 0; 
                }
            }
        }

        // Insert initial position
        //update_hist_3d(&res.pos_3Dhist, atom.pos);
    }

    // Distance from origin
    r = sqrt(r3_inner_product(atom.pos, atom.pos)); // Distance from centre

    //
    // Set initial variables
    res.num_iters = 0; 
    dt = 0;
    dt_c = 0;
    res.time = 0;
    get_values = 0;
    atom.state = 0;

    //print_params(atom, conds, beams_setup, env);

    //
    // Iterations
    while((res.num_iters < conds.i_max) && (r < conds.r_max)){
        // Stimulated absorption
        if(atom.state == 0){
            // Photon absorption
            scatt = photon_absorption(atom, beams_setup, conds, env, B_basis);
        
            //
            // Update time
            if(scatt.R > 0 && scatt.dt < (MAX_dt + dt_c)){
                atom.state = 1;
                dt = fabs(scatt.dt - dt_c);
                dt_c = 0;
            } else {
                dt = 1 / (2*PI*atom.transition.gamma*1e3); // s
                dt_c += dt;
            }

            // Magnetic acceleration
            a_B = magnetic_acceleration(atom, env, B_basis);

            //
            // Update position

            for(i = 0; i < 3; i++) {
                // Previous velocity
                atom.pos[i] += atom.vel[i] * dt;

                // Magnetic acceleration
                atom.pos[i] += (a_B[i] * dt*dt) / 2;
            }

            // Gravity
            if(env.g_bool) 
                atom.pos[2] += -(g * dt*dt) / 2;

            //
            // Update velocity

            // Gravitational acceleration
            if(env.g_bool)
                atom.vel[2] += - g * dt;

            for(i = 0; i < 3; i++){
                // Photonic recoil
                if(atom.state == 1) atom.vel[i] += scatt.vel[i];

                // Magnetic acceleration
                atom.vel[i] += a_B[i] * dt;
            }

            //
            // Update results
            if(get_values == 1){
                if(only_marginals)
                    for(i = 0; i < 3; i++) update_hist(&res.pos_hist[i], atom.pos[i]);
                else update_hist_3d(&res.pos_3Dhist, atom.pos);
            }

            // Distance from origin
            r = r3_mod(atom.pos);

            // Iterations numbers
            res.num_iters += 1;

            //
            // Waiting the equilibrium
            if(res.num_iters > conds.ini_iters){
                res.time += dt;
                get_values = 1;
            }
        
        // Spontaneous emission
        } else if(atom.state == 1){
            scatt = photon_emission(atom);
            atom.state = 0;
            
            //
            // Update time
            dt = scatt.dt;            

            // Magnetic acceleration
            a_B = magnetic_acceleration(atom, env, B_basis);

            //
            // Update position

            for(i = 0; i < 3; i++) {
                // Previous velocity
                atom.pos[i] += atom.vel[i] * dt;

                // Magnetic acceleration
                atom.pos[i] += (a_B[i] * dt*dt) / 2;
            }

            // Gravity
            if(env.g_bool) 
                atom.pos[2] += -(g * dt*dt) / 2;

            //
            // Update velocity

            // Gravitational acceleration
            if(env.g_bool)
                atom.vel[2] += - g * dt;

            for(i = 0; i < 3; i++){
                // Photonic recoil
                atom.vel[i] += scatt.vel[i];

                // Magnetic acceleration
                atom.vel[i] += a_B[i] * dt;
            }

            //
            // Update results
            if(get_values == 1){
                if(only_marginals)
                    for(i = 0; i < 3; i++) update_hist(&res.pos_hist[i], atom.pos[i]);
                else update_hist_3d(&res.pos_3Dhist, atom.pos);
            }

            // Distance from origin
            r = r3_mod(atom.pos);

            // Iterations numbers
            res.num_iters += 1;

            //
            // Waiting the equilibrium
            if(res.num_iters > conds.ini_iters){
                res.time += dt;
                get_values = 1;
            }
        }

        if(r3_inner_product(atom.pos, atom.vel) > 0 && r > env.w)
            break;
        
        /* Print status
        if(res.num_iters % 100 == 0) {
            printf("dt = %f\n", dt);
            printf("R = %f\n", scatt.R);
            printf("tau = %f\n", scatt.dt);
            r3_print(a_B, "a_B");
            printf("a_B*dt = [ ");
            for(i = 0; i < 3; i++)
                printf("%f ", a_B[i]*dt);
            printf(" ]\n");

            print_status(atom, res);
        }*/
    }

    // Print status
    //print_status(atom, res);

    // Release memory
    for(i = 0; i < 3; i++) free(B_basis[i]);
    free(B_basis);
    free(a_B);

    return res;
}

atom_t get_atom(conditions_t conds, char *params_path){
    //
    // Variables
    //

    atom_t atom;

    int i;

    char row[STRING_BUFFER_SIZE];
    char *path = str_concatenate(params_path, "atom.csv");
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
        } else atom.Z = (int) atof(token);
    }

    // Mass
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"mass\" in the file \"%s\"\n", path);
            exit(0);
        } else atom.mass = atof(token);
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
    atom.transition = get_transition(params_path);

    // Close file
    fclose(fp);

    // Release memory
    free(path);

    return atom;
}

transition_t get_transition(char *params_path){
    //
    // Variables
    //

    transition_t transition;
    char row[STRING_BUFFER_SIZE];
    char *path = str_concatenate(params_path, "transition.csv");
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
        } else transition.gamma = atof(str_replace(token, ".", ",")) * 2 * PI;
    }

    // Resonant wave length
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"lambda\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.lambda = atof(str_replace(token, ".", ","));
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
        } else transition.g_gnd = atof(str_replace(token, ".", ","));
    }

    // Landè factor of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_exc\" in the file \"%s\"\n", path);
            exit(0);
        } else transition.g_exc = atof(str_replace(token, ".", ","));
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
    free(path);

    // Return
    return transition;
}

conditions_t get_conditions(char *params_path){
    //
    // Variables
    //

    conditions_t conditions;
    char row[STRING_BUFFER_SIZE];
    char *path = str_concatenate(params_path, "conditions.csv");
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
        } else conditions.T_0 = atof(str_replace(token, ".", ","));
    }

    // Maximum number of iteration
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"i_max\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.i_max = (int) atof(str_replace(token, ".", ","));
    }

    // Maximum distance
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"r_max\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.r_max = atof(str_replace(token, ".", ","));
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
        } else conditions.num_bins = (int) atof(token);
    }

    // Equilibrium iterations
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"ini_iters\" in the file \"%s\"\n", path);
            exit(0);
        } else conditions.ini_iters = (int) atof(token);
    }

    // Close file
    fclose(fp);

    // Release memory
    free(path);

    return conditions;
}

environment_t get_environment(char *params_path){
    //
    // Variables
    //

    int n;
    environment_t env;
    char row[STRING_BUFFER_SIZE];
    char *path = str_concatenate(params_path, "environment.csv");
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

    // Magnetic field gradient
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"B_0\" in the file \"%s\"\n", path);
            exit(0);
        } else env.B_0 = atof(str_replace(token, ".", ","));
    }

    // Magnetic field axial direction
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"B_axial\" in the file \"%s\"\n", path);
            exit(0);
        } else env.B_axial = r3_normalize(get_double_array(token, &n));
    }

    // Local magnetic field gradient
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"local_B\" in the file \"%s\"\n", path);
            exit(0);
        } else env.local_B = atof(str_replace(token, ".", ","));
    }

    // Laser detuning
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"delta\" in the file \"%s\"\n", path);
            exit(0);
        } else { 
            if(token[0] == '-') env.delta= -atof(str_replace(token+1, ".", ",")); 
            else env.delta = atof(str_replace(token, ".", ","));
        }
    }

    // Peak of the saturation parameter
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"s_0\" in the file \"%s\"\n", path);
            exit(0);
        } else {
            if(token[0] == '-') env.s_0 = -atof(str_replace(token+1, ".", ",")); 
            else env.s_0 = atof(str_replace(token, ".", ","));
        }
    }

    // Waist Radius
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"w\" in the file \"%s\"\n", path);
            exit(0);
        } else env.w = atof(str_replace(token, ".", ","));
    }

    // Gravity
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_bool\" in the file \"%s\"\n", path);
            exit(0);
        } else env.g_bool = atoi(token);
    }

    // Close file
    fclose(fp);

    // Release memory
    free(path);

    return env;
}

beams_setup_t get_beams(char *params_path){
    //
    // Variables
    //

    beams_setup_t beams_setup;
    beam_t *beams = (beam_t*) calloc(MAX_BEAMS, sizeof(beam_t));
    beam_t *c_beams;

    int num_beams = 0, n;
    char row[STRING_BUFFER_SIZE];

    char *path = str_concatenate(params_path, "beams.csv");
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
        // Wave vector direction
        rest = row;
        //token = strtok_r(rest, DELIM, &rest); // Id
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"k_dic\" in the file \"%s\"\n", path);
            exit(0);
        } else {
            beams[num_beams].k_dic = r3_normalize(get_double_array(token, &n));
        }

        //
        // Polarization vector
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"eps\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].eps = r3_normalize(get_double_array(token, &n));

        num_beams++;
    }

    c_beams = (beam_t *) calloc(num_beams, sizeof(beam_t));
    for(n = 0; n < num_beams; n++) c_beams[n] = beams[n];

    beams_setup.num = num_beams;
    beams_setup.beams = c_beams;

    // Close file
    fclose(fp);

    // Release memory
    free(beams);
    free(path);

    // Return
    return beams_setup;
}

scattering_t *photonic_recoil(atom_t atom, beams_setup_t beams_setup, conditions_t conds, environment_t env, double **B_basis){
    //
    // Variables
    //

    int i, j;                                   // Iterations variables
    double *B, *eps_probs, *eB, *eK, lambda;    // Electromagnetic fields
    scattering_t *scatts, aux_scatt;            // Scattering options
    double vel_mod, *rd_v;                      // Auxiliary variables
    int *pol_opt, pol;                          // Polarization

    // Get Magnetic field
    B = get_magnetic_field(env, B_basis, atom.pos);

    //
    // Get scattering length
    //

    scatts = (scattering_t*) calloc(beams_setup.num, sizeof(scattering_t));

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
        eps_probs = polarization_probs(beams_setup.beams[i], env, eB);
        pol = random_pick(pol_opt, eps_probs, 3);

        // Compute scattering
        scatts[i].R = scattering_rate(atom, beams_setup.beams[i], conds, env, B, pol); // Hz
        if(scatts[i].R > 0) scatts[i].dt = random_exp(1 / scatts[i].R); // s
        else scatts[i].dt = 0;

        // Absorption event
        vel_mod = 1e4 * h / (lambda * atom.mass * u); // cm / s
        scatts[i].vel = r3_scalar_product(vel_mod, eK); // cm / s

        // Emission event
        rd_v = (double*) calloc(3, sizeof(double));  // Random vector
        for(j = 0; j < 3; j++) {
            rd_v[j] = ((double) rand()) / ((double) RAND_MAX);
            if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[j] = -rd_v[j];
        }
        rd_v = r3_normalize(rd_v); // Normalization
        rd_v = r3_scalar_product(vel_mod, rd_v); // cm / s
        scatts[i].vel = r3_sum(scatts[i].vel, rd_v); // cm / s
    }

    //
    // Sort scatterings
    for(i = 0; i < beams_setup.num - 1; i++){
        for(j = 0; j < beams_setup.num - 1; j++){
            if((scatts[j].dt > scatts[j+1].dt) || (scatts[j].R == 0)){
                aux_scatt = scatts[j];
                scatts[j] = scatts[j+1];
                scatts[j+1] = aux_scatt;
            }
        }
    }

    // Release memory
    free(eK);
    free(eB);
    free(B);
    free(eps_probs);
    free(pol_opt);
    free(rd_v);

    // Return
    return scatts;
}

scattering_t photon_absorption(atom_t atom, beams_setup_t beams_setup, conditions_t conds, environment_t env, double **B_basis){
    //
    // Variables
    //

    int i, j;                                   // Iterations variables
    double *B, *eps_probs, *eB, *eK, lambda;    // Electromagnetic fields
    scattering_t *scatts, aux_scatt;            // Scattering options
    double vel_mod, *rd_v;                      // Auxiliary variables
    int *pol_opt, pol;                          // Polarization

    // Get Magnetic field
    B = get_magnetic_field(env, B_basis, atom.pos);

    //
    // Get scattering length
    //

    scatts = (scattering_t*) calloc(beams_setup.num, sizeof(scattering_t));

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
        eps_probs = polarization_probs(beams_setup.beams[i], env, eB);
        pol = random_pick(pol_opt, eps_probs, 3);

        // Compute scattering
        scatts[i].R = scattering_rate(atom, beams_setup.beams[i], conds, env, B, pol); // Hz
        if(scatts[i].R > 0) scatts[i].dt = random_exp(1 / scatts[i].R); // s
        else scatts[i].dt = 0;

        // Absorption event
        vel_mod = 1e4 * h / (lambda * atom.mass * u); // cm / s
        scatts[i].vel = r3_scalar_product(vel_mod, eK); // cm / s
    }

    //
    // Sort scatterings
    for(i = 0; i < beams_setup.num - 1; i++){
        for(j = 0; j < beams_setup.num - 1; j++){
            if((scatts[j].dt > scatts[j+1].dt) || (scatts[j].R == 0)){
                aux_scatt = scatts[j];
                scatts[j] = scatts[j+1];
                scatts[j+1] = aux_scatt;
            }
        }
    }

    // Release memory
    free(eK);
    free(eB);
    free(B);
    free(eps_probs);
    free(pol_opt);
    free(rd_v);

    // Return
    return scatts[0];
}

scattering_t photon_emission(atom_t atom){
    //
    // Variables
    //

    int j;                  // Iterations variables
    double lambda;          // Electromagnetic fields
    scattering_t scatt;     // Scattering
    double vel_mod, *rd_v;  // Auxiliary variables

    // Get wave length
    lambda = atom.transition.lambda; // nm

    // Speed gain
    vel_mod = 1e4 * h / (lambda * atom.mass * u); // cm / s

    //
    // Random vector
    rd_v = (double*) calloc(3, sizeof(double));
    for(j = 0; j < 3; j++) {
        rd_v[j] = ((double) rand()) / ((double) RAND_MAX);
        if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[j] = -rd_v[j];
    }
    rd_v = r3_normalize(rd_v); // Normalization

    // Time of the process
    scatt.dt = random_exp(1 / (2*PI*atom.transition.gamma*1e3));
    
    // Emission
    scatt.vel = r3_scalar_product(vel_mod, rd_v); // cm / s

    // Return
    return scatt;
}

double *magnetic_acceleration(atom_t atom, environment_t env, double **B_basis){
    // Variables
    int i;
    double *a_B, *del_B, norm;  // Magnetic field
    double g_gnd, mJ_gnd;       // Transition
    double *r_prime;

    // Magnetic field gradient
    if(r3_mod(atom.pos) > 0){
        del_B = (double*) calloc(3, sizeof(double)); // G / cm
        r_prime = (double*) calloc(3, sizeof(double)); // cm

        for(i = 0; i < 3; i++)
            r_prime[i] = r3_inner_product(atom.pos, B_basis[i]);

        norm = sqrt(pow(r_prime[0], 2) * pow(r_prime[1], 2) + 4 * pow(r_prime[2], 2));

        // Anti-helmholtz coils
        for(i = 0; i < 3; i++){
            del_B[i] = 0;
            del_B[i] += B_basis[0][i] * r_prime[0] / (2 * norm);
            del_B[i] += B_basis[1][i] * r_prime[1] / (2 * norm);
            del_B[i] += - B_basis[2][i] * 2  * r_prime[2] / norm;
            del_B[i] = env.B_0 * del_B[i];
        }

        // Local magnetic field
        del_B[2] += env.local_B;

        // Magnetic acceleration
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2

        g_gnd = atom.transition.g_gnd;
        mJ_gnd = -atom.transition.J_gnd;

        for(i = 0; i < 3; i++)
            a_B[i] = - mu_B * g_gnd * mJ_gnd * del_B[i] / (atom.mass * u) * 1e3; // cm / s^2

        // Release memory
        free(del_B);
        free(r_prime);
    } else {
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2
        for(i = 0; i < 3; i++) a_B[i] = 0;
    }

    return a_B;
}

double *get_magnetic_field(environment_t env, double **B_basis, double *r){
    //
    // Variables
    //
    int i;
    double *B;          // Magnetic field vector
    double *r_prime;

    B = (double*) calloc(3, sizeof(double));
    r_prime = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++)
        r_prime[i] = r3_inner_product(r, B_basis[i]);

    // Anti-Helmholtz coils
    for(i = 0; i < 3; i++){
        B[i] = 0;
        B[i] += r_prime[0] * B_basis[0][i] / 2;
        B[i] += r_prime[1] * B_basis[1][i] / 2;
        B[i] += - r_prime[2] * B_basis[2][i];
        B[i] = env.B_0 * B[i];
    }

    // Local Magnetic field
    B[2] += env.local_B * r[2];

    // Release memory
    free(r_prime);

    return B;
}

double *polarization_probs(beam_t beam, environment_t env, double *eB){
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
        free(r3_C[i]);
        free(r3_D[i]);
        free(c3_C[i]);
        free(c3_D[i]);
    }

    free(A1);
    free(A1_i);
    free(A2);
    free(c3_C);
    free(c3_D);
    free(r3_C);
    free(r3_D);
    free(eK);
    free(C_eps);
    free(Cp_eps);
    free(D_eps);
    free(Dp_eps);

    // Return
    return eps_probs;
}

double scattering_rate(atom_t atom, beam_t beam, conditions_t conds, environment_t env, double *B, int pol){
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
    int i;

    // Basis of the beam frame
    C = orthonormal_basis(r3_normalize(beam.k_dic));

    //
    // Saturation parameter
    //

    // Distance from the propagation axis
    r = pow(r3_inner_product(C[0], atom.pos), 2);
    r += pow(r3_inner_product(C[1], atom.pos), 2);

    s = env.s_0;
    s = s * exp(-2 * pow((r / env.w), 2));

    //
    // Detuning (delta / gamma)
    //

    delta = 0;
    gamma = 2*PI*atom.transition.gamma; // kHz
    lambda = atom.transition.lambda;    // nm

    // Laser detuning
    delta += 2 * env.delta;

    // Doppler shift
    doppler_shift = -1e4 * 4 * PI * r3_inner_product(atom.vel, C[2]) / (lambda * gamma);
    delta += doppler_shift;

    // Zeeman shift
    mj_gnd = -atom.transition.J_gnd;        // Ground state
    mj_exc = mj_gnd + pol;                  // Excited state
    g_gnd = atom.transition.g_gnd;          // Landè factor of the ground state
    g_exc = atom.transition.g_exc;          // Landè factor of the excited state

    zeeman_shift = 1e4 * 2 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc) / gamma;  // Hz
    delta += zeeman_shift;

    // Scattering rate
    R = ((1e3 * gamma)/2) * s / (1 + s + delta*delta); // Hz

    // Release memory  
    for(i = 0; i < 3; i++) free(C[i]); 
    free(C);

    // Return
    return R;
}

//
// Utility functions
//

int print_status(atom_t atom, results_t res){
    printf("Simulation status\n--\n");
    printf("number of iterations = %d\n", res.num_iters);
    printf("total time (ms) = %f\n", 1e3*res.time);
    r3_print(atom.pos, "atom position (cm)");
    r3_print(atom.vel, "atom velocity (cm / s)");
    printf("distance from origin (cm) = %f\n", r3_mod(atom.pos));
    printf("atom speed (m / s) = %f\n", 1e-2 * r3_mod(atom.vel));
    printf("\n");

    return 1;
}

int print_results(results_t res){
    printf("Simulation status\n--\n");
    printf("number of iterations = %d\n", res.num_iters);
    printf("total time (ms) = %f\n", res.time*1e3);
    printf("\n");

    return 0;
}

int print_params(atom_t atom, conditions_t conds, beams_setup_t beams, environment_t env){
    printf("Conditions\n--\n");
    printf("T_0 = %f\n", conds.T_0);
    printf("i_max = %d\n", conds.i_max);
    printf("r_max = %f\n", conds.r_max);
    printf("num_bins = %d\n", conds.num_bins);
    printf("ini_iters = %d\n", conds.ini_iters);
    printf("\n");

    printf("Environment\n--\n");
    printf("B_0 = %f\n", env.B_0);
    r3_print(env.B_axial, "B_axial");
    printf("delta = %f\n", env.delta);
    printf("s_0 = %f\n", env.s_0);
    printf("w = %f\n", env.w);
    printf("g_bool = %d\n", env.g_bool);
    printf("\n");
    return 0;
}

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

    // Release memory
    free(token);

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
        if(token[0] == '-') aux_arr[i] = -atof(str_replace(token+1, ".", ",")); 
        else aux_arr[i] = atof(str_replace(token, ".", ","));

        token = strtok_r(str, " ", &str);
        i++;
    }

    arr = (double *) malloc(i * sizeof(double));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    // Release memory
    free(token);

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
    str = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));
    
    for(i = 0; i < ((int) strlen(str1)); i++)
        str[i] = str1[i];

    for(j = 0; j < ((int) strlen(str2)); j++)
        str[i + j] = str2[j];

    str[size+1] = '\0';

    return str;
}

char *concatenate_ROOT_PATH(char *filename){
    int i;
    char *path, *aux_path;

    aux_path = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));
    aux_path[0] = '\0';

    strcat(aux_path, ROOT_PATH);
    strcat(aux_path, filename);

    i = 1;
    while((aux_path[i-1] != '\0') && (i < 124)) i++;

    path = (char*) malloc(i*sizeof(char));
    strcpy(path, aux_path);

    // Release memory
    free(aux_path);

    return path;
}

double random_norm(double mean, double std_dev){
    //
    // Variables
    //

    int i;
    double *v, r, theta;     // Variables for Box-Muller method
    double std_norm;         // Normal(0, 1)
    double norm;             // Adjusted normal

    v = (double*) calloc(2, sizeof(double));

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

    // Release memory
    free(v);

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

int update_hist_3d(histogram_3d_t *hist, double *vals){
    //
    // Variables
    //

    int i, bin, *coord, dim = 3;
    double lower_lim, upper_lim;

    // Coordinates
    coord = (int*) calloc(dim, sizeof(int));

    // Check each coordinate
    for(i = 0; i < dim; i++){
        // Add frequency
        for(bin = 0; bin < (*hist).num_bins[i]; bin++){
            lower_lim = (*hist).coord0[i] + bin * (*hist).bins_size[i];
            upper_lim = lower_lim + (*hist).bins_size[i];

            if((vals[i] >= lower_lim) && (vals[i] < upper_lim)){
                coord[i] = bin;
                break;
            }
        }  
    }

    (*hist).freqs[coord[0]][coord[1]][coord[2]] += 1;

    // Release memory
    free(coord);

    return 1;
}

double **orthonormal_basis(double *v){
    //
    // Variables
    //

    int i;
    double **B;             // Desired basis
    double *v1, *v2, *v3;   // Auxiliary vectors

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
    B = (double **) calloc(3, sizeof(double *));

    for(i = 0; i < 3; i++)
        B[i] = (double *) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++){
        B[0][i] = v1[i];
        B[1][i] = v2[i];
        B[2][i] = v3[i];
    }

    // Free
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

char *str_replace(char *orig, char *rep, char *with){
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep (the string to remove)
    int len_with; // length of with (the string to replace rep with)
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    // sanity checks and initialization
    if (!orig || !rep)
        return NULL;
    len_rep = strlen(rep);
    if (len_rep == 0)
        return NULL; // empty rep causes infinite loop during count
    if (!with)
        with = "";
    len_with = strlen(with);

    // count the number of replacements needed
    ins = orig;
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);

    return result;
}