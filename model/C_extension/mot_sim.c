/**
 * Monte Carlos simulation of a single atom in a Magneto-Optical Trap (MOT)
 * Bruno N. Santos <nicolau.bruno@gmail.com>
 * Version 2.0
 */

//  Header

#include "mot_sim.h"

// --

results_t simulate_atom(char *params_path, int only_marginals, long seed_time){
    //
    // Variables
    int i;
    double r, dt = 0, passed_time = 0;     // Dynamics
    int get_values = 0;             // Simulation dynamics
    results_t res;                  // Results

    // Seed random variable
    srand(time(0) + seed_time);

    // Parameters of the simulation
    conditions_t conds = get_conditions(params_path);
    environment_t env = get_environment(params_path);
    beams_setup_t beams_setup = get_beams(params_path);
    atom_t atom = get_atom(conds, env, params_path);

    //print_params(atom, conds, beams_setup, env);
    //exit(0);

    // Set initial values
    res.time = 0;
    res.transitions = (int*) calloc(3, sizeof(int));
    set_pos_hist(only_marginals, &res, conds);

    // Distance from origin
    r = sqrt(r3_inner_product(atom.pos, atom.pos));

    //
    // Iterations
    while((passed_time < conds.max_time) && (r < conds.r_max)){
        //print_status(atom, res);

        // Move atom
        dt = move(&atom, beams_setup, conds, env);

        // Iterations numbers
        passed_time += dt; // 1 / gamma

        // Distance from origin
        r = r3_mod(atom.pos);

        //
        // Waiting the equilibrium
        if(passed_time > conds.wait_time){
            get_values = 1;
            res.time += dt; // 1 / gamma
        }

        //
        // Update results
        if(get_values == 1 && r < conds.r_max){
            if(only_marginals)
                for(i = 0; i < 3; i++) update_hist(&res.pos_hist[i], atom.pos[i]);
            else update_hist_3d(&res.pos_3Dhist, atom.pos);
        }
    }

    //print_results(atom, res);

    return res;
}

atom_t get_atom(conditions_t conds, environment_t env, char *params_path){
    //
    // Variables
    // --
    int i;
    char row[STRING_BUFFER_SIZE];
    char *path = str_concatenate(params_path, "atom.csv");
    char *token, *rest;
    double std_dev, *rd_v;
    atom_t atom;
    FILE *fp;
    //--

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
        } 

        // Python module
        else if(Py_MODULE) atom.Z = (int) atof(str_replace(token, ".", ","));

        // C program
        else atom.Z = (int) atof(token);
    }

    // Mass
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"mass\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) atom.mass = atof(str_replace(token, ".", ","));

        // C program
        else atom.mass = atof(token);
    }

    //
    // Initial position of the atom
    //--
    atom.pos = (double *) calloc(3, sizeof(double));

    //
    // Random vector
    //--
    rd_v = (double*) calloc(3, sizeof(double));
    
    // Generate a random value for each direction
    for(i = 0; i < 3; i++) {
        rd_v[i] = ((double) rand()) / ((double) RAND_MAX);
        if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[i] = -rd_v[i];
    }

    // Normalization
    rd_v = r3_normalize(rd_v); // Normalization
    //--

    if(env.w < conds.r_max)
        atom.pos = r3_scalar_product(random_norm(0, env.w/2), rd_v); // cm
    else 
        atom.pos = r3_scalar_product(random_norm(0, conds.r_max/2), rd_v); // cm
    //--

    //
    // Initial velocity of the atom
    //--
    atom.vel = (double *) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++){
        // Position
        //atom.pos[i] = 0;

        // Velocity
        std_dev = sqrt(k_B * conds.T_0 / (atom.mass * u)) * 10; // cm / s
        atom.vel[i] = random_norm(0, std_dev); // cm / s
    }
    //--

    // Optical transition
    atom.transition = get_transition(params_path);

    // Initial state
    atom.J = atom.transition.J_gnd;
    atom.mJ = -atom.transition.J_gnd;

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
        } 

        // Python module
        else if(Py_MODULE) transition.gamma = atof(str_replace(token, ".", ","));

        // C program
        else transition.gamma = atof(token);
    }

    // Resonant wave length
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"lambda\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) transition.lambda = atof(str_replace(token, ".", ","));

        // C program
        else transition.lambda = atof(token);
    }

    // Total angular momentum of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"J_gnd\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) transition.J_gnd = (int) atof(str_replace(token, ".", ","));

        // C Program
        else transition.J_gnd = (int) atof(token);        
    }

    // Total angular momentum of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"J_exc\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) transition.J_exc = (int) atof(str_replace(token, ".", ","));

        // C Program
        else transition.J_exc = (int) atof(token);  
    }

    // Landè factor of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_gnd\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) transition.g_gnd = atof(str_replace(token, ".", ","));

        // C Program
        else transition.g_gnd = atof(token);  
    }

    // Landè factor of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"g_exc\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) transition.g_exc = atof(str_replace(token, ".", ","));

        // C Program
        else transition.g_exc = atof(token);  
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
        }

        // Python module
        else if(Py_MODULE) conditions.T_0 = atof(str_replace(token, ".", ","));

        // C program
        else conditions.T_0 = atof(token);
    }

    // Maximum time of simulation
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"max_time\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) conditions.max_time = atof(str_replace(token, ".", ","));

        // C Program
        else conditions.max_time = atof(token);
    }

    // Maximum distance
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"r_max\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) conditions.r_max = atof(str_replace(token, ".", ","));

        // C Program
        else conditions.r_max = atof(token);
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
        }

        // Python module
        else if(Py_MODULE) conditions.num_bins = (int) atof(str_replace(token, ".", ","));

        // C Program
        else conditions.num_bins = (int) atof(token);
    }

    // Time to reach the equilibrium
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"wait_time\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) conditions.wait_time = atof(str_replace(token, ".", ","));

        // C Program
        else conditions.wait_time = atof(token);
    }

    // Time interval
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"max_dt\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) conditions.dt = atof(str_replace(token, ".", ","));

        // C Program
        else conditions.dt = atof(token);
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
        } 

        // Python module
        else if(Py_MODULE) env.B_0 = atof(str_replace(token, ".", ","));

        // C program
        else env.B_0 = atof(token);
    }

    // Magnetic field axial direction
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"B_axial\" in the file \"%s\"\n", path);
            exit(0);
        } else env.B_basis = orthonormal_basis(r3_normalize(get_double_array(token, &n)));
    }

    // Local magnetic field gradient
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"local_B\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) env.local_B = atof(str_replace(token, ".", ","));

        // C program
        else env.local_B = atof(token);
    }

    // Laser detuning
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"delta\" in the file \"%s\"\n", path);
            exit(0);
        } 

        // Python module
        else if(Py_MODULE){
            if(token[0] == '-') env.delta= -atof(str_replace(token+1, ".", ",")); 
            else env.delta = atof(str_replace(token, ".", ","));
        } 

        // C program
        else {
            if(token[0] == '-') env.delta= -atof(token+1); 
            else env.delta = atof(token);
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
        } 

        // Python module
        else if(Py_MODULE) env.s_0 = atof(str_replace(token, ".", ","));

        // C program
        else env.s_0 = atof(token);
    }

    // Waist Radius
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"w\" in the file \"%s\"\n", path);
            exit(0);
        }

        // Python module
        else if(Py_MODULE) env.w = atof(str_replace(token, ".", ","));

        // C program
        else env.w = atof(token);
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
        //--
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        if(!token){
            printf("Invalid parameter \"k_dic\" in the file \"%s\"\n", path);
            exit(0);
        } else {
            beams[num_beams].k_dic = r3_normalize(get_double_array(token, &n));
        }
        //--

        //
        // Polarization vector
        //--
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Invalid parameter \"eps\" in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].eps = r3_normalize(get_double_array(token, &n));
        //--

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

int set_pos_hist(int only_marginals, results_t *res, conditions_t conds){
    // Variables
    int i, j, k;

    //
    // Only marginals
    //--
    if(only_marginals){
        res->pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

        for(i = 0; i < 3; i++){
            res->pos_hist[i].num_bins = conds.num_bins;
            res->pos_hist[i].bin_size = 2 * conds.r_max / res->pos_hist[i].num_bins;
            res->pos_hist[i].coord0 = - conds.r_max;
            res->pos_hist[i].freqs = (int*) calloc(res->pos_hist[i].num_bins, sizeof(int));
        }
    //--

    //
    // Complete histograms
    //--
    } else{
        res->pos_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res->pos_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res->pos_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        for(i = 0; i < 3; i++){
            res->pos_3Dhist.num_bins[i] = conds.num_bins;
            res->pos_3Dhist.bins_size[i] = 2 * conds.r_max / res->pos_3Dhist.num_bins[i];
            res->pos_3Dhist.coord0[i] = - conds.r_max;
        }

        res->pos_3Dhist.freqs = (int***) calloc(res->pos_3Dhist.num_bins[0], sizeof(int**));

        for(i = 0; i < res->pos_3Dhist.num_bins[0]; i++){
            res->pos_3Dhist.freqs[i] = (int**) calloc(res->pos_3Dhist.num_bins[1], sizeof(int*));
            for(j = 0; j < res->pos_3Dhist.num_bins[1]; j++){
                res->pos_3Dhist.freqs[i][j] = (int*) calloc(res->pos_3Dhist.num_bins[2], sizeof(int));
                for(k = 0; k < res->pos_3Dhist.num_bins[2]; k++){
                   res->pos_3Dhist.freqs[i][j][k] = 0; 
                }
            }
        }
    }
    //--

    return 0;
}

double move(atom_t *atom, beams_setup_t beams_setup, conditions_t conds, environment_t env){
    //
    // Variables
    int i, j, chosen_beam;
    double *B, *eB;
    double vel_mod, *a_B, *rd_v;
    double *pol_amp, *probs, p = 0;
    int pol_opt[] = {+1, -1, 0};  
    double dt;
    polarized_beam_t *beams;

    // Time interval
    dt = conds.dt / (atom->transition.gamma*1e3);

    //
    // Direction of the magnetic field magnetic field
    //--
    B = magnetic_field(env, atom->pos);
    if(r3_mod(B) == 0){
        eB = env.B_basis[2];
    } else eB = r3_normalize(B);
    //--

    // Beams
    beams = (polarized_beam_t*) calloc(3*beams_setup.num, sizeof(polarized_beam_t));
    probs = (double*) calloc(3*beams_setup.num, sizeof(double));

    //
    // Get probabilities of each transition happen
    // --
    for(i = 0; i < beams_setup.num; i++){
        // Polarizations amplitude
        pol_amp = polarizations_amplitudes(beams_setup.beams[i], env, eB);

        // Check each polarization
        for(j = 0; j < 3; j++){        
            // Define polarized beam
            beams[3*i+j].s_0 = pol_amp[j] * env.s_0;
            beams[3*i+j].k_dic = beams_setup.beams[i].k_dic;
            beams[3*i+j].eps = pol_opt[j];

            // Probability of the atom to absorb
            probs[3*i+j] = dt * scattering_rate(*atom, beams[3*i+j], conds, env, B);
            p += probs[3*i+j];
        }
    }

    // Probability of the transition does not happen
    probs[3*beams_setup.num] = 1 - p;

    // Chose a beam to be absorbed or not
    chosen_beam = random_pick(probs, 3*beams_setup.num+1);
    //--

    //
    // Movement
    //--

    // Magnetic acceleration
    a_B = magnetic_acceleration(*atom, env);

    //
    // Update position
    // --
    for(i = 0; i < 3; i++) {
        // Previous velocity
        atom->pos[i] += atom->vel[i] * dt;

        // Magnetic acceleration
        atom->pos[i] += (a_B[i] * dt*dt) / 2;
    }

    // Gravity
    if(env.g_bool) 
        atom->pos[2] += -(g * dt*dt) / 2;
    // --

    //
    // Update velocity
    //--
    // Gravitational acceleration
    if(env.g_bool)
        atom->vel[2] += - g * dt;

    // Magnetic acceleration
    for(i = 0; i < 3; i++)
        atom->vel[i] += a_B[i] * dt;

    //
    // Photonic recoil
    //--
    if(chosen_beam < 3*beams_setup.num){
        //
        // Random vector
        //--
        rd_v = (double*) calloc(3, sizeof(double));  // Random vector
        for(j = 0; j < 3; j++) {
            rd_v[j] = ((double) rand()) / ((double) RAND_MAX);
            if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[j] = -rd_v[j];
        }

        // Normalization
        rd_v = r3_normalize(rd_v); // Normalization
        //--

        // Add velocity
        vel_mod = 1e4 * h / (atom->transition.lambda * atom->mass * u); // cm / s
        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, beams[chosen_beam].k_dic)); // cm / s
        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, rd_v)); // cm / s
    }
    //--

    //--

    // Release memory
    free(beams);
    
    // Convert time unit
    dt = dt * (atom->transition.gamma*1e3);
    //printf("dt [ms] = %f\n", dt*1e3);

    return dt;
}

double *magnetic_acceleration(atom_t atom, environment_t env){
    // Variables
    int i;
    double *a_B, *del_B, norm;  // Magnetic field
    double g_lande;             // Transition
    double *r_prime;

    // Magnetic field gradient
    if(r3_mod(atom.pos) > 0){
        del_B = (double*) calloc(3, sizeof(double)); // G / cm
        r_prime = (double*) calloc(3, sizeof(double)); // cm

        for(i = 0; i < 3; i++)
            r_prime[i] = r3_inner_product(atom.pos, env.B_basis[i]);

        norm = sqrt(pow(r_prime[0], 2) * pow(r_prime[1], 2) + 4 * pow(r_prime[2], 2));

        // Anti-helmholtz coils
        for(i = 0; i < 3; i++){
            del_B[i] = 0;
            del_B[i] += env.B_basis[0][i] * r_prime[0] / (2 * norm);
            del_B[i] += env.B_basis[1][i] * r_prime[1] / (2 * norm);
            del_B[i] += - env.B_basis[2][i] * 2  * r_prime[2] / norm;
            del_B[i] = env.B_0 * del_B[i];
        }

        // Local magnetic field
        del_B[2] += env.local_B;

        // Magnetic acceleration
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2

        // Atom state
        if(atom.J == atom.transition.J_gnd)
            g_lande = atom.transition.g_gnd;
        else g_lande = atom.transition.g_exc;

        for(i = 0; i < 3; i++)
            a_B[i] = - mu_B * g_lande * atom.mJ * del_B[i] / (atom.mass * u) * 1e3; // cm / s^2

        // Release memory
        free(del_B);
        free(r_prime);
    } else {
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2
        for(i = 0; i < 3; i++) a_B[i] = 0;
    }

    return a_B;
}

double *magnetic_field(environment_t env, double *r){
    //
    // Variables
    //
    int i;
    double *B;          // Magnetic field vector
    double *r_prime;

    B = (double*) calloc(3, sizeof(double));
    r_prime = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++)
        r_prime[i] = r3_inner_product(r, env.B_basis[i]);

    // Anti-Helmholtz coils
    for(i = 0; i < 3; i++){
        B[i] = 0;
        B[i] += r_prime[0] * env.B_basis[0][i] / 2;
        B[i] += r_prime[1] * env.B_basis[1][i] / 2;
        B[i] += - r_prime[2] * env.B_basis[2][i];
        B[i] = env.B_0 * B[i];
    }

    // Local Magnetic field
    B[2] += env.local_B * r[2];

    // Release memory
    free(r_prime);

    return B;
}

double *polarizations_amplitudes(beam_t beam, environment_t env, double *eB){
    //
    // Variables
    int i, j;
    double *eps_amp;
    double **R1, **R2;
    complex_t **C1, **C2, **A_s_c, **A;
    complex_t *C1_eps, *C2_eps, *C2_s_eps;

    // Bases of the beam frame
    R1 = orthonormal_basis(beam.k_dic);
    C1 = r3_oper_to_c3_oper(R1);

    // Bases of the B frame
    R2 = orthonormal_basis(eB);
    C2 = r3_oper_to_c3_oper(R2);

    //
    // Change-of-Basis matrix of the Spherical basis to the Cartesian basis
    //--
    A_s_c = c3_operator_zeros();

    A_s_c[0][0].re = - 1 / sqrt(2);
    A_s_c[0][1].re = 1 / sqrt(2);

    A_s_c[1][0].im = 1 / sqrt(2);
    A_s_c[1][1].im = 1 / sqrt(2);

    A_s_c[2][2].re = 1;
    //--

    //
    // Change-of-Basis matrix of the Cartesian beam frame to the Cartesian B frame
    //--
    A = c3_operator_zeros();

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            A[i][j] = c3_inner_product(C2[i], C1[j]);
        }
    }
    //--

    //
    // Get polarization amplitudes
    //--
    // Polarization on the Cartesian beam frame
    C1_eps = c3_apply_operator(A_s_c, r3_to_c3(beam.eps));

    // Polarization on the Cartesian B frame
    C2_eps = c3_apply_operator(A, C1_eps);

    // Polarization on the Spherical B frame
    C2_s_eps = c3_apply_operator(c3_operator_dagger(A_s_c), C2_eps);

    // Polarization amplitudes
    eps_amp = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) eps_amp[i] = pow(c_mod(C2_s_eps[i]), 2);
    //--

    //
    // Release memory    
    //--
    for(i = 0; i < 3; i++){
        free(A[i]);
        free(A_s_c[i]);
        free(C1[i]);
        free(C2[i]);
        free(R1[i]);
        free(R2[i]);
    }

    free(A);
    free(A_s_c);
    free(C1);
    free(C1_eps);
    free(C2);
    free(C2_eps);
    free(C2_s_eps);
    free(R1);
    free(R2);
    //--

    return eps_amp;
}

double scattering_rate(atom_t atom, polarized_beam_t beam, conditions_t conds, environment_t env, double *B){
    //
    // Variables
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

    s = beam.s_0;
    s = s * exp(-2 * pow((r / env.w), 2));

    //
    // Detuning
    delta = 0;
    gamma = atom.transition.gamma; // kHz
    lambda = atom.transition.lambda; // nm

    // Laser detuning
    delta += env.delta*gamma; // kHz
    //printf("laser_detuning = %f\n", delta);

    // Doppler shift
    doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz
    delta += doppler_shift;
    //printf("doppler_shift = %f\n", doppler_shift);

    //
    // Zeeman shift
    
    // Ground state
    mj_gnd = - atom.transition.J_gnd;
    
    // Excited state
    mj_exc = mj_gnd + beam.eps;

    // Landè factors
    g_gnd = atom.transition.g_gnd;
    g_exc = atom.transition.g_exc;

    // Compute shift
    zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz
    delta += zeeman_shift;
    //printf("|B| = %f\n", r3_mod(B));
    //printf("zeeman_shift = %f\n", zeeman_shift);

    // Scattering rate
    R = ((gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

    // Release //memory  
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
    printf("total time (ms) = %f\n", 1e3*res.time / (2*PI*1e3 * atom.transition.gamma));
    r3_print(atom.pos, "atom position (cm)");
    r3_print(atom.vel, "atom velocity (cm / s)");
    printf("distance from origin (cm) = %f\n", r3_mod(atom.pos));
    printf("atom speed (cm / s) = %f\n", r3_mod(atom.vel));
    printf("\n");

    return 1;
}

int print_results(atom_t atom, results_t res){
    // Variables
    int i, j;

    printf("Simulation status\n--\n");
    printf("total time [ms] = %f\n\n", res.time / (2*PI*atom.transition.gamma));

    for(i = 0; i < 3;i++){
        printf("dist[%d] = [\n", i+1);
        for(j = 0; j < res.pos_hist[i].num_bins; j++)
            printf("%d ", res.pos_hist[i].freqs[j]);
        printf("]\n\n");
    }
    printf("\n");

    return 0;
}

int print_params(atom_t atom, conditions_t conds, beams_setup_t beams_setup, environment_t env){
    // Variable
    int i;

    // Atom
    printf("Atom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass [u] = %f\n", atom.mass);
    r3_print(atom.pos, "pos [cm/s]");
    r3_print(atom.vel, "vel [cm/s]");
    printf("J = %d\n", atom.J);
    printf("mJ = %d\n", atom.mJ);
    printf("\n");

    // Transition
    printf("Transition\n--\n");
    printf("gamma [kHz] = %f\n", 2*PI*atom.transition.gamma);
    printf("lambda [nm] = %f\n", atom.transition.lambda);
    printf("J_gnd = %d\n", atom.transition.J_gnd);
    printf("g_gnd = %f\n", atom.transition.g_gnd);
    printf("J_exc = %d\n", atom.transition.J_exc);
    printf("g_exc = %f\n", atom.transition.g_exc);
    printf("\n");

    // Conditions
    printf("Conditions\n--\n");
    printf("T_0 = %f\n", conds.T_0);
    printf("max_time [ms] = %f\n", conds.max_time / (2*PI*atom.transition.gamma));
    printf("r_max = %f\n", conds.r_max);
    printf("num_bins = %d\n", conds.num_bins);
    printf("wait_time [ms] = %f\n", conds.wait_time / (2*PI*atom.transition.gamma));
    printf("time interval [ms] = %f\n", conds.dt / (2*PI*atom.transition.gamma));
    printf("\n");

    // Environment
    printf("Environment\n--\n");
    printf("B_0 = %f\n", env.B_0);
    r3_operator_print(env.B_basis, "B_basis");
    printf("local_B = %f\n", env.local_B);
    printf("delta = %f\n", env.delta);
    printf("s_0 = %f\n", env.s_0);
    printf("w = %f\n", env.w);
    printf("g_bool = %d\n", env.g_bool);
    printf("\n");

    // Beams
    printf("Beams\n--\n");
    for(i = 0; i < beams_setup.num; i++){
        printf("Beam %d\n", i+1);
        r3_print(beams_setup.beams[i].k_dic, "k");
        r3_print(beams_setup.beams[i].eps, "eps");
        printf("--\n");
    }

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
        if(token[0] == '-') aux_arr[i] = -atoi(token+1); 
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

double random_norm(double mean, double std_dev){
    //
    // Variables
    int i;
    double *v, r, theta;     // Variables for Box-Muller method
    double std_norm;         // Normal(0, 1)
    double norm;             // Adjusted normal

    v = (double*) calloc(2, sizeof(double));

    //
    // Box-Muller transform
    //--
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
    //--

    // Release memory
    free(v);

    return norm;
}

double random_exp(double mean){
    return (- mean * log(1 - ((double) rand()) / ((double) RAND_MAX)));
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

int random_pick(double *probs, int size){
    //
    // Variables
    int i, picked = 0;
    double module;
    double rd_n;          // Random number
    double *cum_probs;    // Cumulative probabilities
    int *idx;

    // Indexes
    idx = (int*) calloc(size, sizeof(int));
    for(i = 0; i < size; i++) idx[i] = i;

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
            picked = idx[i];
            break;
        }
    }

    // Release memory
    free(cum_probs);
    free(idx);

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