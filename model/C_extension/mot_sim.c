/**
 * Monte Carlos simulation of a single atom in a Magneto-Optical Trap (MOT)
 * Bruno N. Santos <nicolau.bruno@gmail.com>
 * Version 2.0
 */

//  Header

#include "mot_sim.h"

// --

results_t simulate_atom(char *params_path, int opt, long seed_time){
    // Variables
    int i;
    double r, v, dt = 0;
    double progress;
    int last_progress = 0;
    results_t res;

    // Seed random variable
    srand(time(0) + seed_time);

    // Parameters of the simulation
    initial_conditions_t ini_conds = get_initial_conditions(params_path);
    performance_t perform = get_performance(params_path);
    magnetic_field_t B_params = get_magnetic_field(params_path);
    beams_setup_t beams_setup = get_beams(params_path);
    atom_t atom = get_atom(ini_conds, perform, beams_setup, opt, params_path);

    // Print parameters
    //printf("opt = %d\n", opt);
    //print_initial_conditions(ini_conds);
    //print_performance(perform);
    //print_magnetic_field(B_params);
    //print_beams(beams_setup);
    //print_atom(atom);
    //exit(0);

    // Set initial values
    res.time = 0;
    res.trapped_atom = 0;
    set_hist(opt, &res, perform);

    // Distance from origin
    r = r3_mod(atom.pos);

    // Speed of the atom
    v = r3_mod(atom.vel);

    //
    // Iterations
    //--
    while((res.time < perform.max_time) && (r <= perform.max_r) && (v <= perform.max_v)){
        //print_status(atom, res);

        // Move atom
        dt = move(&atom, beams_setup, perform, B_params);

        // Iterations numbers
        res.time += dt; // 1 / gamma

        // Distance from origin
        r = r3_mod(atom.pos);

        // Speed of the atom
        v = r3_mod(atom.vel);

        //
        // Waiting the equilibrium
        if(res.time > perform.wait_time && res.trapped_atom == 0){
            res.trapped_atom = 1;
            if(opt == 2) break;
        }

        //
        // Update results
        if(res.trapped_atom == 1 && r < perform.max_r){
            //
            // Update position and velocity
            //--
            // Marginal histogram
            if(opt == 1){
                for(i = 0; i < 3; i++) {
                    update_hist(&res.pos_hist[i], atom.pos[i]);
                    if(v < perform.max_v) update_hist(&res.vel_hist[i], atom.vel[i]);
                }
            }

            // 3D-Histograms
            else {
                update_hist_3d(&res.pos_3Dhist, atom.pos);
                if(v < perform.max_v) update_hist_3d(&res.vel_3Dhist, atom.vel);
            }
            //--
        }

        progress = (100*res.time / perform.max_time);
        if((((int) progress) % 5) == 0 && last_progress < ((int) progress)){
            //print_status(atom, res);
            //printf("res.time [1/Gamma] = %f\n", res.time);
            //printf("res.time [ms] = %f\n", res.time / (2 * PI * atom.transition.gamma));
            //printf("progress = %f\n\n", progress);
            last_progress = (int) progress;
        }
    }
    //--

    //print_status(atom, res);
    //printf("res.time [1/Gamma] = %f\n", res.time);
    //printf("res.time [ms] = %f\n", res.time / (2 * PI * atom.transition.gamma));
    //printf("progress = %f\n\n", progress);
    //printf("trapped_atom = %d\n\n", res.trapped_atom);
    //print_results(res, atom, opt);

    return res;
}

performance_t get_performance(char *params_path){
    //
    // Variables
    //

    int i = 0;
    char *token, *saveptr, *path, **rows;
    performance_t perform;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "performance.csv");
    rows = read_lines(path);

    //
    // Maximum time of simulation
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.max_time = atof(str_replace(token, ".", ","));

    // C Program
    else perform.max_time = atof(token);
    //--

    //
    // Maximum distance
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.max_r = atof(str_replace(token, ".", ","));

    // C Program
    else perform.max_r = atof(token);
    //--

    //
    // Maximum speed
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.max_v = atof(str_replace(token, ".", ","));

    // C Program
    else perform.max_v = atof(token);
    //--

    // Skip parameter (number of simulations)
    i += 1;

    //
    // Number of bins
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.num_bins = (int) atof(str_replace(token, ".", ","));

    // C Program
    else perform.num_bins = (int) atof(token);
    //--

    //
    // Waiting time (time to reach the equilibrium)
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.wait_time = atof(str_replace(token, ".", ","));

    // C Program
    else perform.wait_time = atof(token);
    //--

    //
    // Time interval
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.dt = atof(str_replace(token, ".", ","));

    // C Program
    else perform.dt = atof(token);
    //--

    // Release memory
    free(path);

    return perform;
}

initial_conditions_t get_initial_conditions(char *params_path){
    //
    // Variables
    //

    int i = 0, n;
    char *token, *saveptr, *path, **rows;
    initial_conditions_t ini_conds;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "initial_conditions.csv");
    rows = read_lines(path);

    //
    // Initial temperature 
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.T_0 = atof(str_replace(token, ".", ","));

    // C program
    else ini_conds.T_0 = atof(token);
    //--

    //
    // Module of the initial velocity of the atoms
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.v_0 = atof(str_replace(token, ".", ","));

    // C program
    else ini_conds.v_0 = atof(token);
    //--

    //
    // Direction of the initial velocity of the atoms
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.v_0_dir = get_double_array(token, &n);

    // C program
    else ini_conds.v_0_dir = get_double_array(token, &n);
    //--

    //
    // Gravity
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.g_bool = atoi(str_replace(token, ".", ","));

    // C program
    else ini_conds.g_bool = atoi(token);
    //--

    // Release memory
    free(path);

    return ini_conds;
}

magnetic_field_t get_magnetic_field(char *params_path){
    // Variables
    int i = 0, n;
    char *token, *saveptr, *path, **rows;
    magnetic_field_t B_params;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "magnetic_field.csv");
    rows = read_lines(path);

    //
    // Magnetic field gradient
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value
    
    // Python module
    if(Py_MODULE) B_params.B_0 = atof(str_replace(token, ".", ","));

    // C program
    else B_params.B_0 = atof(token);
    //--

    //
    // Magnetic field axial direction
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value
        
    B_params.B_basis = orthonormal_basis(r3_normalize(get_double_array(token, &n)));
    //--

    //
    // Bias of magnetic field
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) B_params.B_bias = get_double_array(token, &n);

    // C program
    else B_params.B_bias = get_double_array(token, &n);
    //--

    //
    // Linear magnetic field gradient
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) B_params.B_lin_grad = get_double_array(token, &n);

    // C program
    else B_params.B_lin_grad = get_double_array(token, &n);
    //--

    // Release memory
    free(path);

    return B_params;
}

beams_setup_t get_beams(char *params_path){
    // Variables
    int i = 0, num_beams = 0, n;
    char *token, *saveptr, *path, **rows;
    beams_setup_t beams_setup;
    double s_0, delta, w;
    beam_t *beams = (beam_t*) calloc(MAX_BEAMS, sizeof(beam_t));
    beam_t *c_beams;

    //
    // Get main parameters of the beams
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/main.csv");
    rows = read_lines(path);
    i = 0;

    //
    // Laser detuning
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE){
        if(token[0] == '-') delta= -atof(str_replace(token+1, ".", ",")); 
        else delta = atof(str_replace(token, ".", ","));
    } 

    // C program
    else {
        if(token[0] == '-') delta= -atof(token+1); 
        else delta = atof(token);
    }
    //--

    //
    // Peak of the saturation parameter
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) s_0 = atof(str_replace(token, ".", ","));

    // C program
    else s_0 = atof(token);
    //--

    //
    // Waist Radius
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) w = atof(str_replace(token, ".", ","));

    // C program
    else w = atof(token);
    //--
    //--

    //
    // Get all beams in the setup
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/setup.csv");
    rows = read_lines(path);

    //
    // Get beams
    //--
    for(i = 1; !(rows[i] == NULL); i++){
        //
        // Wave vector direction
        //--
        token = strtok_r(rows[i], DELIM, &saveptr);
        beams[num_beams].k_dir = r3_normalize(get_double_array(token, &n));
        //--

        //
        // Polarization vector
        //--
        token = strtok_r(NULL, DELIM, &saveptr); // Value
        beams[num_beams].pol_amp = r3_normalize(get_double_array(token, &n));
        //--

        // Add main parameters
        beams[num_beams].s_0 = s_0;
        beams[num_beams].delta = delta;
        beams[num_beams].w = w;

        num_beams++;
    }
    //--

    c_beams = (beam_t *) calloc(num_beams, sizeof(beam_t));
    for(n = 0; n < num_beams; n++) c_beams[n] = beams[n];

    beams_setup.num = num_beams;
    beams_setup.beams = c_beams;
    //--

    //
    // Get sidebands
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/sidebands.csv");
    rows = read_lines(path);
    i = 0;

    //
    // Number of sidebands
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE)
        beams_setup.sidebands.num = atoi(str_replace(token, ".", ","));

    // C program
    else beams_setup.sidebands.num = atoi(token);
    //--

    //
    // Resonant frequency of the sidebands
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) beams_setup.sidebands.freq = atof(str_replace(token, ".", ","));

    // C program
    else beams_setup.sidebands.freq = atof(token);
    //--
    //--

    // Release memory
    free(beams);
    free(path);

    // Return
    return beams_setup;
}

transition_t get_transition(char *params_path){
    // Variables
    int i = 0;
    char *token, *saveptr, *path, **rows;
    transition_t transition;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "transition.csv");
    rows = read_lines(path);

    //
    // Transition rate
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.gamma = atof(str_replace(token, ".", ","));

    // C program
    else transition.gamma = atof(token);
    //--

    //
    // Resonant wave length
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    if(Py_MODULE) transition.lambda = atof(str_replace(token, ".", ","));

    // C program
    else transition.lambda = atof(token);
    //--

    //
    // Total angular momentum of the ground state
    //--    
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.J_gnd = (int) atof(str_replace(token, ".", ","));

    // C Program
    else transition.J_gnd = (int) atof(token);  
    //--

    //
    // Total angular momentum of the excited state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.J_exc = (int) atof(str_replace(token, ".", ","));

    // C Program
    else transition.J_exc = (int) atof(token);  
    //--

    //
    // Landè factor of the ground state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.g_gnd = atof(str_replace(token, ".", ","));

    // C Program
    else transition.g_gnd = atof(token);  
    //--

    //
    // Landè factor of the excited state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.g_exc = atof(str_replace(token, ".", ","));

    // C Program
    else transition.g_exc = atof(token); 
    //--

    //
    // Check values
    //--
    if((transition.J_exc - transition.J_gnd) < 0){
        printf("J_exc must be grater than J_gnd.\n");
        exit(0);
    } else if(transition.J_exc < 0 || transition.J_gnd < 0){
        printf("J_exc and J_gnd must be positive values.\n");
        exit(0);
    }
    //--

    //
    // Release memory
    //--
    free(path);
    //--

    // Return
    return transition;
}

atom_t get_atom(initial_conditions_t ini_conds, performance_t perform, beams_setup_t beams_setup, int opt, char *params_path){
    //
    // Variables
    // --
    int i = 0;
    char *path, **rows;
    char *token, *saveptr;
    double std_dev, *rd_v;
    atom_t atom;
    //--

    // Read lines from the CSV file
    path = str_concatenate(params_path, "atom.csv");
    rows = read_lines(path);

    //
    // Symbol
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    atom.symbol = (char *) malloc(strlen(token) * sizeof(char));
    strcpy(atom.symbol, token);
    //--

    //
    // Atomic number
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) atom.Z = (int) atof(str_replace(token, ".", ","));

    // C program
    else atom.Z = (int) atof(token);
    //--

    //
    // Mass
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) atom.mass = atof(str_replace(token, ".", ","));

    // C program
    else atom.mass = atof(token);
    //--

    //
    // Initial position
    //--
    if(opt < 2){
        //
        // Random vector
        //--
        rd_v = (double*) calloc(3, sizeof(double));

        for(i = 0; i < 3; i++) {
            rd_v[i] = ((double) rand()) / ((double) RAND_MAX);
            if((((double) rand()) / ((double) RAND_MAX)) < 0.5) rd_v[i] = -rd_v[i];
        }

        rd_v = r3_normalize(rd_v); // Normalization
        //--

        if(perform.max_r < beams_setup.beams[0].w) std_dev = (perform.max_r / 2);
        else std_dev = (beams_setup.beams[0].w / 2);
        
        atom.pos = r3_scalar_product(random_norm(0, std_dev), rd_v); // cm
    
    } else if(opt == 2)
        atom.pos = r3_scalar_product(-perform.max_r, r3_normalize(ini_conds.v_0_dir));

    //--

    //
    // Initial velocity of the atom
    //--
    if(opt < 2){
        atom.vel = (double *) calloc(3, sizeof(double));
        for(i = 0; i < 3; i++){
            std_dev = sqrt(k_B * ini_conds.T_0 / (atom.mass * u)) * 10; // cm / s
            atom.vel[i] = random_norm(0, std_dev); // cm / s
        }  

    } else if(opt == 2)
        atom.vel = r3_scalar_product(ini_conds.v_0, r3_normalize(ini_conds.v_0_dir));
    
    //--

    // Optical transition
    atom.transition = get_transition(params_path);

    // Initial state
    atom.J = atom.transition.J_gnd;
    atom.mJ = -atom.transition.J_gnd;

    // Release memory
    free(path);

    return atom;
}

int set_hist(int only_marginals, results_t *res, performance_t perform){
    // Variables
    int i, j, k;

    //
    // Only marginals
    //--
    if(only_marginals){
        res->pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));
        res->vel_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

        for(i = 0; i < 3; i++){
            // Position
            res->pos_hist[i].num_bins = perform.num_bins;
            res->pos_hist[i].bin_size = 2 * perform.max_r / res->pos_hist[i].num_bins;
            res->pos_hist[i].coord0 = - perform.max_r;
            res->pos_hist[i].freqs = (int*) calloc(res->pos_hist[i].num_bins, sizeof(int));

            // Velocity
            res->vel_hist[i].num_bins = perform.num_bins;
            res->vel_hist[i].bin_size = 2 * perform.max_v / res->vel_hist[i].num_bins;
            res->vel_hist[i].coord0 = - perform.max_v;
            res->vel_hist[i].freqs = (int*) calloc(res->vel_hist[i].num_bins, sizeof(int));
        }
    //--

    //
    // Complete histograms
    //--
    } else{
        // Position
        res->pos_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res->pos_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res->pos_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        // Velocity
        res->vel_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res->vel_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res->vel_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        for(i = 0; i < 3; i++){
            // Position
            res->pos_3Dhist.num_bins[i] = perform.num_bins;
            res->pos_3Dhist.bins_size[i] = 2 * perform.max_r / res->pos_3Dhist.num_bins[i];
            res->pos_3Dhist.coord0[i] = - perform.max_r;

            // Velocity
            res->vel_3Dhist.num_bins[i] = perform.num_bins;
            res->vel_3Dhist.bins_size[i] = 2 * perform.max_v / res->vel_3Dhist.num_bins[i];
            res->vel_3Dhist.coord0[i] = - perform.max_v;
        }

        //
        // Position
        //--
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
        //--

        //
        // Velocity
        //--
        res->vel_3Dhist.freqs = (int***) calloc(res->vel_3Dhist.num_bins[0], sizeof(int**));

        for(i = 0; i < res->vel_3Dhist.num_bins[0]; i++){
            res->vel_3Dhist.freqs[i] = (int**) calloc(res->vel_3Dhist.num_bins[1], sizeof(int*));
            for(j = 0; j < res->vel_3Dhist.num_bins[1]; j++){
                res->vel_3Dhist.freqs[i][j] = (int*) calloc(res->vel_3Dhist.num_bins[2], sizeof(int));
                for(k = 0; k < res->vel_3Dhist.num_bins[2]; k++){
                   res->vel_3Dhist.freqs[i][j][k] = 0; 
                }
            }
        }
        //--
    }
    //--

    return 0;
}

double move(atom_t *atom, beams_setup_t beams_setup, performance_t perform, magnetic_field_t B_params){
    //
    // Variables
    int i, j;                                   // Auxiliary variables
    double dt, aux_dt;                          // Time interval
    double R_max = -1, *R;                      // Scattering rates
    double *probs;                              // Probability to absorb a beam
    int chosen_beam = 0;                        // Absorption variables
    double *a_B;                                // Magnetic acceleration
    double *rd_v, vel_mod;                      // Photonic recoil

    // Allocate variables
    probs = (double*) calloc(beams_setup.num + 1, sizeof(double));

    // Get scattering rates
    R = get_scatt_rate(beams_setup, B_params, *atom);

    // Time interval
    //--
    aux_dt = perform.dt / (2 * PI * atom->transition.gamma*1e3);
    
    // Get minimum scattering rate
    for(i = 0; i < beams_setup.num; i++)
        if(R[i] > R_max || R_max < 0) R_max = R[i];

    if(aux_dt > (1 / (2*R_max))) dt = 1 / (4*R_max);
    else dt = aux_dt;
    //--

    // Get probabilities to absorb a beam
    for(i = 0; i < beams_setup.num; i++){
        probs[i + 1] = R[i] * dt;
        probs[0] += probs[i + 1];
    }

    // Probability of the atom does not absorb a beam
    probs[0] = 1 - probs[0];

    // Pick a beam
    chosen_beam = random_pick(probs, beams_setup.num + 1);
    //printf("chosen_beam = %d\n", chosen_beam);
    //if(chosen_beam > 0) r3_print(beams_setup.beams[chosen_beam-1].k_dir, "k");

    //
    // Movement
    //--

    // Magnetic acceleration
    a_B = magnetic_acceleration(*atom, B_params);

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
    if(perform.g_bool) 
        atom->pos[2] += -(g * dt*dt) / 2;
    // --

    //
    // Update velocity
    //--
    // Gravitational acceleration
    if(perform.g_bool) 
        atom->vel[2] += - g * dt;

    // Magnetic acceleration
    for(i = 0; i < 3; i++)
        atom->vel[i] += a_B[i] * dt;

    //
    // Photonic recoil
    //--
    if(chosen_beam > 0){
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
        //r3_print(r3_scalar_product(vel_mod, beams_setup.beams[chosen_beam-1].k_dir), "Photon momentum");
        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, beams_setup.beams[chosen_beam-1].k_dir)); // cm / s
        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, rd_v)); // cm / s

        // Release memory
        free(rd_v);
    }
    //--
    //--
    
    // Convert time unit
    dt = dt * (2 * PI * atom->transition.gamma*1e3);

    // Release memory
    free(probs);

    return dt;
}

double *magnetic_field(magnetic_field_t B_params, double *r){
    //
    // Variables
    //
    int i;
    double *B;          // Magnetic field vector
    double *r_prime;

    B = (double*) calloc(3, sizeof(double));
    r_prime = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++)
        r_prime[i] = r3_inner_product(r, B_params.B_basis[i]);

    // Anti-Helmholtz coils
    for(i = 0; i < 3; i++){
        B[i] = 0;

        // Linear magnetic field
        B[i] += r_prime[0] * B_params.B_basis[0][i] / 2;
        B[i] += r_prime[1] * B_params.B_basis[1][i] / 2;
        B[i] += - r_prime[2] * B_params.B_basis[2][i];
        B[i] = B_params.B_0 * B[i];

        // Linear magnetic field gradient
        B[i] += B_params.B_lin_grad[i]*r[i];
    }

    // Bias magnetic field
    B = r3_sum(B_params.B_bias, B);

    // Release memory
    free(r_prime);

    return B;
}

double *magnetic_acceleration(atom_t atom, magnetic_field_t B_params){
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
            r_prime[i] = r3_inner_product(atom.pos, B_params.B_basis[i]);

        norm = sqrt(pow(r_prime[0], 2) * pow(r_prime[1], 2) + 4 * pow(r_prime[2], 2));

        // Anti-helmholtz coils
        for(i = 0; i < 3; i++){
            del_B[i] = 0;
            del_B[i] += B_params.B_basis[0][i] * r_prime[0] / (2 * norm);
            del_B[i] += B_params.B_basis[1][i] * r_prime[1] / (2 * norm);
            del_B[i] += - B_params.B_basis[2][i] * 2  * r_prime[2] / norm;
            del_B[i] = B_params.B_0 * del_B[i];

            // Linear magnetic field gradient
            del_B[i] += B_params.B_lin_grad[i];
        }

        // Magnetic field Bias
        del_B = r3_sum(del_B, B_params.B_bias);

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

double *get_scatt_rate(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom){
    // Variables          
    int i, j, l, m;                                             // Auxiliary variables
    double *B, *eB;                                             // Magnetic field     
    double *R;                                                  // Scattering rates
    double s, s_0, r;                                           // Saturation parameter
    double zeeman_shift, doppler_shift, laser_detuning, delta;  // Detuning
    double lambda, gamma, g_gnd, g_exc;                         // Transition
    int mj_gnd, mj_exc;                                         // Transition
    double **C;                                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                                    // All polarizations
    beam_t beam;                                                // Beam     

    // Allocate memory
    R = (double*) calloc(beams_setup.num, sizeof(double));

    //
    // Magnetic field
    //--
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0){
        eB = B_params.B_basis[2];
    } else eB = r3_normalize(B);
    //--

    //
    // Check each beam
    //--
    for(i = 0; i < beams_setup.num; i++){
        // Get beam
        beam = beams_setup.beams[i];
        //r3_print(beam.k_dir, "k");

        // Polarizations
        set_polarizations_amplitudes(&beam, eB);

        //
        // Initial Saturation parameter
        //--
        // Basis of the beam frame
        C = orthonormal_basis(r3_normalize(beam.k_dir));

        // Distance from the propagation axis
        r = pow(r3_inner_product(C[0], atom.pos), 2);
        r += pow(r3_inner_product(C[1], atom.pos), 2);

        s_0 = beam.s_0;
        s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
        //--

        // Transition
        gamma = atom.transition.gamma; // kHz / 2pi
        lambda = atom.transition.lambda; // nm

        // Doppler shift
        doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz / 2pi

        //
        // Check all possible transitions
        //--
        // Polarizations
        for(j = 0; j < 3; j++){
            //
            // Zeeman shift
            //--
            // Ground state
            mj_gnd = - atom.transition.J_gnd;
            
            // Excited state
            mj_exc = mj_gnd + pol[j];

            // Landè factors
            g_gnd = atom.transition.g_gnd;
            g_exc = atom.transition.g_exc;

            // Compute shift
            zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz / 2pi
            //--

            // Saturation parameter considering sidebands
            s = beam.pol_amp[j] * s_0 / (2*beams_setup.sidebands.num + 1);

            // Main beam
            laser_detuning = beam.delta*gamma; // kHz / 2pi
            delta = laser_detuning + zeeman_shift + doppler_shift;
            R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

            // Get all scattering rates due to the sidebands
            for(m = 0; m < beams_setup.sidebands.num; m++){
                // Right sideband
                laser_detuning = (beam.delta + (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

                // Left sideband
                laser_detuning = (beam.delta - (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            }
        }

        // Release memory  
        for(l = 0; l < 3; l++) free(C[l]); 
        free(C);

        //printf("R[%d] = %f\n\n", i + 1, R[i]);
    }
    //--

    return R;    
}

int set_polarizations_amplitudes(beam_t *beam, double *eB){
    //
    // Variables
    int i, j;
    double **R1, **R2;
    double *pol_amp;
    complex_t **C1, **C2, **A_s_c, **A;
    complex_t *C1_eps, *C2_eps, *C2_s_eps;

    // Bases of the beam frame
    R1 = orthonormal_basis(beam->k_dir);
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
    C1_eps = c3_apply_operator(A_s_c, r3_to_c3(beam->pol_amp));

    // Polarization on the Cartesian B frame
    C2_eps = c3_apply_operator(A, C1_eps);

    // Polarization on the Spherical B frame
    C2_s_eps = c3_apply_operator(c3_operator_dagger(A_s_c), C2_eps);

    // Set Polarization amplitudes
    pol_amp = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) pol_amp[i] = pow(c_mod(C2_s_eps[i]), 2);
    beam->pol_amp = pol_amp;
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

    return 1;
}

/*

double *get_probs(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom, double dt){
    // Variables          
    int i, j, l, m;                                             // Auxiliary variables
    double *B, *eB;                                             // Magnetic field     
    double *probs;                                              // Probabilities
    double s, s_0, r;                                           // Saturation parameter
    double zeeman_shift, doppler_shift, laser_detuning, delta;  // Detuning
    double lambda, gamma, g_gnd, g_exc;                         // Transition
    int mj_gnd, mj_exc;                                         // Transition
    double **C;                                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                                    // All polarizations
    beam_t beam;                                                // Beam     

    // Allocate memory
    probs = (double*) calloc(beams_setup.num + 1, sizeof(double));

    //
    // Magnetic field
    //--
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0){
        eB = B_params.B_basis[2];
    } else eB = r3_normalize(B);
    //--

    //
    // Check each beam
    //--
    for(i = 0; i < beams_setup.num; i++){
        // Get beam
        beam = beams_setup.beams[i];
        //r3_print(beam.k_dir, "k");

        // Polarizations
        set_polarizations_amplitudes(&beam, eB);

        //
        // Initial Saturation parameter
        //--
        // Basis of the beam frame
        C = orthonormal_basis(r3_normalize(beam.k_dir));

        // Distance from the propagation axis
        r = pow(r3_inner_product(C[0], atom.pos), 2);
        r += pow(r3_inner_product(C[1], atom.pos), 2);

        s_0 = beam.s_0;
        s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
        //--

        // Transition
        gamma = atom.transition.gamma; // kHz / 2pi
        lambda = atom.transition.lambda; // nm

        // Doppler shift
        doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz / 2pi

        //
        // Check all possible transitions
        //--
        // Polarizations
        for(j = 0; j < 3; j++){
            //
            // Zeeman shift
            //--
            // Ground state
            mj_gnd = - atom.transition.J_gnd;
            
            // Excited state
            mj_exc = mj_gnd + pol[j];

            // Landè factors
            g_gnd = atom.transition.g_gnd;
            g_exc = atom.transition.g_exc;

            // Compute shift
            zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz / 2pi
            //--

            // Saturation parameter considering sidebands
            s = beam.pol_amp[j] * s_0 / (2*beams_setup.sidebands.num + 1);

            // Main beam
            laser_detuning = beam.delta*gamma; // kHz / 2pi
            delta = laser_detuning + zeeman_shift + doppler_shift;
            probs[i+1] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

            // Get all scattering rates due to the sidebands
            for(m = 0; m < beams_setup.sidebands.num; m++){
                // Right sideband
                laser_detuning = (beam.delta + (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                probs[i+1] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

                // Left sideband
                laser_detuning = (beam.delta - (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                probs[i+1] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            }
        }

        // Release memory  
        for(l = 0; l < 3; l++) free(C[l]); 
        free(C);

        probs[i+1] = probs[i+1] * dt;
        probs[0] += probs[i+1];

        //printf("prob[%d] = %f\n\n", i + 1, probs[i+1]);
    }

    // Probability of do not absorb a beam
    probs[0] = 1 - probs[0];
    //printf("\nprob[0] = %f\n", probs[0]);

    //--

    return probs;    
}

double *get_all_scatt_rate(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom){
    // Variables          
    int i, j, k = 0, l, m;                                      // Auxiliary variables
    double *B, *eB;                                             // Magnetic field     
    double *R;                                                  // Scattering rates
    double s, s_0, r;                                           // Saturation parameter
    double zeeman_shift, doppler_shift, laser_detuning, delta;  // Detuning
    double lambda, gamma, g_gnd, g_exc;                         // Transition
    int mj_gnd, mj_exc;                                         // Transition
    double **C;                                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                                    // All polarizations
    beam_t beam;                                                // Beam     

    // Allocate memory
    //R = (double*) calloc(beams_setup.num*3*(2*beams_setup.sidebands.num + 1), sizeof(double));
    R = (double*) calloc(beams_setup.num, sizeof(double));

    //
    // Magnetic field
    //--
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0){
        eB = B_params.B_basis[2];
    } else eB = r3_normalize(B);
    //--

    //
    // Check each beam
    //--
    for(i = 0; i < beams_setup.num; i++){
        // Get beam
        beam = beams_setup.beams[i];
        //r3_print(beam.k_dir, "k");

        // Polarizations
        set_polarizations_amplitudes(&beam, eB);

        //
        // Initial Saturation parameter
        //--
        // Basis of the beam frame
        C = orthonormal_basis(r3_normalize(beam.k_dir));

        // Distance from the propagation axis
        r = pow(r3_inner_product(C[0], atom.pos), 2);
        r += pow(r3_inner_product(C[1], atom.pos), 2);

        s_0 = beam.s_0;
        s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
        //--

        // Transition
        gamma = atom.transition.gamma; // kHz / 2pi
        lambda = atom.transition.lambda; // nm

        // Doppler shift
        doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz / 2pi
        //printf("doppler_shift = %f\n", doppler_shift);

        //
        // Check all possible transitions
        //--
        // Polarizations
        for(j = 0; j < 3; j++){
            //
            // Zeeman shift
            //--
            // Ground state
            mj_gnd = - atom.transition.J_gnd;
            
            // Excited state
            mj_exc = mj_gnd + pol[j];

            // Landè factors
            g_gnd = atom.transition.g_gnd;
            g_exc = atom.transition.g_exc;

            // Compute shift
            //printf("B = %f\n", r3_mod(B));
            zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz / 2pi
            //--

            // Saturation parameter considering sidebands
            s = beam.pol_amp[j] * s_0 / (2*beams_setup.sidebands.num + 1);
            //printf("s = %f\n", s);

            // Main beam
            laser_detuning = beam.delta*gamma; // kHz / 2pi
            delta = laser_detuning + zeeman_shift + doppler_shift;
            R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

            //printf("pol = %d (eps = %f)\n", pol[j], beam.pol_amp[j]);
            //printf("sideband %d\n", 0);
            //printf("Laser detuning [kHz / 2pi] = %f\n", laser_detuning);
            //printf("Doppler shift [kHz / 2pi] = %f\n", doppler_shift);
            //printf("Zeeman shift [kHz / 2pi] = %f\n", zeeman_shift);
            //printf("|B| = %f\n", r3_mod(B));
            //printf("R[%d] = %f\n\n", k, R[k]);

            //k += 1;

            // Get all scattering rates due to the sidebands
            for(m = 0; m < beams_setup.sidebands.num; m++){
                // Right sideband
                laser_detuning = (beam.delta + (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

                //r3_print(beam.k_dir, "k");
                //printf("pol = %d (eps = %f)\n", pol[j], beam.pol_amp[j]);
                //printf("sideband %d (Right)\n", m+1);
                //printf("Laser detuning [kHz / 2pi] = %f\n", laser_detuning);
                //printf("Doppler shift [kHz / 2pi] = %f\n", doppler_shift);
                //printf("Zeeman shift [kHz / 2pi] = %f\n", zeeman_shift);
                //printf("R[%d] = %f\n\n", k, R[k]);

                //k += 1;

                // Left sideband
                laser_detuning = (beam.delta - (m+1)*beams_setup.sidebands.freq)*gamma; // kHz / 2pi
                delta = laser_detuning + zeeman_shift + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

                //r3_print(beam.k_dir, "k");
                //printf("pol = %d (eps = %f)\n", pol[j], beam.pol_amp[j]);
                //printf("sideband %d (Left)\n", m+1);
                //printf("Laser detuning [kHz / 2pi] = %f\n", laser_detuning);
                //printf("Doppler shift [kHz / 2pi] = %f\n", doppler_shift);
                //printf("Zeeman shift [kHz / 2pi] = %f\n", zeeman_shift);
                //printf("R[%d] = %f\n\n", k, R[k]);

                //k += 1;
            }
        }

        // Release memory  
        for(l = 0; l < 3; l++) free(C[l]); 
        free(C);

        //printf("R[%d] = %f\n\n", i+1, R[i]);

    }
    //--

    //exit(0);

    return R;
}

double *get_scatt_rates(beam_t beam, sidebands_t sidebands, double *B, atom_t atom){
    //
    // Variables
    double *R;                                  // Scattering rates
    double s, s_0, r;                           // Saturation parameter
    double zeeman_shift, doppler_shift, delta;  // Detuning
    double lambda, gamma, g_gnd, g_exc;         // Transition
    int mj_gnd, mj_exc;                         // Transition
    double **C;                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                    // All polarizations
    int i, j, k;

    // Allocate memory
    R = (double*) calloc(3*(2*sidebands.num + 1), sizeof(double));

    //
    // Initial Saturation parameter
    //--
    // Basis of the beam frame
    C = orthonormal_basis(r3_normalize(beam.k_dir));

    // Distance from the propagation axis
    r = pow(r3_inner_product(C[0], atom.pos), 2);
    r += pow(r3_inner_product(C[1], atom.pos), 2);

    s_0 = beam.s_0;
    s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
    //--

    // Transition
    gamma = atom.transition.gamma; // kHz / 2pi
    lambda = atom.transition.lambda; // nm

    // Doppler shift
    doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz / 2pi
    //printf("doppler_shift = %f\n", doppler_shift);

    //
    // Check all possible transitions
    //--
    // Indexes
    k = 0;

    // Polarizations
    for(i = 0; i < 3; i++){
        //
        // Zeeman shift
        //--
        // Ground state
        mj_gnd = - atom.transition.J_gnd;
        
        // Excited state
        mj_exc = mj_gnd + pol[i];

        // Landè factors
        g_gnd = atom.transition.g_gnd;
        g_exc = atom.transition.g_exc;

        // Compute shift
        //printf("B = %f\n", r3_mod(B));
        zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz / 2pi
        //--

        // Saturation parameter considering sidebands
        s = beam.pol_amp[i] * s_0 / (2*sidebands.num + 1);
        //printf("s = %f\n", s);

        // Main beam
        delta = beam.delta*gamma; // kHz / 2pi
        delta += zeeman_shift + doppler_shift;
        //r3_print(beam.k_dir, "k");
        //r3_print(beam.pol_amp, "eps");
        //printf("Laser detuning [kHz / 2pi] = %f\n", beam.delta*gamma);
        //printf("Doppler shift [kHz / 2pi] = %f\n", doppler_shift);
        //printf("Zeeman shift [kHz / 2pi] = %f\n\n", zeeman_shift);
        R[k] = ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
        printf("R = %f\n\n", R[k]);
        k += 1;

        // Get all scattering rates due to the sidebands
        for(j = 0; j < sidebands.num; j++){
            // Right sideband
            delta = beam.delta + (j+1)*sidebands.freq; // kHz / 2pi
            delta += zeeman_shift + doppler_shift;
            R[k] = ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            printf("R = %f\n\n", R[k]);
            k += 1;

             // Left sideband
            delta = beam.delta - (j+1)*sidebands.freq; // kHz / 2pi
            delta += zeeman_shift + doppler_shift;
            R[k] = ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            printf("R = %f\n\n", R[k]);
            k += 1;
        }
    }

    // Release memory  
    for(i = 0; i < 3; i++) free(C[i]); 
    free(C);

    // Return
    return R;
}

double get_total_scatt_rate(beam_t beam, sidebands_t sidebands, double *B, atom_t atom){
    //
    // Variables
    double R = 0;                               // Total Scattering rate
    double s, s_0, r;                           // Saturation parameter
    double zeeman_shift, doppler_shift, delta;  // Detuning
    double lambda, gamma, g_gnd, g_exc;         // Transition
    int mj_gnd, mj_exc;                         // Transition
    double **C;                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                    // All polarizations
    int i, j, k;

    //
    // Initial Saturation parameter
    //--
    // Basis of the beam frame
    C = orthonormal_basis(r3_normalize(beam.k_dir));

    // Distance from the propagation axis
    r = pow(r3_inner_product(C[0], atom.pos), 2);
    r += pow(r3_inner_product(C[1], atom.pos), 2);

    s_0 = beam.s_0;
    s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
    //--

    // Transition
    gamma = atom.transition.gamma; // kHz / 2pi
    lambda = atom.transition.lambda; // nm

    // Doppler shift
    doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // kHz / 2pi
    //printf("doppler_shift = %f\n", doppler_shift);

    //
    // Check all possible transitions
    //--
    // Polarizations
    //r3_print(beam.pol_amp, "pol_amp");
    for(i = 0; i < 3; i++){
        //
        // Zeeman shift
        //--
        // Ground state
        mj_gnd = - atom.transition.J_gnd;
        
        // Excited state
        mj_exc = mj_gnd + pol[i];

        // Landè factors
        g_gnd = atom.transition.g_gnd;
        g_exc = atom.transition.g_exc;

        // Compute shift
        //printf("B = %f\n", r3_mod(B));
        zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * (g_gnd * mj_gnd - g_exc * mj_exc);  // kHz / 2pi
        //--

        // Saturation parameter considering sidebands
        s = beam.pol_amp[i] * s_0 / (2*sidebands.num + 1);
        //printf("s = %f\n", s);

        // Main beam
        delta = beam.delta*gamma; // kHz / 2pi
        delta += zeeman_shift + doppler_shift;
        //r3_print(beam.k_dir, "k");
        //r3_print(beam.pol_amp, "eps");
        //printf("Laser detuning [kHz / 2pi] = %f\n", beam.delta*gamma);
        //printf("Doppler shift [kHz / 2pi] = %f\n", doppler_shift);
        //printf("Zeeman shift [kHz / 2pi] = %f\n\n", zeeman_shift);
        R += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
        printf("R = %f\n\n", R);

        // Get all scattering rates due to the sidebands
        for(j = 1; j < (sidebands.num + 1); j++){

            for(k = -1; k < 2; k+=2){
                // Detuning
                delta = beam.delta + i*k*sidebands.freq; // kHz / 2pi
                delta += zeeman_shift + doppler_shift;

                R += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            }
        }
    }

    // Release memory  
    for(i = 0; i < 3; i++) free(C[i]); 
    free(C);

    // Return
    return R;
}

*/

//
// Utility functions
//

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

char **read_lines(char *path){
    //
    // Variables
    int i = 0, j = 0, k = 0;
    FILE *fp;
    char *row, *aux_row;
    char **rows;

    // Alocate memory
    rows = (char**) calloc(MAX_LINES, sizeof(char*));

    // Open file
    if((fp = fopen(path, "r"))){
        // Read lines
        while(!feof(fp)){
            // Allocate memory
            aux_row = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));
            row = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));

            // Try to read line
            if(fgets(row, STRING_BUFFER_SIZE, fp) != NULL){
                // Get row length and remove \n
                k = 0;
                for(j = 0; row[j] != '\0'; j++) {
                    if(row[j] != '\n'){
                        aux_row[k] = row[j];
                        k += 1;
                    }
                }

                // Copy the read row
                rows[i] = (char*) calloc(k+1, sizeof(char));
                strcpy(rows[i], aux_row);

                i += 1;
            } else break;

            // Release memory
            free(row);
            free(aux_row);
        }

        // Indicate the file ending
        rows[i+1] = (char*) calloc(3, sizeof(char));
        rows[i+1] = NULL;

        // Close file
        fclose(fp); 
    } else {
        printf("File \"%s\" does not exist\n", path);
        exit(0);
    }

    return rows;
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
        // Py MODULE
        if(Py_MODULE){
            if(token[0] == '-') aux_arr[i] = -atof(str_replace(token+1, ".", ",")); 
            else aux_arr[i] = atof(str_replace(token, ".", ","));
        }

        // C Program
        else {
            if(token[0] == '-') aux_arr[i] = -atof(token+1); 
            else aux_arr[i] = atof(token);
        }

        token = strtok_r(str, " ", &str);
        i++;
    }

    arr = (double *) malloc(i * sizeof(double));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    // Release memory
    free(token);

    return arr;
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

double random_exp(double mean){
    return (- mean * log(1 - ((double) rand()) / ((double) RAND_MAX)));
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

int update_hist(histogram_t *hist, double val){
    //
    // Variables
    //

    int bin;
    double lower_lim, upper_lim;

    // Add frequency
    for(bin = 0; bin < hist->num_bins; bin++){
        lower_lim = hist->coord0 + bin * hist->bin_size;
        upper_lim = lower_lim + hist->bin_size;

        if((val >= lower_lim) && (val < upper_lim)){
            hist->freqs[bin] += 1;
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

//
// Debug

int print_performance(performance_t perform){
    printf("Performance\n--\n");
    printf("max_time [1/Gamma] = %f\n", perform.max_time);
    printf("max_r [cm] = %f\n", perform.max_r);
    printf("max_v [cm/s] = %f\n", perform.max_v);
    printf("num_bins = %d\n", perform.num_bins);
    printf("wait_time [1/Gamma] = %f\n", perform.wait_time);
    printf("time interval [1/Gamma] = %f\n", perform.dt);
    printf("\n");

    return 1;
}

int print_initial_conditions(initial_conditions_t ini_conds){
    printf("Initial Conditions\n--\n");
    printf("T_0 [uK] = %f\n", ini_conds.T_0);
    printf("v_0 [m/s] = %f\n", ini_conds.v_0);
    r3_print(ini_conds.v_0_dir, "v_0_dir");
    printf("g_bool = %d\n", ini_conds.g_bool);
    printf("\n");

    return 1;
}

int print_magnetic_field(magnetic_field_t B_params){
    printf("Magnetic field\n--\n");
    printf("B_0 [G/cm] = %f\n", B_params.B_0);
    r3_operator_print(B_params.B_basis, "B_basis");
    r3_print(B_params.B_bias, "B_bias [G]");
    printf("\n");

    return 1;
}

int print_atom(atom_t atom){
    printf("Atom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass [u] = %f\n", atom.mass);
    r3_print(atom.pos, "pos [cm/s]");
    printf("r [cm] = %f\n", r3_mod(atom.pos));
    r3_print(atom.vel, "vel [cm/s]");
    printf("v [cm/s] = %f\n", r3_mod(atom.vel));
    printf("J = %d\n", atom.J);
    printf("mJ = %d\n", atom.mJ);
    printf("\n");

    return 1;
}

int print_beams(beams_setup_t beams_setup){
    // Variables
    int i;

    printf("\nSidebands\n");
    printf("num = %d\n", beams_setup.sidebands.num);
    printf("freq = %f\n", beams_setup.sidebands.freq);
    printf("\n");
    printf("All Beams\n--\n");
    for(i = 0; i < beams_setup.num; i++){
        printf("Beam %d\n", i+1);
        printf("s_0 = %f\n", beams_setup.beams[i].s_0);
        printf("delta [1/Gamma] = %f\n", beams_setup.beams[i].delta);
        printf("w [cm] = %f\n", beams_setup.beams[i].w);
        r3_print(beams_setup.beams[i].k_dir, "k");
        r3_print(beams_setup.beams[i].pol_amp, "pol_amp");
        printf("--\n");
    }
    printf("\n");

    return 1;
}

int print_results(results_t res, atom_t atom, int opt){
    // Variables
    int i, j;

    printf("Simulation status\n--\n");
    printf("total time [ms] = %f\n", res.time / (atom.transition.gamma));
    printf("total time [1/Gamma] = %f\n", res.time);
    printf("Atom trapped [ms] = %d\n", res.trapped_atom);
    printf("\n");

    if(opt == 1){
        //
        // Position frequencies
        //-
        for(i = 0; i < 3;i++){
            printf("pos[%d] = [\n", i+1);
            for(j = 0; j < res.pos_hist[i].num_bins; j++)
                printf("%d ", res.pos_hist[i].freqs[j]);
            printf("]\n\n");
        }
        printf("\n");
        //--

        //
        // Velocity frequencies
        //--
        for(i = 0; i < 3;i++){
            printf("vel[%d] = [\n", i+1);
            for(j = 0; j < res.vel_hist[i].num_bins; j++)
                printf("%d ", res.vel_hist[i].freqs[j]);
            printf("]\n\n");
        }
        printf("\n");
        //--
    }

    return 1;
}

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


/**

int print_params(atom_t atom, performance_t perform, beams_setup_t beams_setup, environment_t env){
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
    printf("T_0 [uK] = %f\n", perform.T_0);
    printf("max_time [ms] = %f\n", perform.max_time / (atom.transition.gamma));
    printf("max_r [cm] = %f\n", perform.max_r);
    printf("max_v [cm/s] = %f\n", perform.max_v);
    printf("num_bins = %d\n", perform.num_bins);
    printf("wait_time [ms] = %f\n", perform.wait_time / (atom.transition.gamma));
    printf("time interval [ms] = %f\n", perform.dt / (atom.transition.gamma));
    printf("\n");

    // Environment
    printf("Environment\n--\n");
    printf("B_0 [G/cm] = %f\n", env.B_0);
    printf("local_B [G/cm] = %f\n", env.local_B);
    printf("delta [1/gamma] = %f\n", env.delta);
    printf("s_0 = %f\n", env.s_0);
    printf("w [cm] = %f\n", env.w);
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
**/