
//  Header

#include "mot_sim.h"

// --

// Execute simulation for a single atom
int C_run(){
    //
    // Variables
    //

    // Simple variables
    int i;

    // Parameters of the simulation
    int *active_params_ids = C_get_active_params();
    transition_t transition = C_get_transition(active_params_ids[0]);
    environment_t environment = C_get_environment(active_params_ids[1]);
    initial_conditions_t initial_conditions = C_get_initial_conditions(active_params_ids[2]);
    settings_t settings = C_get_settings(active_params_ids[3]);

    //
    //
    //

    return 0;
}

// Get transition
transition_t C_get_transition(int id){
    //
    // Variables
    //
    transition_t transition;
    int row_cter = 0, num_row_header = 12;
    char *rest;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/transitions.csv";
    char path[] = "../parameters/transitions.csv";
    char *token;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    // Get transition
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            if(atoi(token) == id){
                //
                // Atom (atom_id)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid atom_id in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);
                
                } else transition.atom = C_get_atom(atoi(token));

                //
                // Transition rate (gamma)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid gamma in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.gamma = atof(token);

                //
                // Resonant wave length (lambda)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid lambda in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.lambda = atof(token);

                //
                // Total angular momentum of the ground state (J_gnd)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(token == NULL){
                    printf("Invalid J_gnd in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.J_gnd = atoi(token);

                //
                // Total angular momentum of the excited state (J_exc)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid J_exc in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.J_exc = atoi(token);
                
                //
                // Landè factor of the ground state (g_gnd)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid g_gnd in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.g_gnd = atof(token);
                
                //
                // Landè factor of the excited state (g_exc)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid g_exc in the transition id=%d in the file \"%s\"!\n", id, path);
                    exit(0);

                } else transition.g_exc = atof(token);

                // Stop looping
                break;
            }  
        }        

        row_cter++;
    }

    // Check if there are added transition in the file
    if(row_cter < num_row_header) {
        printf("Entry with transitions in the file \"%s\"!\n", path);
        exit(0);
    }

    return transition;
}

// Get atom
atom_t C_get_atom(int id){
    //
    // Variables
    //

    atom_t atom;
    int row_cter = 0, num_row_header = 9;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/atoms.csv";
    char path[] = "../parameters/atoms.csv";
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    // Get transition
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            if(atoi(token) == id){
                //
                // Atom symbol (symbol)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid symbol in the atom id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else atom.symbol = token;

                //
                // Atom name (name)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid name in the atom id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } 

                //
                // Atomic number (Z)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid Z in the atom id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else atom.Z = atoi(token);

                //
                // Atom mass (mass)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid mass in the atom id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else atom.mass = atof(token);

                // Stop looping
                break;
            }  
        }        
        
        row_cter++;
    }

    // Check if there are added atoms in the file
    if(row_cter < num_row_header) {
        printf("Entry with atoms in the file \"%s\"\n", path);
        exit(0);
    }

    return atom;
}

// Get environment
environment_t C_get_environment(int id){
    //
    // Variables
    //

    environment_t environment;
    int *beams_id = (int*) malloc(sizeof(int));
    int row_cter = 0, num_row_header = 8, num_beams, i;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/environments.csv";
    char path[] = "../parameters/environments.csv";
    char *token, *token2, *rest, *str_beams_id;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    // Get transition
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            if(atoi(token) == id){
                //
                // Name
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid name in the environment id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                }

                //
                // Beams setup
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid bemas_setup in the environment id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else {
                    // Get IDs of the beams    
                    num_beams = C_get_int_array(token, MAX_BEAMS, beams_id);

                    // Get beams from CSV file
                    environment.beams_setup = C_get_beams_setup(num_beams, beams_id);
                }

                // Stop looping
                break;
            }  
        }        
        
        row_cter++;
    }

    // Check if there are added atoms in the file
    if(row_cter < num_row_header) {
        printf("Entry with atoms in the file \"%s\"\n", path);
        exit(0);
    }

    return environment;
}

// Get beams setup
beams_setup_t C_get_beams_setup(int num_beams, int *beams_id){
    //
    // Variables
    //

    beams_setup_t beams_setup;
    beam_t beams[num_beams];

    int row_cter = 0, num_row_header = 10, i, beam_idx;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/beams.csv";
    char path[] = "../parameters/beams.csv";
    char *token, *rest;
    FILE *fp;

    // Number of beams
    beams_setup.num = num_beams;
    beams_setup.beams = (beam_t*) malloc(num_beams * sizeof(beam_t));

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    //
    // Get beams
    //

    beam_idx = 0;

    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL) && (beam_idx < num_beams)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            for(i = 0; i < num_beams; i++){
                if(atoi(token) == beams_id[i]){
                    //
                    // Detuning (delta)
                    //

                    token = strtok_r(rest, DELIM, &rest);

                    if(!token){
                        printf("Invalid delta in the beam id=%d in the file \"%s\"\n", beams_id[i], path);
                        exit(0);
                    } else {
                        if(token[0] == '-') beams_setup.beams[beam_idx].delta = -atof((token + 1)); 
                        else beams_setup.beams[beam_idx].delta = atof(token);
                    }

                    //
                    // Wave vector direction (k_dic)
                    //

                    token = strtok_r(rest, DELIM, &rest);

                    if(!token){
                        printf("Invalid k_dic in the beam id=%d in the file \"%s\"\n", beams_id[i], path);
                        exit(0);
                    } else {
                        beams_setup.beams[beam_idx].k_dic = (float*) malloc(sizeof(float));
                        C_get_float_array(token, 2, beams_setup.beams[beam_idx].k_dic);
                    }

                    //
                    // Polarization vector (eps)
                    //

                    token = strtok_r(rest, DELIM, &rest);

                    if(!token){
                        printf("Invalid eps in the beam id=%d in the file \"%s\"\n", beams_id[i], path);
                        exit(0);
                    } else {
                        beams_setup.beams[beam_idx].eps = (float*) malloc(sizeof(float));
                        C_get_float_array(token, 3, beams_setup.beams[beam_idx].eps);
                    }

                    //
                    // Peak of the saturation parameter (s_0)
                    //

                    token = strtok_r(rest, DELIM, &rest);

                    if(!token){
                        printf("Invalid s_0 in the beam id=%d in the file \"%s\"\n", beams_id[i], path);
                        exit(0);
                    } else {
                        if(token[0] == '-') beams_setup.beams[beam_idx].s_0 = -atof((token + 1)); 
                        else beams_setup.beams[beam_idx].s_0 = atof(token);
                    }

                    //
                    // Waist radius (w)
                    //

                    token = strtok_r(rest, DELIM, &rest);

                    if(!token){
                        printf("Invalid w in the beam id=%d in the file \"%s\"\n", beams_id[i], path);
                        exit(0);
                    } else {
                        if(token[0] == '-') beams_setup.beams[beam_idx].w = -atof((token + 1)); 
                        else beams_setup.beams[beam_idx].w = atof(token);
                    }

                    // Stop looping
                    beam_idx++;
                    break;
                }
            }  
        }      
        
        row_cter++;
    }

    // Check if there are added atoms in the file
    if(row_cter < num_row_header) {
        printf("Entry with atoms in the file \"%s\"\n", path);
        exit(0);
    }

    // Check if was found all beams
    if(beam_idx != num_beams) {
        printf("It was not found all beams in the file \"%s\"\n", path);
        exit(0);
    }

    return beams_setup;
}

// Get initial conditions
initial_conditions_t C_get_initial_conditions(int id){
    //
    // Variables
    //

    initial_conditions_t initial_conditions;
    int row_cter = 0, num_row_header = 7;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/initial_conditions.csv";
    char path[] = "../parameters/initial_conditions.csv";
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    // Get initial conditions
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            if(atoi(token) == id){
                //
                // Initial temperature (T_0)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid T_0 in the initial conditions set with id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else initial_conditions.T_0 = atoi(token);

                //
                // Gravity (g_bool)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid g_bool in the initial conditions set with id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else initial_conditions.g_bool = atoi(token);

                // Stop looping
                break;
            }  
        }        
        
        row_cter++;
    }

    // Check if there are added atoms in the file
    if(row_cter < num_row_header) {
        printf("Entry with initial conditions set in the file \"%s\"\n", path);
        exit(0);
    }

    return initial_conditions;
}

// Get settings
settings_t C_get_settings(int id){
    //
    // Variables
    //

    settings_t settings;
    int row_cter = 0, num_row_header = 9;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/settings.csv";
    char path[] = "../parameters/settings.csv";
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\".\n", path);
        exit(0);
    }

    // Get initial conditions
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest);

        // Skip head
        if(row_cter >= num_row_header){
            // Check line
            if(atoi(token) == id){
                //
                // Maximum number of iteration (i_max)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid i_max in the settings id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else settings.i_max = (int) atof(token);

                //
                // Maximum distance (r_max)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid r_max in the settings id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else settings.r_max = atof(token);

                //
                // Number of simulations (num_sim)
                //

                token = strtok_r(rest, DELIM, &rest);
                if(!token){
                    printf("Invalid num_sim in the settings id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                }

                //
                // Number of bins in the histogram (num_bins)
                //

                token = strtok_r(rest, DELIM, &rest);

                if(!token){
                    printf("Invalid num_bins in the settings id=%d in the file \"%s\"\n", id, path);
                    exit(0);
                } else settings.num_bins = atoi(token);

                // Stop looping
                break;
            }  
        }        
        
        row_cter++;
    }

    // Check if there are added atoms in the file
    if(row_cter < num_row_header) {
        printf("Entry with initial conditions set in the file \"%s\"\n", path);
        exit(0);
    }

    return settings;
}

// Get active parameters
int * C_get_active_params(){
    //
    // Variables
    //

    static int ids[NUM_ACTIVE_PARAMS];
    int row_cter = 0;
    char row[STRING_BUFFER_SIZE];
    //char path[] = "model/parameters/active.csv";
    char path[] = "../parameters/active.csv";
    char *token;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the the file \"%s\".\n", path);
        exit(0);
    }

    // Read first line (header)
    if(fgets(row, STRING_BUFFER_SIZE, fp) == NULL) {
        printf("The file \"%s\" is empty.\n", path);
        exit(0);
    }

    // Read second line (active parameters IDs)
    if(fgets(row, STRING_BUFFER_SIZE, fp) == NULL) {
        printf("Entry with active parameters in the file \"%s\"!\n", path);
        exit(0);

    // Parse row
    } else {
        token = strtok(row, DELIM);

        ids[0] = atoi(token);
        row_cter = 1;

        while(!(token == NULL) && (row_cter < NUM_ACTIVE_PARAMS)){
            ids[row_cter] = atoi(token);
            row_cter++;

            token = strtok(NULL, DELIM);
        }        
    }

    return ids;
}

// Get int array from string in the format [i1 i2 ... in]
int C_get_int_array(char *str, int max_size, int *arr){
    //
    // Variables
    //
    int i, j;
    char *token;
    int aux_arr[max_size];

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

    free(arr);
    arr = (int *) malloc(i * sizeof(int));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    return i;
}

// Get float array from string in the format [f1 f2 ... fn]
int C_get_float_array(char *str, int max_size, float *arr){
    //
    // Variables
    //
    int i, j;
    char *token;
    float aux_arr[max_size];

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

    free(arr);
    arr = (float *) malloc(i * sizeof(float));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    return i;
}
