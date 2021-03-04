
//  Header

#include "mot_sim.h"

// --

// Execute simulation for a single atom
int C_run(){
    //
    // Variables
    //

    // Simple variables
    //char *results_path = concatenate_ROOT_PATH("results/");
    int i;

    // Parameters of the simulation
    atom_t atom = C_get_atom();
    transition_t transition = C_get_transition();
    conditions_t conditions = C_get_conditions();
    beams_setup_t beams_setup = C_get_beams();
    constants_t cts = C_get_constants();

    //
    // Show all parameters
    //

    /*
    printf("\nAtom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass = %f\n", atom.mass);

    printf("\nTransition\n--\n");
    printf("gamma = %f\n", transition.gamma);
    printf("J_gnd = %d\n", transition.J_gnd);
    printf("J_exc = %d\n", transition.J_exc);
    printf("g_gnd = %f\n", transition.g_gnd);
    printf("g_exc = %f\n", transition.g_exc);

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
        printf("k_dic = [%f, %f]\n", beams_setup.beams[i].k_dic[0], beams_setup.beams[i].k_dic[1]);
        printf("eps = [%f, %f, %f]\n", beams_setup.beams[i].eps[0], beams_setup.beams[i].eps[1], beams_setup.beams[i].eps[2]);
        printf("s_0 = %f\n", beams_setup.beams[i].s_0);
        printf("w = %f\n", beams_setup.beams[i].w);
    }

    printf("\nConstants\n--\n");
    printf("h = %f\n", cts.h);
    printf("e = %f\n", cts.e);
    printf("c = %f\n", cts.c);
    printf("k_B = %f\n", cts.k_B);
    printf("mu_B = %f\n", cts.mu_B);
    */

    return 0;
}

// Get atom
atom_t C_get_atom(){
    //
    // Variables
    //

    atom_t atom;
    int row_cter = 0;
    char row[STRING_BUFFER_SIZE];
    char *path = concatenate_ROOT_PATH("parameters/atom.csv");
    char *token, *rest;
    FILE *fp;

    // Open file
    fp = fopen(path, "r");

    if (fp == NULL) {
        printf("Error to access the file \"%s\"\n", path);
        exit(0);
    }

    // Skip header
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        if(row[0] == '#') row_cter++;
        else break;
    }

    // Symbol
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable symbol is invalid in the file \"%s\"\n", path);
            exit(0);
        } else {
            if(token[2] == '\n') token[2] = '\0';

            atom.symbol = (char *) malloc(strlen(token) * sizeof(char));
            strcpy(atom.symbol, token);
        };
    }

    // Name
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable name is invalid in the file \"%s\"\n", path);
            exit(0);
        }
    }

    // Atomic number
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable Z is invalid in the file \"%s\"\n", path);
            exit(0);
        } else atom.Z = atoi(token);
    }

    // Mass
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable mass is invalid in the file \"%s\"\n", path);
            exit(0);
        } else atom.mass = atoi(token);
    }

    fclose(fp);
    return atom;
}

// Get transition
transition_t C_get_transition(){
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
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        if(row[0] == '#') row_cter++;
        else break;
    }

    // Transition rate
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable gamma is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.gamma = atof(token);
    }

    // Resonant wave length
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable lambda is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.lambda = atof(token);
    }

    // Total angular momentum of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable J_gnd is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.J_gnd = atoi(token);
    }

    // Total angular momentum of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable J_exc is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.J_exc = atoi(token);
    }

    // Landè factor of the ground state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable g_gnd is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.g_gnd = atof(token);
    }

    // Landè factor of the excited state
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable g_exc is invalid in the file \"%s\"\n", path);
            exit(0);
        } else transition.g_exc = atof(token);
    }

    fclose(fp);
    return transition;
}

// Get conditions
conditions_t C_get_conditions(){
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
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        if(row[0] == '#') row_cter++;
        else break;
    }

    // Initial temperature 
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable T_0 is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.T_0 = atof(token);
    }

    // Magnetic field gradient
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable B_0 is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.B_0 = atof(token);
    }

    // Gravity
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable g_bool is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.g_bool = atoi(token);
    }

    // Maximum number of iteration
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable i_max is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.i_max = atoi(token);
    }

    // Maximum distance
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable r_max is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.r_max = atof(token);
    }

    // Number of simulations
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable num_sim is invalid in the file \"%s\"\n", path);
            exit(0);
        }
    }

    // Number of bins
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable num_bins is invalid in the file \"%s\"\n", path);
            exit(0);
        } else conditions.num_bins = atoi(token);
    }

    fclose(fp);
    return conditions;
}

// Get beams setup
beams_setup_t C_get_beams(){
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
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        if(row[0] == '#') row_cter++;
        else break;
    }

    // Get beams
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        //
        // Laser detuning
        //

        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable num_bins is invalid in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].delta = atof(token);

        //
        // Wave vector direction
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable k_dic is invalid in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].k_dic = C_get_float_array(token, &n);

        //
        // Polarization vector
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable eps is invalid in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].eps = C_get_float_array(token, &n);

        //
        // Peak of the saturation parameter
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable s_0 is invalid in the file \"%s\"\n", path);
            exit(0);
        } else beams[num_beams].s_0 = atof(token);

        //
        // Waist Radius
        //

        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable w is invalid in the file \"%s\"\n", path);
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
constants_t C_get_constants(){
    //
    // Variables
    //

    constants_t cts;
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
    while(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        if(row[0] == '#') row_cter++;
        else break;
    }

    // Planck constant
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable h is invalid in the file \"%s\"\n", path);
            exit(0);
        } else cts.h = atof(token);
    }

    // Elementary charge
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable e is invalid in the file \"%s\"\n", path);
            exit(0);
        } else cts.e = atof(token);
    }

    // Speed of light
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable c is invalid in the file \"%s\"\n", path);
            exit(0);
        } else cts.c = atof(token);
    }

    // Boltzmann constant
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable k_B is invalid in the file \"%s\"\n", path);
            exit(0);
        } else cts.k_B = atof(token);
    }

    // Bohr magneton
    if(!(fgets(row, STRING_BUFFER_SIZE, fp) == NULL)){
        rest = row;
        token = strtok_r(rest, DELIM, &rest); // Variable name
        token = strtok_r(rest, DELIM, &rest); // Value

        if(!token){
            printf("Variable mu_B is invalid in the file \"%s\"\n", path);
            exit(0);
        } else cts.mu_B = atof(token);
    }

    fclose(fp);
    return cts;
}

//
// Utility functions
//

// Get int array (max 124 numbers) from string in the format [i1 i2 ... in]
int *C_get_int_array(char *str, int *size){
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
float *C_get_float_array(char *str, int *size){
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