
struct energy_parameters {
    int  natoms ;
    float temperature ;
    float sigma_11  ;
    float sigma_22 ;
    float sigma_12  ;
    float epsilon_a_11 ;
    float epsilon_a_22 ;
    float epsilon_a_12 ;
    float epsilon_r_11 ;
    float epsilon_r_22 ;
    float epsilon_r_12 ;
    float r_a_11 ;
    float r_a_22 ;
    float r_a_12 ;
    float r_r_11 ;
    float r_r_22 ;
    float r_r_12 ;
    float energy ;
    float beta ;
} ;


void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, const char *dcdfile_name, int number_of_steps, int ncell, int ncell_1d, float cell_length, float delta, energy_parameters parameters)  ;

