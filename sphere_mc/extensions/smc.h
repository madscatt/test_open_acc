
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
    float r_cutoff ;
    float max_displacement ;
    float contrast_1 ;
    float contrast_2 ;
} ;

struct system_parameters {

    const int number_of_steps ;
    float cell_length ; 
    float delta ;
    int ncell_1d ;
    int ncell ;
    const char * dcdfile_name ;
} ;




void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, energy_parameters parameters, system_parameters system_parameters)  ;

float energy(float *x_array, float *y_array, float *z_array, int *atom_id,  energy_parameters p, int atom) ;

int get_my_cell(float x, float y, float z, system_parameters system_parameters) ;

float linked_list_energy(float *x_array, float *y_array, float *z_array, int *atom_id,  int *linked_list, int *head_of_chain_list, energy_parameters parameters, system_parameters system_parameters, int atom) ; 

float pair_energy(float xi, float yi, float zi, float xj, float yj, float zj, int atom_i, int atom_j, int id_i, int id_j, energy_parameters p) ;
