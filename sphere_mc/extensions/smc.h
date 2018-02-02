
struct energy_parameters {
    int natoms ;
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
} ; // end of struct energy_parameters

struct system_parameters {

    const int number_of_steps ;
    float cell_length ; 
    float delta ;
    int ncell_1d ;
    int ncell ;
    int mapsize ;
    const char * dcdfile_name ;
} ; // end of struct system_parameters


void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, energy_parameters ep, system_parameters sp)  ;

float energy(float *x_array, float *y_array, float *z_array, int *atom_id,  energy_parameters p, int atom) ;

int get_my_cell(system_parameters sp, float x, float y, float z) ;

float linked_list_energy(float *x_array, float *y_array, float *z_array, int *atom_id,  int *linked_list, int *head_of_chain_list, energy_parameters ep, system_parameters sp, int atom) ; 

float pair_energy(energy_parameters ep, float xi, float yi, float zi, float xj, float yj, float zj, int atom_i, int atom_j, int id_i, int id_j) ;
