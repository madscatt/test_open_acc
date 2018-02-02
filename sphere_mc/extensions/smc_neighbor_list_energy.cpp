#include <iostream>
#include <math.h>

#ifndef SMC_H
#define SMC_H

#include "smc.h"

#endif

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

float linked_list_energy(float *x_array, float *y_array, float *z_array, int *atom_id,  int *linked_list, int *head_of_chain_list, energy_parameters ep, system_parameters sp, int atom){ 

    int i, j, icel ;
    float xi, yi, zi, xj, yj, zj ;

    float energy = 0.0 ;
    float overlap = 1E99 ;

    xi = x_array[atom] ; yi = y_array[atom] ; zi = z_array[atom] ;

    // find the cell that atom is in

    icel = get_my_cell(sp, xi, yi, zi) ;

    // get the index of the head of chain for the atom in that cell
   
    i = head_of_chain_list[icel] ;
   
    // return if head of chain is exhausted ... not sure this should happen 
    if(i < 0){
        std::cout << "for atom = " << atom << "hoc = " << i << "not sure this should happen" << std::endl ;
        std::cout << "x[" << atom << "] = " << xi << std::endl ;
        std::cout << "y[" << atom << "] = " << yi << std::endl ;
        std::cout << "z[" << atom << "] = " << zi << std::endl ;
        std::cout << "for atom = " << atom << "energy = " << energy << "not sure this should happen" << std::endl ;
        return overlap ; // a hack
        
    } ; // end of if i < 0

    // grab the atom in the link list for atoms in that cell
    j = linked_list[i] ;

    // return if linked list is exhausted 
    while(j < 0){

        xj = x_array[j] ; yj = y_array[j] ; zj = z_array[j] ;
   
        energy += pair_energy(ep, xi, yi, zi, xj, yj, zj, atom, j, atom_id[atom], atom_id[j]) ;
    
        j = linked_list[j] ;

    } ; // end of while j < 0

    return overlap ;

} ; // end of linked_list_energy 

float pair_energy(energy_parameters ep, float xi, float yi, float zi, float xj, float yj, float zj, int atom_i, int atom_j, int i_id, int j_id) { 

    /*
        method to calculate energy.  
        returns the sum of the hard-sphere and long-range energy
    */ 
    float overlap = 1E99 ;
    float u_long_range = 0.0 ;
    float u_long_range_1 = 0.0 ;
    float u_long_range_2 = 0.0 ;

    float dx, dy, dz, dx2, dy2, dz2, r ;
     
    if(atom_j != atom_i){
  
        dx = xj- xi ; dy = yj - yi ; dz = zj - zi ;  
        dx2 = dx * dx ; dy2 = dy * dy ; dz2 = dz * dz ;
        r = sqrt(dx2 + dy2 + dz2) ;

        if (r < ep.r_cutoff){

            if(i_id == j_id) {
                if(i_id == 0){
                // 11 
                    if (r < ( ep.sigma_11 + ep.sigma_11)) {
                        return overlap ;
                    } // end of overlap check

                    u_long_range_1 = -ep.epsilon_a_11 * pow(( ep.sigma_11 / ep.r_a_11 ),2.0) * exp(-(r/ep.r_a_11)) ;
                    u_long_range_2 =  ep.epsilon_r_11 * pow(( ep.sigma_11 / ep.r_r_11 ),2.0) * exp(-(r/ep.r_r_11)) ;

                }else {
                // 22
                    if (r < ( ep.sigma_22 + ep.sigma_22)) {
                        return overlap ;
                    } // end of overlap check

                    u_long_range_1 = -ep.epsilon_a_22 * pow(( ep.sigma_22 / ep.r_a_22 ),2.0) * exp(-(r/ep.r_a_22)) ;
                    u_long_range_2 =  ep.epsilon_r_22 * pow(( ep.sigma_22 / ep.r_r_22 ),2.0) * exp(-(r/ep.r_r_22)) ;

                } // end of if/else

            }else {
            //12
                if (r < ( ep.sigma_11 + ep.sigma_22)) {
                    return overlap ;
                } // end of overlap check

                u_long_range_1 = -ep.epsilon_a_12 * pow(( ep.sigma_12 / ep.r_a_12 ),2.0) * exp(-(r/ep.r_a_12)) ;
                u_long_range_2 =  ep.epsilon_r_12 * pow(( ep.sigma_12 / ep.r_r_12 ),2.0) * exp(-(r/ep.r_r_12)) ;
            
            } // end of if/else (starting with i_id == 0)
           
            u_long_range += u_long_range_1 + u_long_range_2 ; 
        } // end if (r < ep.r_cutoff) 
    } // end of if(atom_j != atom_i)

    return u_long_range ;

} ; // end of pair_energy

