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

float linked_list_energy(float *x_array, float *y_array, float *z_array, int *atom_id,  energy_parameters p, int atom) {

    /*
        method to calculate energy.  
        returns the sum of the hard-sphere and long-range energy
    */ 
    int i, j ;
    int i_id, j_id ;
    float overlap = 1E99 ;
    float u_long_range = 0.0 ;
    float u_long_range_1, u_long_range_2 ;

    float xi, yi, zi, xf, yf, zf, dx, dy, dz, dx2, dy2, dz2, r ;
     
    xi = x_array[atom];
    yi = y_array[atom];
    zi = z_array[atom];
   
    i_id = atom_id[atom] ;
 
    for (j = 0 ; j < p.natoms ; ++j)
    {
        if(j != atom){
            xf = x_array[j];
            yf = y_array[j];
            zf = z_array[j];
  
            dx = xf - xi ; dy = yf - yi ; dz = zf - zi ;  
            dx2 = dx * dx ; dy2 = dy * dy ; dz2 = dz * dz ;
            r = sqrt(dx2 + dy2 + dz2) ;

            j_id = atom_id[j] ;
 
            if(i_id == j_id) {
                if(i_id == 0){
                // 11 
                if (r < ( p.sigma_11 + p.sigma_11)) {
                    return overlap ;
                }
                u_long_range_1 = -p.epsilon_a_11 * pow(( p.sigma_11 / p.r_a_11 ),2.0) * exp(-(r/p.r_a_11)) ;
                u_long_range_2 =  p.epsilon_r_11 * pow(( p.sigma_11 / p.r_r_11 ),2.0) * exp(-(r/p.r_r_11)) ;


                }else {
                // 22
                if (r < ( p.sigma_22 + p.sigma_22)) {
                    return overlap ;
                }
                u_long_range_1 = -p.epsilon_a_22 * pow(( p.sigma_22 / p.r_a_22 ),2.0) * exp(-(r/p.r_a_22)) ;
                u_long_range_2 =  p.epsilon_r_22 * pow(( p.sigma_22 / p.r_r_22 ),2.0) * exp(-(r/p.r_r_22)) ;

                }

            }else {
            //12
            if (r < ( p.sigma_11 + p.sigma_22)) {
                    return overlap ;
            }
            u_long_range_1 = -p.epsilon_a_12 * pow(( p.sigma_12 / p.r_a_12 ),2.0) * exp(-(r/p.r_a_12)) ;
            u_long_range_2 =  p.epsilon_r_12 * pow(( p.sigma_12 / p.r_r_12 ),2.0) * exp(-(r/p.r_r_12)) ;
            
            }
           
            u_long_range += u_long_range_1 + u_long_range_2 ; 
        }
    } // end of j-loop


    return u_long_range ;

} ; //end of energy

