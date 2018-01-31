#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"
#include "oacc_smc.h"
#include <vector> 

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <random>

#include <stdio.h>
#include "dcdio.h"

#include "smc.h"

using namespace std;

extern void get_distances();

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

float energy(float *x_array, float *y_array, float *z_array, int *atom_id,  energy_parameters p, int atom) {

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
 
    for (j = 0 ; j < p.natoms ; j++)
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
}

void update_cell_list(int *linked_list, int *head_of_chain_list, float *x_array, float *y_array, float *z_array, float delta, float cell_length, int natoms, int ncell_1d, int ncell) {

    int i ; 

    int icell_x, icell_y, icell_z, icel ;

    for (i = 0; i < natoms ; i++) {
        linked_list[i] = 0 ;
    }

    for (i = 0; i < ncell ; i++) {
        head_of_chain_list[i] = -1 ;
    }

    for (i = 0 ; i < natoms ; i++){

        icell_x = int((x_array[0] + (0.5 * delta))/cell_length) ; 
        icell_y = int((y_array[1] + (0.5 * delta))/cell_length) ;
        icell_z = int((z_array[2] + (0.5 * delta))/cell_length) ;

        icel = icell_x + (ncell_1d * icell_y) + (ncell_1d * ncell_1d * icell_z) ;

        linked_list[i] = head_of_chain_list[icel] ;
        head_of_chain_list[icel] =  i ;

    }

    return ;

} // end of update_cell_list


int surface_move(float *x_array, float * y_array, float *z_array, int *atom_id, int i, energy_parameters parameters, int *linked_list, int *head_of_chain_list, int *nmap, int ncell) {

    float x, y, z, r, dx, dy, dz, tx, ty, tz ;
    float max_disp, norm ;
    float u_long_range, boltz, delta_energy, ran ; 
    float initial_energy ;

    std::random_device rd ;
    std::random_device rd2 ;
    std::mt19937 mt(rd()) ;
    std::mt19937 mt2(rd2()) ;
    std::uniform_real_distribution<float> dist(-1.0,1.0) ;
    std::uniform_real_distribution<float> dist2(0.0,1.0) ;

    bool accepted = false ;

    x = x_array[i] ; 
    y = y_array[i] ; 
    z = z_array[i] ; 

    tx = x ; ty = y ; tz = z ;
    //std::cout << "c : x = " << x << std::endl ;
    
    r = sqrt(x*x+y*y+z*z) ;

    max_disp = 0.05 ;
    //max_disp = 1.0 ;
    
    dx = max_disp * dist(mt);
    dy = max_disp * dist(mt);
    dz = max_disp * dist(mt);

    x += dx ; y += dy ; z += dz ;
    norm = sqrt(x*x + y*y + z*z) ;
    x *= r/norm ;
    y *= r/norm ;
    z *= r/norm ;

    x_array[i] = x ;
    y_array[i] = y ;
    z_array[i] = z ;

    u_long_range = energy(x_array, y_array, z_array, atom_id, parameters, i) ;

    if (u_long_range < parameters.energy) {
        accepted = true ;
    } else {
        delta_energy = u_long_range - parameters.energy ;
        boltz = exp(-parameters.beta * delta_energy) ;
        ran = dist2(mt2) ;
        if (ran < boltz) {
            accepted = true ;
        } 
    }
    if (accepted) {
        parameters.energy = u_long_range ;
        return 1 ;
    } else {
        x_array[i] = tx ;
        y_array[i] = ty ;
        z_array[i] = tz ;
        return 0 ;
    }
}


void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, int *nmap, char *filename, int number_of_steps, int ncell, int ncell_1d, float cell_lenght, float delta, energy_parameters parameters) {

    int i, j ;

    int linked_list[parameters.natoms] ;
    int head_of_chain_list[ncell] ;
    int number_accepted = 0 ;

    // Open DCD file for writing 

    FILE * filepointer  ;
    filepointer = open_dcd_write(filename) ;

    int istart = 0 ;
    int nsavc = 1 ;
    int nset = 1 ;
    double dcddelta = 1.0 ;
    int headerresult ;
    int stepresult ;
    headerresult = write_dcdheader(filepointer, filename, parameters.natoms, nset, istart, nsavc, dcddelta) ;

    std::cout << "write dcd header result = " << headerresult << std::endl ;

    for(i = 0 ; i < number_of_steps ; i++) {
     //   std::cout << "c : x_array[0] = " << x_array[0] << std::endl ;
        std::cout << i << " " << std::flush ;   

        for(j = 0 ; j < parameters.natoms ; j++){ 
            update_cell_list(linked_list, head_of_chain_list, x_array, y_array, z_array, delta, cell_length, parameters.natoms, ncell_1d, ncell) ;

            number_accepted += surface_move(x_array, y_array, z_array, atom_id, j, parameters, linked_list, head_of_chain_list, nmap, ncell) ;
        }

        stepresult = write_dcdstep(filepointer, parameters.natoms, x_array, y_array, z_array, i) ;

    } // end of i loop

    stepresult = close_dcd_write(filepointer) ;

    return ;

}


