#include <iostream>
#include <math.h>
#include <random>

#ifndef SMC_H
#define SMC_H

#include "smc.h"

#endif

#ifndef DCDIO_H
#define DCDIO_H

#include "dcdio.h"

#endif


/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

int get_icell(int ix, int iy, int iz, int ncell_1d) {
    
    int xpart, ypart, zpart ;

    xpart = ((ix + ncell_1d) % ncell_1d) ; 
    ypart = ((iy + ncell_1d) % ncell_1d) * ncell_1d ; 
    zpart = ((iz + ncell_1d) % ncell_1d) * ncell_1d * ncell_1d ;

    return xpart + ypart + zpart ; 

} ; // end of get_icell 

void  make_neighbor_map(int *map, int ncell_1d){

    int ix, iy, iz, icell, imap ;

    for(iz = 0; iz < ncell_1d ; ++iz){

        for(iy = 0; iy < ncell_1d ; ++iy){

            for(ix = 0; ix < ncell_1d ; ++ix){

                icell = get_icell(ix, iy, iz, ncell_1d) ;
                imap = icell * 26 ;

                map[imap] = get_icell(ix, iy, iz + 1, ncell_1d) ;                       // 1
                map[imap + 1] = get_icell(ix, iy, iz - 1, ncell_1d) ;                       // 2

                map[imap + 2] = get_icell(ix, iy + 1, iz, ncell_1d) ;                       // 3
                map[imap + 3] = get_icell(ix, iy - 1, iz, ncell_1d) ;                       // 4

                map[imap + 4] = get_icell(ix + 1, iy, iz, ncell_1d) ;                       // 5
                map[imap + 5] = get_icell(ix - 1, iy, iz, ncell_1d) ;                       // 6

                map[imap + 6] = get_icell(ix + 1, iy + 1, iz, ncell_1d) ;                   // 7
                map[imap + 7] = get_icell(ix + 1, iy - 1, iz, ncell_1d) ;                   // 8
                map[imap + 8] = get_icell(ix - 1, iy + 1, iz, ncell_1d) ;                   // 9
                map[imap + 9] = get_icell(ix - 1, iy - 1, iz, ncell_1d) ;                   // 10

                map[imap + 10] = get_icell(ix + 1, iy, iz + 1, ncell_1d) ;                   // 11
                map[imap + 11] = get_icell(ix + 1, iy, iz - 1, ncell_1d) ;                   // 12
                map[imap + 12] = get_icell(ix - 1, iy, iz + 1, ncell_1d) ;                   // 13
                map[imap + 13] = get_icell(ix - 1, iy, iz - 1, ncell_1d) ;                   // 14

                map[imap + 14] = get_icell(ix + 1, iy + 1, iz + 1, ncell_1d) ;               // 15
                map[imap + 15] = get_icell(ix + 1, iy + 1, iz - 1, ncell_1d) ;               // 16
                map[imap + 16] = get_icell(ix + 1, iy - 1, iz + 1, ncell_1d) ;               // 17
                map[imap + 17] = get_icell(ix + 1, iy - 1, iz - 1, ncell_1d) ;               // 18
                map[imap + 18] = get_icell(ix - 1, iy + 1, iz + 1, ncell_1d) ;               // 19
                map[imap + 19] = get_icell(ix - 1, iy + 1, iz - 1, ncell_1d) ;               // 20
                map[imap + 20] = get_icell(ix - 1, iy - 1, iz + 1, ncell_1d) ;               // 21
                map[imap + 21] = get_icell(ix - 1, iy - 1, iz - 1, ncell_1d) ;               // 22

                map[imap + 22] = get_icell(ix, iy + 1, iz + 1, ncell_1d) ;                   // 23
                map[imap + 23] = get_icell(ix, iy + 1, iz - 1, ncell_1d) ;                   // 24
                map[imap + 24] = get_icell(ix, iy - 1, iz + 1, ncell_1d) ;                   // 25
                map[imap + 25] = get_icell(ix, iy - 1, iz - 1, ncell_1d) ;                   // 26


                if (ix == 0 and iy == 0 and iz == 0){
                    std::cout <<  "map0 = " << map[4] << std::endl ; 
                    std::cout <<  "map0 = " << map[6] << std::endl ; 
                    std::cout <<  "map0 = " << map[2] << std::endl ; 
                    std::cout <<  "map0 = " << map[8] << std::endl ; 
                } // end of if debugging
            } // end loop over ix
        } // end of loop over iy
    } // end of loop over iz

    return ; 
} ; // end of make_neighbor_map


int get_my_cell(float x, float y, float z, float delta, float cell_length, int ncell_1d) {

        int icell_x, icell_y, icell_z, icel ;

        icell_x = int((x + (0.5 * delta))/cell_length) ; 
        icell_y = int((y + (0.5 * delta))/cell_length) ;
        icell_z = int((z + (0.5 * delta))/cell_length) ;

        icel = icell_x + (ncell_1d * icell_y) + (ncell_1d * ncell_1d * icell_z) ;

        return icel ;

} ; // end of get_my_cell


void update_cell_list(int *linked_list, int *head_of_chain_list, float *x_array, float *y_array, float *z_array, float delta, float cell_length, int natoms, int ncell_1d, int ncell) {

    int i, icel ; 

    for (i = 0; i < natoms ; ++i) {
        linked_list[i] = -1 ;

    } // end of loop over natoms

    for (i = 0; i < ncell ; ++i) {
        head_of_chain_list[i] = -1 ;

    } // end of loop over ncell

    for (i = 0 ; i < natoms ; ++i){

        icel = get_my_cell(x_array[i], y_array[i], z_array[i], delta, cell_length, ncell_1d) ;
        linked_list[i] = head_of_chain_list[icel] ;
        head_of_chain_list[icel] =  i ;

    } // end of loop over natoms

    return ;

} // end of update_cell_list


int surface_move(float *x_array, float * y_array, float *z_array, int *atom_id, int i, energy_parameters parameters, int *linked_list, int *head_of_chain_list, int *map, int ncell) {

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

} ; // end of surface move

FILE *open_dcd( const char *dcdfile_name, int natoms) {

    char * filename = (char*)dcdfile_name ;

    FILE * filepointer  ;
    filepointer = open_dcd_write(filename) ;
    int istart = 0 ;
    int nsavc = 1 ;
    int nset = 1 ;
    double dcddelta = 1.0 ;
    int headerresult ;
    headerresult = write_dcdheader(filepointer, filename, natoms, nset, istart, nsavc, dcddelta) ;

    std::cout << "write dcd header result = " << headerresult << std::endl ;

    return filepointer  ;

} ; // end of open_dcd

 
void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, const char *dcdfile_name, int number_of_steps, int ncell, int ncell_1d, float cell_length, float delta, energy_parameters parameters) {

    int i, j ;
    int mapsize, stepresult ;

    int linked_list[parameters.natoms] ;
    int head_of_chain_list[ncell] ;
    int number_accepted = 0 ;

    FILE * filepointer ;
    filepointer = open_dcd(dcdfile_name, parameters.natoms) ;

    mapsize = int(pow(ncell_1d,3)) * 26 ;
    int map[mapsize] ;

    for(i = 0 ; i < mapsize ; ++i){
        map[i] = 0 ; 
    } // end of i loop

    make_neighbor_map(map, ncell_1d) ;

    for(i = 0 ; i < number_of_steps ; ++i) {
        std::cout << i << " " << std::flush ;   

        // need to pick a random atom, not in succession

        for(j = 0 ; j < parameters.natoms ; ++j){ 
            update_cell_list(linked_list, head_of_chain_list, x_array, y_array, z_array, delta, cell_length, parameters.natoms, ncell_1d, ncell) ;

            number_accepted += surface_move(x_array, y_array, z_array, atom_id, j, parameters, linked_list, head_of_chain_list, map, ncell) ;
        }

        stepresult = write_dcdstep(filepointer, parameters.natoms, x_array, y_array, z_array, i) ;

    } // end of i loop

    stepresult = close_dcd_write(filepointer) ;

    return ;

} ; // end of smc_core


