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

int get_icell(system_parameters sp, int ix, int iy, int iz) {
    
    int xpart, ypart, zpart ;

    xpart = ((ix + sp.ncell_1d) % sp.ncell_1d) ; 
    ypart = ((iy + sp.ncell_1d) % sp.ncell_1d) * sp.ncell_1d ; 
    zpart = ((iz + sp.ncell_1d) % sp.ncell_1d) * sp.ncell_1d * sp.ncell_1d ;

    return xpart + ypart + zpart ; 

} ; // end of get_icell 

void  make_neighbor_map(int *map, system_parameters sp){

    int i, ix, iy, iz, icell, imap ;

    for(i = 0 ; i < sp.mapsize ; ++i){
        map[i] = 0 ; 
    } // end of i loop over sp.mapsize

    for(iz = 0; iz < sp.ncell_1d ; ++iz){

        for(iy = 0; iy < sp.ncell_1d ; ++iy){

            for(ix = 0; ix < sp.ncell_1d ; ++ix){

                icell = get_icell(sp, ix, iy, iz) ;
                imap = icell * 26 ;

                map[imap] = get_icell(sp, ix, iy, iz + 1) ;                 // 1
                map[imap + 1] = get_icell(sp, ix, iy, iz - 1) ;             // 2

                map[imap + 2] = get_icell(sp, ix, iy + 1, iz) ;             // 3
                map[imap + 3] = get_icell(sp, ix, iy - 1, iz) ;             // 4

                map[imap + 4] = get_icell(sp, ix + 1, iy, iz) ;             // 5
                map[imap + 5] = get_icell(sp, ix - 1, iy, iz) ;             // 6

                map[imap + 6] = get_icell(sp, ix + 1, iy + 1, iz) ;         // 7
                map[imap + 7] = get_icell(sp, ix + 1, iy - 1, iz) ;         // 8
                map[imap + 8] = get_icell(sp, ix - 1, iy + 1, iz) ;         // 9
                map[imap + 9] = get_icell(sp, ix - 1, iy - 1, iz) ;         // 10

                map[imap + 10] = get_icell(sp, ix + 1, iy, iz + 1) ;        // 11
                map[imap + 11] = get_icell(sp, ix + 1, iy, iz - 1) ;        // 12
                map[imap + 12] = get_icell(sp, ix - 1, iy, iz + 1) ;        // 13
                map[imap + 13] = get_icell(sp, ix - 1, iy, iz - 1) ;        // 14

                map[imap + 14] = get_icell(sp, ix + 1, iy + 1, iz + 1) ;    // 15
                map[imap + 15] = get_icell(sp, ix + 1, iy + 1, iz - 1) ;    // 16
                map[imap + 16] = get_icell(sp, ix + 1, iy - 1, iz + 1) ;    // 17
                map[imap + 17] = get_icell(sp, ix + 1, iy - 1, iz - 1) ;    // 18
                map[imap + 18] = get_icell(sp, ix - 1, iy + 1, iz + 1) ;    // 19
                map[imap + 19] = get_icell(sp, ix - 1, iy + 1, iz - 1) ;    // 20
                map[imap + 20] = get_icell(sp, ix - 1, iy - 1, iz + 1) ;    // 21
                map[imap + 21] = get_icell(sp, ix - 1, iy - 1, iz - 1) ;    // 22

                map[imap + 22] = get_icell(sp, ix, iy + 1, iz + 1) ;        // 23
                map[imap + 23] = get_icell(sp, ix, iy + 1, iz - 1) ;        // 24
                map[imap + 24] = get_icell(sp, ix, iy - 1, iz + 1) ;        // 25
                map[imap + 25] = get_icell(sp, ix, iy - 1, iz - 1) ;        // 26


                if (ix == 0 and iy == 0 and iz == 0){
                    std::cout <<  "map0 = " << map[4] << std::endl ; 
                    std::cout <<  "map0 = " << map[6] << std::endl ; 
                    std::cout <<  "map0 = " << map[2] << std::endl ; 
                    std::cout <<  "map0 = " << map[8] << std::endl ; 
                } // end of if ix, iy, iz debugging
            } // end ix loop over sp.ncell_1d
        } // end of iy loop over sp.ncell_1d
    } // end of iz loop over sp.ncell_1d

    return ; 

} ; // end of make_neighbor_map


int get_my_cell(system_parameters sp, float x, float y, float z) {

        int icell_x, icell_y, icell_z, icel ;

        icell_x = int((x + (0.5 * sp.box_length))/sp.cell_length) ; 
        icell_y = int((y + (0.5 * sp.box_length))/sp.cell_length) ;
        icell_z = int((z + (0.5 * sp.box_length))/sp.cell_length) ;

        icel = icell_x + (sp.ncell_1d * icell_y) + (sp.ncell_1d * sp.ncell_1d * icell_z) ;

        return icel ;

} ; // end of get_my_cell


void update_cell_list(float *x_array, float *y_array, float *z_array, int *linked_list, int *head_of_chain_list, energy_parameters ep, system_parameters sp) {

    int i, icel ; 

    for (i = 0; i < ep.natoms ; ++i) {
        linked_list[i] = -1 ;

    } // end of i loop over ep.natoms

    for (i = 0; i < sp.ncell ; ++i) {
        head_of_chain_list[i] = -1 ;

    } // end of loop over sp.ncell

    for (i = 0 ; i < ep.natoms ; ++i){

        icel = get_my_cell(sp, x_array[i], y_array[i], z_array[i]) ;
        linked_list[i] = head_of_chain_list[icel] ;
        head_of_chain_list[icel] =  i ;

    } // end of i loop over ep.natoms

    return ;

} // end of update_cell_list


int surface_move(float *x_array, float * y_array, float *z_array, int *atom_id, int *linked_list, int *head_of_chain_list, int *map, energy_parameters ep, system_parameters sp, int atom) {

    float x, y, z, r, dx, dy, dz, tx, ty, tz ;
    float norm, final_energy, initial_energy, boltz, delta_energy, ran ; 

    std::random_device rd ;
    std::random_device rd2 ;
    std::mt19937 mt(rd()) ;
    std::mt19937 mt2(rd2()) ;
    std::uniform_real_distribution<float> dist(-1.0,1.0) ;
    std::uniform_real_distribution<float> dist2(0.0,1.0) ;

    bool accepted = false ;

    x = x_array[atom] ;
    y = y_array[atom] ;
    z = z_array[atom] ;

    tx = x ; ty = y ; tz = z ;
    
    r = sqrt(x*x+y*y+z*z) ;

    dx = ep.max_displacement * dist(mt);
    dy = ep.max_displacement * dist(mt);
    dz = ep.max_displacement * dist(mt);

    x += dx ; y += dy ; z += dz ;
    norm = sqrt(x*x + y*y + z*z) ;
    x *= r/norm ;
    y *= r/norm ;
    z *= r/norm ;

    // get initial energy

    initial_energy = linked_list_energy(x_array, y_array, z_array, atom_id, linked_list, head_of_chain_list, ep, sp, atom) ;

    x_array[atom] = x ;
    y_array[atom] = y ;
    z_array[atom] = z ;

    // get final energy
    final_energy = linked_list_energy(x_array, y_array, z_array, atom_id, linked_list, head_of_chain_list, ep, sp, atom) ;

    if (final_energy < initial_energy) {
        accepted = true ;
    } else {
        delta_energy = final_energy - initial_energy ;
        boltz = exp(-ep.beta * delta_energy) ;
        ran = dist2(mt2) ;
        if (ran < boltz) {
            accepted = true ;
        } // end if ran 
    } // end of if/else
    if (accepted) {
        ep.energy = final_energy ;
        return 1 ;
    } else {
        x_array[atom] = tx ;
        y_array[atom] = ty ;
        z_array[atom] = tz ;
        return 0 ;
    } // end of if/else

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

 
void smc_core(float *x_array, float *y_array, float *z_array, int *atom_id, energy_parameters ep, system_parameters sp) {

    int i, j, atom ;
    int step_result ;

    std::random_device rd ;
    std::mt19937 mt(rd()) ;
    std::uniform_int_distribution<int> random_atom(0,ep.natoms - 1) ;

    int linked_list[ep.natoms] ;
    int head_of_chain_list[sp.ncell] ;
    int number_accepted = 0 ;

    FILE * filepointer ;
    filepointer = open_dcd(sp.dcdfile_name, ep.natoms) ;

    sp.mapsize = int(pow(sp.ncell_1d,3)) * 26 ;
    int map[sp.mapsize] ;

    make_neighbor_map(map, sp) ;

    for(i = 0 ; i < sp.number_of_steps ; ++i) {
        std::cout << i << " " << std::flush ;   

        for(j = 0 ; j < ep.natoms ; ++j){ 

            atom  = random_atom(mt) ;

            update_cell_list(x_array, y_array, z_array, linked_list, head_of_chain_list, ep, sp) ;

            number_accepted += surface_move(x_array, y_array, z_array, atom_id, linked_list, head_of_chain_list, map, ep, sp, atom) ;
        } // end of j loop over ep.natoms

        step_result = write_dcdstep(filepointer, ep.natoms, x_array, y_array, z_array, i) ;
        
       // if(step_result != 1){
        //    std::cout << "error writing dcd frame: " << sp.dcdfile_name << "step = " << i << std::endl ;
       // } // end if step_result

    } // end of i loop over sp.number_of_steps

    step_result = close_dcd_write(filepointer) ;

    //if(step_result != 1){
     //   std::cout << "error closing dcd file: " << sp.dcdfile_name << std::endl ;
    //} // end if step_result

    return ;

} ; // end of smc_core


