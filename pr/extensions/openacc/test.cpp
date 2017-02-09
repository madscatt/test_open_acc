
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>

#include "oacc_pr.h"

//int get_distances(double ***coor, const int nframes, const int natoms) {

int main(){

    const int nframes=10;
    const int natoms=10;
    //double coor[2][3][4] = { { {1., 2., 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4} },
//                     { {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4} } };

    double*** coor ;
    coor = new double**[nframes] ;

    for(int x = 0 ; x < nframes ; x++){
        coor[x] = new double*[natoms] ;
        for(int y = 0 ; y < natoms ; y++){
            coor[x][y] = new double[3] ;
            for(int z = 0 ; z < 3 ; z++){
                coor[x][y][z] = 0.0 ;
            }
        }
    }
    coor[0][0][0] = 1.002 ;

    printf("test: nframes = %d\n", nframes) ;   
    printf("test: coor[0][0][0] = %f\n", coor[0][0][0]) ;   
    printf("calling get_distances\n");
    get_distances(coor, nframes, natoms) ;
    printf("done with get_distances\n");
 

    return 0 ;
}
