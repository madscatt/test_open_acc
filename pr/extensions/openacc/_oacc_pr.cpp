#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "oacc_pr.h"

int get_distances(double ***coor, const int nframes, const int natoms) {

    int i,j,k,d,pair,npairs ;
    double x1, y1, z1, x2, y2, z2, sdist ;
    double dx2, dy2, dz2 ;
    double dist[(natoms * (natoms -1))/2] ;
    printf("oacc: %d\n", nframes) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: coor[0][0][0] = %f\n", coor[0][0][0]) ;
    npairs = (natoms * (natoms - 1))/2 ;
   
    for(d=0 ; d < npairs ; d++) {
        dist[d] = 0.0 ;
    }
    //#pragma acc parallel loop
    #pragma acc enter data create(dist[npairs],coor[nframes][natoms][3])
    for(i=0 ; i < nframes ; i++){
        pair = 0 ;
        for(j=0 ; j < natoms-1 ; j++){
            x1 = coor[i][j][0] ;
            y1 = coor[i][j][1] ;
            z1 = coor[i][j][2] ;
            #pragma acc update self(dist[npairs])
            for(k=j+1 ; k < natoms ; k++){
                x2 = coor[i][k][0] ;
                y2 = coor[i][k][1] ;
                z2 = coor[i][k][2] ;
                dx2 = (x1 - x2) * (x1 - x2) ;
                dy2 = (y1 - y2) * (y1 - y2) ;
                dz2 = (z1 - z2) * (z1 - z2) ;
                sdist = sqrt(dx2 + dy2 + dz2) ;
                dist[pair] += sdist ; 
                pair++ ;
            }
        }
    }

    printf("dist[0] = %f\n", dist[0]/double(nframes)) ;
    printf("dist[npairs-1] = %f\n", dist[npairs-1]/double(nframes)) ;

    //delete [] dist ;

    return 0 ;
}
