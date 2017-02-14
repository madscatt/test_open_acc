#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "oacc_pr.h"

int get_distances(double ***coor, const int nframes, const int natoms) {

    int i,j,k ;
    unsigned long long d, npairs ;
    //int d, npairs ;
    double x1, y1, z1, x2, y2, z2, sdist ;
    double dx2, dy2, dz2 ;
    //double *restrict dist[(natoms * (natoms -1))/2] ;
    double dist[(natoms * (natoms -1))/2] ;
    //int count[(natoms * (natoms -1))/2] ;
    //int count_pairs ;
    unsigned long long local_count ; 
    //int local_count ; 

    printf("oacc: %d\n", nframes) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: coor[0][0][0] = %f\n", coor[0][0][0]) ;
    npairs = (natoms * (natoms - 1))/2 ;

    //count_pairs = 0 ; 

    for(d=0 ; d < npairs ; d++) {
        dist[d] = 0.0 ;
       // count[d] = 0 ;
    }
    #pragma acc data copyin(coor[nframes][natoms][3]) copy(dist[npairs])
    {
    for(i=0 ; i < nframes ; i++){
        #pragma acc parallel loop
        {
        for(j=0 ; j < natoms-1 ; j++){
            x1 = coor[i][j][0] ;
            y1 = coor[i][j][1] ;
            z1 = coor[i][j][2] ;
            #pragma acc loop 
            {
            for(k=j+1 ; k < natoms ; k++){
                x2 = coor[i][k][0] ;
                y2 = coor[i][k][1] ;
                z2 = coor[i][k][2] ;
                dx2 = (x1 - x2) * (x1 - x2) ;
                dy2 = (y1 - y2) * (y1 - y2) ;
                dz2 = (z1 - z2) * (z1 - z2) ;
                sdist = sqrt(dx2 + dy2 + dz2) ;
                local_count = ((j*natoms)-((j*(j+1))/2))+k-(j+1) ;
                dist[local_count] += sdist ; 
        //        count[local_count] += 1 ; 
                //if(i==0 && j<2){ 
                //if(i==0){
                 //   printf("%i\t%i\t%i\t%i\ttdist[c] = %f\n", i,j,k,local_count,dist[local_count]) ;
    //           }
               // count_pairs++;
            }
            }
        //if(i==0) printf("\n") ;
        }
        }
    }
    }

    //for(i=0 ; i < 30 ; i++){
     //   printf("dist[%i] = %f\n", i,dist[i]/double(nframes)) ;
     ////   printf("dist[%i] = %f\n", i,dist[i]/double(count[i]+1)) ;
    //} 

    //for(i=0 ; i < npairs ; i++){
     //   printf("count[%i] = %i\n", i,count[i]) ;
      //  //printf("count[%i] = %f\n", i,count[i]/double(nframes)) ;
   // }

    //printf("count_pairs = %i\n", count_pairs) ;
    //printf("npairs * nframes = %llu\n", npairs*nframes) ;
    printf("npairs * nframes = %i\n", npairs*nframes) ;

    //delete [] dist ;

    return 0 ;
}
