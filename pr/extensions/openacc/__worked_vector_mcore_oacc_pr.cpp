#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "oacc_pr.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <typeinfo>

using namespace std;

void get_distances(double ***coor, const int nframes, const int natoms, std::vector<double>& dist) {

    int i,j,k ;
    unsigned long long npairs ;
    double x1, y1, z1, x2, y2, z2, sdist ;
    double dx2, dy2, dz2 ;
    unsigned long long local_count ; 
    unsigned long long z ; 

    /// temporary printing for debugging
    //
    //
    std::ostringstream sstream;
    std::string remark = "#hello";
    const std::string filename = "dum_oacc.txt";
    
    std::ofstream outfile(filename.c_str()) ;
    outfile << remark << std::endl;

    //
    
    printf("oacc: %d\n", nframes) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: coor[0][0][0] = %f\n", coor[0][0][0]) ;
    npairs = (natoms * (natoms - 1))/2 ;
    unsigned long long lc ;

    #pragma acc data copyin(coor[nframes][natoms][3]) copy(dist[npairs])
    {
    for(i=0 ; i < nframes ; i++){
        #pragma acc parallel loop
        lc=0;
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
                lc++;
            } // end of k-loop
            } // end of pragma acc loop
        } // end of j-loop
        } // end of pragma acc parallel loop
    } // end of i-loop
    } // end of pragma acc data 

    int nbins = 21 ;
    double bin_width = 1.0 ;
    double this_low_bin, this_high_bin ;
    std::vector<int> hist(nbins,0);

    for(z=0 ; z < npairs ; z++){
        for(i=0 ; i < nbins ; i++){
            this_low_bin = double(i)*bin_width ;
            this_high_bin = this_low_bin + bin_width ;
           // if(z==0){
            //    std::cout << "tlb = " << this_low_bin << "\tthb = " << this_high_bin << std::endl ;
           // }
            if(dist[z]/double(nframes) > this_low_bin && dist[z]/double(nframes) <= this_high_bin){
                hist[i] += 1 ;
                break ;
            }
        }
    }

    for(i=0 ; i < nbins ; i++){
        std::cout << hist[i] << std::endl ;
    }
/*
    std::cout << typeid(coor[0][0][0]).name() << '\n';
    std::cout << typeid(x1).name() << '\n';
    std::cout << typeid(sdist).name() << '\n';
    std::cout << typeid(dist[0]).name() << '\n';

    std::cout << "local count = " << local_count << std::endl ; 
    std::cout << "local count + 1 = " << local_count + 1 << std::endl ; 
    std::cout << "lc = " << lc << std::endl ; 
    std::cout << "npairs = " << npairs << std::endl ; 
    printf("lc = %llu\n", lc) ;

    //for(i=0; i < nframes ; i++){
     //for(j=0; j < natoms ; j++){
      //  printf("%f\t%f\t%f\n", coor[i][j][0],coor[i][j][1],coor[i][j][2]);
     //}
   // }
    for (z=0 ; z<npairs ; z++) {
     //   if(double(nframes) != 0.0){
            printf("%f\n", dist[z]/double(nframes));
      //      printf("%f\n", dist[z]/double(nframes)) ;
     //   }
        //sstream << dist[z]/double(nframes) << endl;
    }
    //outfile << sstream ;
 */     

    return ;
}
