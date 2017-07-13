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

void get_distances(double ***coor, const int nframes, const int natoms, std::vector<int>& hist, const int nbins, const double bin_width) {
//void get_distances(double ***coor, const int nframes, const int natoms, std::vector<double>& dist) {
//void get_distances(double ***coor, const int nframes, const int natoms, double* dist) {

    int i,j,k,l ;
    unsigned long long npairs ;
    double x1, y1, z1, x2, y2, z2, sdist ;
    double dx2, dy2, dz2 ;
    double this_low_bin, this_high_bin ;
    unsigned long long local_count ; 
    unsigned long long local_hist[nbins] ; 
    unsigned long long z ; 

    /// temporary printing for debugging
    //
    //
    std::ostringstream sstream;
    std::string remark = "#hello";
    const std::string filename = "dum_oacc.txt";
    
    //std::ofstream output_file(filename.c_str()) ;
    std::ofstream output_file ;
    std::ofstream outfile(filename.c_str()) ;
    outfile << remark << std::endl;

    //
  
    printf("oacc: %d\n", nframes) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: coor[0][0][0] = %f\n", coor[0][0][0]) ;
    npairs = (natoms * (natoms - 1))/2 ;
    double dist[npairs] ;
    double local_dist[npairs] ; 
    
    std::cout << "nframes * npairs = " << nframes*npairs << std::endl ;
    std::cout << "nf * np = " << nframes * npairs << std::endl ;
    std::cout << "nf = " << nframes << std::endl ;
    std::cout << "np = " << npairs << std::endl ;

    output_file.open("dum.txt") ;

    std::cout << "starting parallel loops" << std::endl ; 
    #pragma acc data copyin(coor[nframes][natoms][3]) copy(dist[npairs], local_dist[npairs], local_hist[nbins])
    {
    for(i=0 ; i < nframes ; i++){
        for(z=0; z < npairs ; z++) { local_dist[z] = 0.0 ; }
        for(j=0; j < nbins ; j++) { local_hist[j] = 0 ; }
        std::cout << "." << std::flush ;
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
                local_dist[local_count] = sdist ; 

                #pragma acc loop 
                for(l=0; l<nbins; l++){
                    this_low_bin = double(l)*bin_width ;
                    this_high_bin = this_low_bin + bin_width ;
                    if(sdist > this_low_bin && sdist <= this_high_bin){
                        local_hist[l] += 1 ;
                        break;
                    }
                } // end of l-loop

            } // end of k-loop
            } // end of pragma acc loop
        } // end of j-loop
        } // end of pragma acc parallel loop
        
        //for(j=0 ; j < npairs ; j++){
         //   for(k=0 ; k < nbins ; k++){
          //      this_low_bin = double(k)*bin_width ;
           //     this_high_bin = this_low_bin + bin_width ;
            //    if(local_dist[j] > this_low_bin && local_dist[j] <= this_high_bin){
             //       hist[k] += 1 ;
              //      local_hist[k] += 1 ;
               //     break ;
                //}
           // }
        //}
   
        if(i<nframes){ 
            //outfile << sstream ;
            for(k=0 ; k < nbins ; k++){
                output_file << local_hist[k] << "\n" ;
            //sstream << dist[z]/double(nframes) << endl;
            }
            output_file << "#\n" ;
            output_file.flush() ;
        }

    } // end of i-loop
    } // end of pragma acc data 

    std::cout << std::endl ;

    output_file.close();
    
    //std::vector<int> hist(nbins,0);

    std::cout << "creating histogram: 1" << std::endl ;
    for(i=0 ; i < nframes ; i++){

        std::cout << " i = " << i << std::endl ; 
 
        for(j=0 ; j < npairs ; j++){
            z = (i * npairs) + j ;
            std::cout << " z = " << z << std::endl ; 
            std::cout << " npairs = " << npairs << std::endl ; 
            for(k=0 ; k < nbins ; k++){
                this_low_bin = double(k)*bin_width ;
                this_high_bin = this_low_bin + bin_width ;
                if(local_dist[z] > this_low_bin && local_dist[z] <= this_high_bin){
                    hist[k] += 1 ;
                    break ;
                }
            }
        }
    }
    std::cout << "creating histogram: 1" << std::endl ;
    std::cout << "creating histogram: 1" << std::endl ;

    //if(nframes != 0){
    //for(k=0 ; k < nbins ; k++){
     //   hist[k] /= double(nframes) ;
    //}
    //}

//    std::cout << "creating histogram: 2" << std::endl ;

 //   //#pragma acc data copyin(dist[npairs]) copy(hist[nbins])
//    //#pragma acc parallel loop
 //   for(z=0 ; z < npairs ; z++){
  //      //#pragma acc for private(this_low_bin, this_high_bin)
   //     for(i=0 ; i < nbins ; i++){
    //        this_low_bin = double(i)*bin_width ;
     //       this_high_bin = this_low_bin + bin_width ;
      //      if(dist[z]/double(nframes) > this_low_bin && dist[z]/double(nframes) <= this_high_bin){
       //         hist[i] += 1 ;
        //        break ;
         //   }
        //}
    //}

    //for(i=0 ; i < nbins ; i++){
     //   std::cout << hist[i] << std::endl ;
    //}
    std::cout << "leaving oacc" << std::endl ;

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

    //delete[] local_dist ;
   // free(local_dist) ;

    return ;
}
