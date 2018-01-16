#include <vector>

void get_distances(double ***coor, const int nframes, const int natoms, std::vector<std::vector<int> > &hist, const int nbins, const double bin_width) ;

void get_distances_0(double ***coor, const int nframes, const int natoms, std::vector<int>& hist, const int nbins, const double bin_width) ;

void get_distances_1(double ***coor, const int nframes, const int natoms, double* dist) ;

;
void get_distances_2(double ***coor, const int nframes, const int natoms, std::vector<double>& dist);

