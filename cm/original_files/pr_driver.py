import sasmol.sasmol as sasmol
import numpy
import pr_parallel as pr_parallel
import math
import sys

try:
    import prc as prc
    flag = True
except:
    print 'could not import prc'
    flag = False

def get_number_of_bins_from_structure(m, input_nbins):

    minmax = m.calc_minmax_all_steps(dcdfile)
 
    dx = minmax[1][0] - minmax[0][0]
    dy = minmax[1][1] - minmax[0][1]
    dz = minmax[1][2] - minmax[0][2]

    rmax = math.sqrt(dx*dx + dy*dy + dz*dz)

    structure_nbins = int(math.ceil(rmax / bin_width))

    if structure_nbins > input_nbins:
        nbins = structure_nbins
    else:
        nbins = input_nbins

    return nbins

def calculate_pr(pdbfile, dcdfile, input_nbins, bin_width):

    m = sasmol.SasMol(0)
    m.read_pdb(pdbfile)

    nbins = get_number_of_bins_from_structure(m, input_nbins)
    
    print 'nbins = ', nbins
    
    m.read_dcd(dcdfile)
    nf = m.number_of_frames()
    print 'nf = ', nf

    coor = m.coor()
    print 'coor[0][0] = ', coor[0][0]
    print 'coor[0][-1] = ', coor[0][-1]

    natoms = m.natoms()
    print 'natoms = ', natoms

    if flag:
        npairs = (natoms * (natoms - 1))/2
        all_distances = numpy.zeros(npairs, numpy.float32)
        for i in xrange(nf): 
            distances = prc.prc(coor[i])
            all_distances += numpy.array(distances)
            if i==0:
                print 'one_distance[0] = ',distances[0]
            print '.',
     
        print 'all_distances[0] = ', all_distances[0]/nf

    print 'calling pr_parallel'

    print 'python: coor[0][0][0] = ', coor[0][0][0]
    dist = pr_parallel.pr_parallel(coor,nf,natoms,nbins,bin_width)
    
    print '\nback in python\n\n'
    
    outfile = open('dist.txt','w')
    for val in dist:
        outfile.write('%f\n' % val) 

    outfile.close()

if __name__ == "__main__":

    nbins = 200
    bin_width = 1.0

    #pdbfile = 'ten_mer.pdb'
    #pdbfile = 'n.pdb'
    pdbfile = 'nist_mab.pdb'
    
    #dcdfile = 'ten_mer.dcd'
    #dcdfile = 'ten_mer_4591.dcd'
    #dcdfile = 'n1000.dcd'
    #dcdfile = 'n10000.dcd'
    #dcdfile = 'n200.dcd'
    dcdfile = 'xray_x2_lt_55.dcd'

    import time
    start_time = time.time()
    calculate_pr(pdbfile, dcdfile, nbins, bin_width)
    elapsed_time = time.time() - start_time
    print 'elapsed time = ', elapsed_time
    nframes = 1000
    print 'seconds per frame (1000) = ', elapsed_time/nframes

