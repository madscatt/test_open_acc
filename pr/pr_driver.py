import sasmol.sasmol as sasmol
import numpy
import pr_parallel as pr_parallel

try:
    import prc as prc
    flag = True
except:
    print 'could not import prc'
    flag = False

def calculate_pr(pdbfile, dcdfile):

    m = sasmol.SasMol(0)
    m.read_pdb(pdbfile)
    m.read_dcd(dcdfile)

    coor = m.coor()
    nf = m.number_of_frames()

    print 'nf = ', nf
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
    dist = pr_parallel.pr_parallel(coor,nf,natoms)
    
    print '\nback in python\n\n'
    
    outfile = open('dist.txt','w')
    for val in dist:
        outfile.write('%f\n' % val) 

    outfile.close()

if __name__ == "__main__":

    #pdbfile = 'ten_mer.pdb'
    pdbfile = 'n.pdb'
    #dcdfile = 'ten_mer.dcd'
    #dcdfile = 'n1000.dcd'
    #dcdfile = 'n10000.dcd'
    dcdfile = 'n200.dcd'

    import time
    start_time = time.time()
    calculate_pr(pdbfile, dcdfile)
    elapsed_time = time.time() - start_time
    print 'elapsed time = ', elapsed_time
    nframes = 1000
    print 'seconds per frame (1000) = ', elapsed_time/nframes

