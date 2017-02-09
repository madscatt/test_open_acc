import sasmol.sasmol as sasmol
import numpy
import pr_parallel as pr_parallel

try:
    import prc as prc
except:
    print 'could not import prc'

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
    dum = pr_parallel.pr_parallel(coor,nf,natoms)


if __name__ == "__main__":

    pdbfile = 'ten_mer.pdb'
    dcdfile = 'ten_mer.dcd'

    calculate_pr(pdbfile, dcdfile)