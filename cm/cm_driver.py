import sasmol.sasmol as sasmol
import numpy
import cm_parallel as cm_parallel
import math
import sys

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

def get_domains(m, basis_1, basis_2):

    frame = 0

    error, domain_1_mask = m.get_subset_mask(basis_1)
    error, domain_2_mask = m.get_subset_mask(basis_2)

    domain_1_mol = sasmol.SasMol(0)
    domain_2_mol = sasmol.SasMol(0)

    error = m.copy_molecule_using_mask(domain_1_mol,domain_1_mask,frame) 
    error = m.copy_molecule_using_mask(domain_2_mol,domain_2_mask,frame) 

    return domain_1_mask, domain_2_mask, domain_1_mol, domain_2_mol

def calculate_cm(pdbfile, dcdfile, input_nbins, bin_width, basis_1, basis_2, cutoff):

    m = sasmol.SasMol(0)
    m.read_pdb(pdbfile)

    nbins = get_number_of_bins_from_structure(m, input_nbins)
    
    print 'nbins = ', nbins
    
    m.read_dcd(dcdfile)
    nf = m.number_of_frames()
    print 'nf = ', nf

    coor = m.coor()
    natoms = m.natoms()

    domain_1_mask, domain_2_mask, domain_1_mol, domain_2_mol = get_domains(m, basis_1, basis_2)

    coor_1 = domain_1_mol.coor()
    coor_2 = domain_2_mol.coor()
    
    natoms_1 = domain_1_mol.natoms()
    natoms_2 = domain_2_mol.natoms()

    print 'natoms_1 = ', natoms_1
    print 'natoms_2 = ', natoms_2

    print 'python: coor_1[0][0][0] = ', coor_1[0][0][0]
    print 'python: coor_2[0][0][0] = ', coor_2[0][0][0]

    print 'calling cm_parallel'
    dist = cm_parallel.cm_parallel(coor_1,coor_2,nf,natoms_1,natoms_2,nbins,bin_width,cutoff)
    #dist = cm_parallel.cm_parallel(coor,nf,natoms,nbins,bin_width)
    
    print '\nback in python\n\n'
    
    #outfile = open('dist.txt','w')
    #for val in dist:
    #    outfile.write('%f\n' % val) 
#
#    outfile.close()


if __name__ == "__main__":

    input_nbins = 200
    bin_width = 1.0

    cutoff = 4.0 

    pdbfile = 'nist_mab.pdb'
    
    dcdfile = 'nist_mab_10_frames.dcd'
    #dcdfile = 'xray_x2_lt_55.dcd'

    basis_1 = '(name[i][0] != "H") and (segname[i] == "HC1" and resid[i] < 219) or (segname[i] == "LC1")' # fab_1
    basis_2 = '(name[i][0] != "H") and (segname[i] == "HC2" and resid[i] < 219) or (segname[i] == "LC2")' # fab_2
    
    basis_2 = '(name[i][0] != "H") and (segname[i] == "HC1" and resid[i] > 233) or (segname[i] == "HC2" \
                   and resid[i] > 233) or (segname[i][:3] == "SUG")'               # fc

    #basis_2 = '(name[i][0] != "H") and (segname[i] == "HC1" and resid[i] < 219) or (segname[i] == "HC1" and resid[i] > 233)' # heavy chain 1 no sugar
    
    import time
    start_time = time.time()

    calculate_cm(pdbfile, dcdfile, input_nbins, bin_width, basis_1, basis_2, cutoff)

    elapsed_time = time.time() - start_time
    print 'elapsed time = ', elapsed_time
    nframes = 1000
    print 'seconds per frame (1000) = ', elapsed_time/nframes

