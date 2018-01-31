import sasmol.sasmol as sasmol
import numpy
import math
import random
import os
import sys
import smc_parallel

sys.path.append("./")

import make_ring as make_ring
import cell_list

class parameters():

    def __init__(self):

        self.number_of_steps = 2
        #self.temperature = 1.0
        self.temperature = 0.3
        self.rho = 0.65
        self.natoms = 7352
        self.radius = 30.0

        self.sigma_11 = 1.0
        self.sigma_22 = 1.0
        self.sigma_12= 1.0

        self.epsilon_a_11 = 1.0
        self.epsilon_a_22 = 1.0
        self.epsilon_a_12 = 1.0

        self.epsilon_r_11 = 1.0
        self.epsilon_r_22 = 1.0
        self.epsilon_r_12 = 1.0

        self.r_a_11 = 1.0
        self.r_a_22 = 1.0
        self.r_a_12 = 1.0

        self.r_r_11 = 2.0
        self.r_r_22 = 2.0
        self.r_r_12 = 2.0

        self.energy = 1E99
        self.beta = 1.0     # Boltzman temperature

        self.u_long_range_1 = -self.epsilon_a_11 * pow(( self.sigma_11 / self.r_a_11 ),2.0) * math.exp(-(self.sigma_11/self.r_a_11)) ;
        self.u_long_range_2 = self.epsilon_r_11 * pow(( self.sigma_11 / self.r_r_11 ),2.0) * math.exp(-(self.sigma_11/self.r_r_11)) ;
        
        self.beta = self.u_long_range_1 + self.u_long_range_2

        self.contrast_1 = 1E-6
        self.contrast_2 = -1E-7


        self.smidge = 1.25
        self.r_cutoff = 12.0

        self.dcdfile_name = 'test.dcd'


def get_atom_id_list(mol):

    atom_id_list = []
    name = mol.name()
    for atom in mol.name():
        if atom == 'H':
            atom_id_list.append(0)
        else:
            atom_id_list.append(1)
        
    return atom_id_list

def mc_run(restartpdb):

    p = parameters()

    if restartpdb:
        print '>> reading initial configuration from pdb file'
        mol = sasmol.SasMol(0)
        mol.read_pdb(restartpdb)
    else:
        print '>> setting up initial configuration'
        mol = make_ring.ring(p.radius*2,p.natoms,p.sigma_11*2.0)

	p.natoms = mol.natoms()
	mol.center(0)
	print '>> natoms = ', p.natoms

    map, linked_list, head_of_chain, cell_length, delta, ncell_1d, ncell = cell_list.get_cell_list(mol, p.r_cutoff, p.smidge)

    coor = mol.coor()#[0]

    atom_id_list = get_atom_id_list(mol)

    print 'atom_id_list[0] = ', atom_id_list[0]
    print 'atom_id_list[1] = ', atom_id_list[1]

    smc_parallel.smc_parallel(coor,\
                              atom_id_list,\
                              cell_length,\
                              delta,\
                              ncell_1d,\
                              ncell,\
                              p.number_of_steps,\
                              p.natoms,\
                              p.temperature,\
                              p.sigma_11,\
                              p.sigma_22,\
                              p.sigma_12,\
                              p.epsilon_a_11,\
                              p.epsilon_a_22,\
                              p.epsilon_a_12,\
                              p.epsilon_r_11,\
                              p.epsilon_r_22,\
                              p.epsilon_r_12,\
                              p.r_a_11,\
                              p.r_a_22,\
                              p.r_a_12,\
                              p.r_r_11,\
                              p.r_r_22,\
                              p.r_r_12,\
                              p.beta,\
                              p.contrast_1,\
                              p.contrast_2,\
                              p.dcdfile_name)

    return

if __name__ == "__main__":

    restartpdb = False
    #restartpdb = 'initial_structure.pdb'
    import time
    start_time = time.time()
    mc_run(restartpdb)
    elapsed_time = time.time() - start_time
    print 'elapsed time = ', elapsed_time
   
    # hmm, replae this in the future; means nothing for now

    nframes = 1000
    print 'seconds per frame (1000) = ', elapsed_time/nframes

    print '\n\n>>> DONE <<<\n\n'
