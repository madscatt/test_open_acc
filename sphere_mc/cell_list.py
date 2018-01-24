import sasmol.sasmol as sasmol
import math
import numpy


def get_cell_list(pdbfile, r_cutoff):

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdbfile)

    smidge = 1.5

    dimensions = mol.calcminmax()

    print 'dimensions = ', dimensions

    x_min = math.floor(dimensions[0][0]) ; x_max = math.ceil(dimensions[1][0])
    y_min = math.floor(dimensions[0][1]) ; y_max = math.ceil(dimensions[1][1])
    z_min = math.floor(dimensions[0][2]) ; z_max = math.ceil(dimensions[1][2])

    print 'x_min = ', x_min, ' : x_max = ', x_max
    print 'y_min = ', y_min, ' : y_max = ', y_max
    print 'z_min = ', z_min, ' : z_max = ', z_max

    dx = x_max - x_min ; dy = y_max - y_min ; dz = z_max - z_min

    dx = dx * smidge ; dy = dy * smidge ; dz = dz * smidge

    print 'dx = ', dx, ' : dy = ', dy, ' : dz = ', dz

    cell_length_x = dx /int(dx/r_cutoff) 
    cell_length_y = dy /int(dy/r_cutoff) 
    cell_length_z = dz /int(dz/r_cutoff) 

    print 'r_cutoff = ', r_cutoff
    print 'cell_length_x = ', cell_length_x
    print 'cell_length_y = ', cell_length_y
    print 'cell_length_z = ', cell_length_z

    ncell_x = int(math.ceil(dx/cell_length_x)) ; print 'ncell_x = ', ncell_x
    ncell_y = int(math.ceil(dy/cell_length_y)) ; print 'ncell_y = ', ncell_y
    ncell_z = int(math.ceil(dz/cell_length_z)) ; print 'ncell_z = ', ncell_z

    ncell = ncell_x * ncell_y * ncell_z

    hoc = numpy.zeros(ncell)

    coor = mol.coor()[0]
    acell = []

    cell_occupancy = numpy.zeros(ncell,numpy.int)
    atom_number = 0

    atoms_in_cell = [ [] for _ in range(ncell)]
    not_found = 0 
    found = 0 
    for atom in coor:
   
        continue_flag = True 
        cell_number = 0 
        while continue_flag:
            for i in xrange(ncell_x):
                this_min_x = x_min + cell_length_x * i ; #print this_min
                this_max_x = this_min_x + cell_length_x ; #print this_max
                for j in xrange(ncell_y): 
                    this_min_y = y_min + cell_length_y * j ; #print this_min
                    this_max_y = this_min_y + cell_length_y ; #print this_max
                    for k in xrange(ncell_z): 
                        this_min_z = z_min + cell_length_z * k ; #print this_min
                        this_max_z = this_min_z + cell_length_z ; #print this_max
                        if atom[0] > this_min_x and atom[0] < this_max_x:
                            if atom[1] > this_min_y and atom[1] < this_max_y:
                                if atom[2] > this_min_z and atom[2] < this_max_z:
                                    acell.append(cell_number)
                                    cell_occupancy[cell_number] += 1
                                    atoms_in_cell[cell_number].append(atom_number)
                                    found += 1
                                    continue_flag = False
                                    break
                        cell_number += 1
                    if not continue_flag:
                        break
                if not continue_flag:
                    break
            #print 'uh oh ! no home for an atom : atom number ', atom_number
            continue_flag = False
            not_found += 1
            break
        atom_number += 1

    print 'natoms = ', mol.natoms()
    print 'not found = ', not_found
    print 'found = ', found

    cell_number = 0
    for cell in atoms_in_cell:
        if len(cell) > 0:
            tmol = sasmol.Molecule_Maker(len(cell))
            this_pdb = 'temp/cell_'+str(cell_number)+'.pdb'
            tcoor=numpy.zeros((1,len(cell),3),numpy.float)
            for i in xrange(len(cell)):
                tcoor[0][i][0] = mol.coor()[0][cell[i]][0]    
                tcoor[0][i][1] = mol.coor()[0][cell[i]][1]    
                tcoor[0][i][2] = mol.coor()[0][cell[i]][2]    
            tmol.setCoor(tcoor)
            tmol.write_pdb(this_pdb,0,'w')
            cell_number += 1 


    

    return


if __name__ == "__main__":

    pdbfile = 'initial_structure.pdb'

    r_cutoff = 10.0

    get_cell_list(pdbfile, r_cutoff)


