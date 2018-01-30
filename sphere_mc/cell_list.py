import sasmol.sasmol as sasmol
import math
import numpy
import sys

def get_icell(ix, iy, iz, ncell_1d):

    xpart = ((ix + ncell_1d) % ncell_1d) 
    ypart = ((iy + ncell_1d) % ncell_1d) * ncell_1d
    zpart = ((iz + ncell_1d) % ncell_1d) * ncell_1d * ncell_1d

    return xpart + ypart + zpart

def make_neighbor_map(ncell_1d):

    mapsize = int(math.pow(ncell_1d,3)) * 26

    print 'setting up neighbor maps'
    print 'mapsize = ', mapsize

    map = numpy.zeros(mapsize, numpy.int)

    for iz in xrange(ncell_1d):

        for iy in xrange(ncell_1d):

            for ix in xrange(ncell_1d):

                icell = get_icell(ix, iy, iz, ncell_1d)
                imap = icell * 26

                map[imap] = get_icell(ix, iy, iz + 1, ncell_1d)                       # 1
                map[imap + 1] = get_icell(ix, iy, iz - 1, ncell_1d)                       # 2
    
                map[imap + 2] = get_icell(ix, iy + 1, iz, ncell_1d)                       # 3
                map[imap + 3] = get_icell(ix, iy - 1, iz, ncell_1d)                       # 4
    
                map[imap + 4] = get_icell(ix + 1, iy, iz, ncell_1d)                       # 5
                map[imap + 5] = get_icell(ix - 1, iy, iz, ncell_1d)                       # 6

                map[imap + 6] = get_icell(ix + 1, iy + 1, iz, ncell_1d)                   # 7
                map[imap + 7] = get_icell(ix + 1, iy - 1, iz, ncell_1d)                   # 8
                map[imap + 8] = get_icell(ix - 1, iy + 1, iz, ncell_1d)                   # 9
                map[imap + 9] = get_icell(ix - 1, iy - 1, iz, ncell_1d)                   # 10

                map[imap + 10] = get_icell(ix + 1, iy, iz + 1, ncell_1d)                   # 11
                map[imap + 11] = get_icell(ix + 1, iy, iz - 1, ncell_1d)                   # 12
                map[imap + 12] = get_icell(ix - 1, iy, iz + 1, ncell_1d)                   # 13
                map[imap + 13] = get_icell(ix - 1, iy, iz - 1, ncell_1d)                   # 14

                map[imap + 14] = get_icell(ix + 1, iy + 1, iz + 1, ncell_1d)               # 15
                map[imap + 15] = get_icell(ix + 1, iy + 1, iz - 1, ncell_1d)               # 16
                map[imap + 16] = get_icell(ix + 1, iy - 1, iz + 1, ncell_1d)               # 17
                map[imap + 17] = get_icell(ix + 1, iy - 1, iz - 1, ncell_1d)               # 18
                map[imap + 18] = get_icell(ix - 1, iy + 1, iz + 1, ncell_1d)               # 19
                map[imap + 19] = get_icell(ix - 1, iy + 1, iz - 1, ncell_1d)               # 20
                map[imap + 20] = get_icell(ix - 1, iy - 1, iz + 1, ncell_1d)               # 21
                map[imap + 21] = get_icell(ix - 1, iy - 1, iz - 1, ncell_1d)               # 22

                map[imap + 22] = get_icell(ix, iy + 1, iz + 1, ncell_1d)                   # 23
                map[imap + 23] = get_icell(ix, iy + 1, iz - 1, ncell_1d)                   # 24
                map[imap + 24] = get_icell(ix, iy - 1, iz + 1, ncell_1d)                   # 25
                map[imap + 25] = get_icell(ix, iy - 1, iz - 1, ncell_1d)                   # 26


                if (ix == 0 and iy == 0 and iz == 0):
                    print 'map0 = ',map[4]
                    print 'map1 = ',map[6]
                    print 'map2 = ',map[2]
                    print 'map3 = ',map[8]

    return map


def find_box(atom, ncell_1d, x_min, y_min, z_min, cell_length, acell, cell_occupancy, atoms_in_cell, atom_number):

    found = 0 
    not_found = 0

    cell_number = 0

    for i in xrange(ncell_1d):
        this_min_x = x_min + cell_length * i ; 
        this_max_x = this_min_x + cell_length ; 
        for j in xrange(ncell_1d): 
            this_min_y = y_min + cell_length * j ; 
            this_max_y = this_min_y + cell_length ; 
            for k in xrange(ncell_1d): 
                this_min_z = z_min + cell_length * k ; 
                this_max_z = this_min_z + cell_length ; 
                if atom[0] > this_min_x and atom[0] <= this_max_x:
                    if atom[1] > this_min_y and atom[1] <= this_max_y:
                        if atom[2] > this_min_z and atom[2] <= this_max_z:
                            acell.append(cell_number)
                            cell_occupancy[cell_number] += 1
                            atoms_in_cell[cell_number].append(atom_number)
                            found += 1
                            return found, not_found
                cell_number += 1
            cell_number += 1
        cell_number += 1
    
    not_found += 1

    print "did not find a cell for atom = ", atom_number
    print 'atom[0] = ', atom[0]
    print 'atom[1] = ', atom[1]
    print 'atom[2] = ', atom[2]

    return found, not_found


def set_up_cells(mol, r_cutoff, smidge ):

    dimensions = mol.calcminmax()

    print 'dimensions = ', dimensions

    '''
    Get dimensions from coordinates, then make min and max cover a range +/- r_cutoff
    then make it larger again by "smidge" to make sure that all occupied cells are
    surrounded by cells.

    One should have a symmetrical system, but if not, the size of the system in regards
    to cell size is based on the maximum cell length encountered in x, y or z.

    '''

    x_min = math.floor(dimensions[0][0]) ; x_max = math.ceil(dimensions[1][0])
    y_min = math.floor(dimensions[0][1]) ; y_max = math.ceil(dimensions[1][1])
    z_min = math.floor(dimensions[0][2]) ; z_max = math.ceil(dimensions[1][2])

    print ; print 'real dimensions : '
    print 'x_min = ', x_min, ' : x_max = ', x_max
    print 'y_min = ', y_min, ' : y_max = ', y_max
    print 'z_min = ', z_min, ' : z_max = ', z_max

    x_min -= r_cutoff ; x_max += r_cutoff
    y_min -= r_cutoff ; y_max += r_cutoff
    z_min -= r_cutoff ; z_max += r_cutoff

    print ; print 'dimensions +/- r_cutoff : '
    print 'x_min = ', x_min, ' : x_max = ', x_max
    print 'y_min = ', y_min, ' : y_max = ', y_max
    print 'z_min = ', z_min, ' : z_max = ', z_max

    x_min *= smidge ; x_max *= smidge
    y_min *= smidge ; y_max *= smidge
    z_min *= smidge ; z_max *= smidge

    print ; print 'dimensions * smidge = ', smidge,' :' 
    print 'x_min = ', x_min, ' : x_max = ', x_max
    print 'y_min = ', y_min, ' : y_max = ', y_max
    print 'z_min = ', z_min, ' : z_max = ', z_max

    dx = x_max - x_min ; dy = y_max - y_min ; dz = z_max - z_min

    delta = max([dx,dy,dz])

    cell_length = delta/(int(delta/r_cutoff))

    print ; print 'dx = ', dx, ' : dy = ', dy, ' : dz = ', dz

    print 'r_cutoff = ', r_cutoff
    print 'cell_length = ', cell_length

    ncell_1d = int(math.ceil(delta/cell_length)) ; print 'ncell_1d = ', ncell_1d

    ncell = int(math.pow(ncell_1d, 3))

    return ncell_1d, ncell, x_min, y_min, z_min, x_max, y_max, z_max, delta, cell_length

def get_cell_list(mol, r_cutoff, smidge, **kwargs):

    debug = False

    if 'debug' in kwargs:
        debug = True

    ncell_1d, ncell, x_min, y_min, z_min, x_max, y_max, z_max, delta, cell_length = set_up_cells(mol, r_cutoff,  smidge)

    #    ncell_1d = 5 ; # use this number of cells to validate return values
 
    map = make_neighbor_map(ncell_1d)

    hoc = numpy.ones(ncell, numpy.int) * -1

    coor = mol.coor()[0]
    acell = []

    cell_occupancy = numpy.zeros(ncell,numpy.int)
    atom_number = 0

    atoms_in_cell = [ [] for _ in range(ncell)]
    found = 0 ; not_found = 0

    ll = numpy.zeros(mol.natoms(), numpy.int)

    for atom in coor:
   
        if debug:
            local_found, local_not_found = find_box(atom, ncell_1d, x_min, y_min, z_min, cell_length, acell, cell_occupancy, atoms_in_cell, atom_number)

        # for testing only

        if atom_number == 0:
            atom[0] = -(delta/2.0) + 0.1
            atom[1] = -(delta/2.0) + 0.1
            atom[2] = -(delta/2.0) + 0.1
            atom[0] = (delta/2.0) - 0.1 
            atom[1] = (delta/2.0)  - 0.1
            atom[2] = (delta/2.0)  - 0.1

        # allen & tildsley style (assumes cell length = 1)

        #icell_x = int((atom[0] + (0.5 * delta)) * 1.0/ncell_1d)
        #icell_y = int((atom[1] + (0.5 * delta)) * 1.0/ncell_1d)
        #icell_z = int((atom[2] + (0.5 * delta)) * 1.0/ncell_1d)

        # frenkel & smit style

        icell_x = int((atom[0] + (0.5 * delta))/cell_length) 
        icell_y = int((atom[1] + (0.5 * delta))/cell_length) 
        icell_z = int((atom[2] + (0.5 * delta))/cell_length) 

        icel = icell_x + (ncell_1d * icell_y) + (ncell_1d * ncell_1d * icell_z)

        if atom_number == 0:
            print 'atom[0] = ', atom[0]
            print 'atom[1] = ', atom[1]
            print 'atom[2] = ', atom[2]
            print 'atom 0 is in cell: ', icel

        if icel > ncell:
            print 'icel = ', icel
            print 'atom[0] = ', atom[0]
            print 'atom[1] = ', atom[1]
            print 'atom[2] = ', atom[2]

        ll[atom_number] = hoc[icel]
        hoc[icel] = atom_number 
         
        atom_number += 1
 
        if debug:
            found += local_found ; not_found += local_not_found

    number_of_empty_cells = 0
    number_of_occupied_cells = 0
    for i in xrange(ncell):
        if hoc[i] > 0:
            #print 'hoc[i] = ', hoc[i],
            number_of_occupied_cells += 1
        else:
            #print 'cell : ', i, ' is empty'
            number_of_empty_cells += 1

    print 'number of occupied cells = ', number_of_occupied_cells
    print 'number of empty cells = ', number_of_empty_cells


    if debug:
        for i in xrange(mol.natoms()):
            print ll[i],

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

    return map, ll, hoc


if __name__ == "__main__":

    pdbfile = 'initial_structure.pdb'

    r_cutoff = 10.0
    smidge = 1.25 # percentage to increase box size so that all occupied cells have neighbors

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdbfile)

    map, ll, hoc = get_cell_list(mol, r_cutoff, smidge, debug=True)


