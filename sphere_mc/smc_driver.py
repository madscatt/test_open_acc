import sasmol.sasmol as sasmol
import numpy,math,random
import os,sys

#import smc_parallel

sys.path.append("./")

import make_ring as make_ring

class parameters():

        def __init__(self):
		self.number_of_steps = 100
		self.temperature = 3.0
		self.rho = 0.65
		self.natoms = 701
		self.sigma_1 = 1.0
		self.sigma_2 = 1.0
		self.radius = 30.0

		self.epsilon_ab_a = 1.0
		self.epsilon_ab_r = 1.0
		self.Rab_a = 1.0	
		self.Rab_r = 2.0	
		self.sigma_ab = 1.0
		self.energy = 1E99
                self.beta = 1.0 # Boltzman temperature

def calc_dist(coor1,coor2):

	dx = coor1[0] - coor2[0]
	dy = coor1[1] - coor2[1]
	dz = coor1[2] - coor2[2]

	r = math.sqrt(dx*dx + dy*dy + dz*dz)

	return r

def energy(mol,p):

        ''' 
        method to calculate energy.  If particles overlap then it returns
        a boolean False, otherwise it returns the sum of the long-range energy
        '''

	u_long_range = 0.0 

	for i in xrange(p.natoms-1):
#		print i+1,	
#		sys.stdout.flush()

		coor_i = mol.coor()[0][i][:]	
			
		for j in xrange(i+1,p.natoms):

			coor_j = mol.coor()[0][j][:]

			r = calc_dist(coor_i,coor_j)

                        if r < p.sigma_1 + p.sigma_2:
	                    return False
                
			u_long_range_1 = -p.epsilon_ab_a * ( ( p.sigma_ab / p.Rab_a )**2.0 ) * math.exp(-(r/p.Rab_a))
			u_long_range_2 =  p.epsilon_ab_r * ( ( p.sigma_ab / p.Rab_r )**2.0 ) * math.exp(-(r/p.Rab_r))

			u_long_range += u_long_range_1 + u_long_range_2


	return u_long_range

def surface_move(mol,p,i):

#		x = R sin(th) cos(phi)
#		y = R sin(th) sin(phi)
#		z = R cos(th)

	max_dtheta = 0.01
	max_dphi = 0.01

	x = mol.coor()[0][i][0]
	y = mol.coor()[0][i][1]
	z = mol.coor()[0][i][2]

	r = math.sqrt(x*x+y*y+z*z)

#	theta = math.acos(z/r)
#	phi = math.atan(y/x)
	
#	ran_theta = max_dtheta*random.uniform(-1.0,1.0)
#	ran_phi = max_dphi*random.uniform(-1.0,1.0)

#	dx = r * math.sin(theta+ran_theta)*math.cos(phi+ran_phi)
#	dy = r * math.sin(theta+ran_theta)*math.sin(phi+ran_phi)
#	dz = r * math.cos(theta+ran_theta)

	max_disp = 0.05
	dx = max_disp * random.uniform(-1.0,1.0)
	dy = max_disp * random.uniform(-1.0,1.0)
	dz = max_disp * random.uniform(-1.0,1.0)

	x += dx ; y += dy ; z += dz
	norm = math.sqrt(x*x + y*y + z*z)
	x *= r/norm
	y *= r/norm
	z *= r/norm

        accepted = False

        u_long_range = energy(mol, p)
        
        if u_long_range:
            if u_long_range < p.energy:
                accepted = True
            else:
                delta_energy = u_long_range - p.energy
                boltz = math.exp(-p.beta * delta_energy)
                ran = random.random()
                if ran < boltz:
                    accepted = True

        if accepted:
            p.energy = u_long_range
	    mol.coor()[0][i][0] = x
	    mol.coor()[0][i][1] = y
	    mol.coor()[0][i][2] = z
	
        return

def mc_run(restartpdb):

	p = parameters()
        
        if restartpdb:
	    print '>> reading initial configuration from pdb file'
            if restartpdb:
                mol = sasmol.SasMol(0)
                mol.read_pdb(restartpdb)
        else:
            print '>> setting up initial configuration'
        
	    mol = make_ring.ring(p.radius*2,p.natoms,p.sigma_1*2.0)

	p.natoms = mol.natoms()
	mol.center(0)
	print '>> natoms = ',p.natoms

	dcdoutfile = mol.open_dcd_write("traj.dcd")

#	u_long_range = energy(mol,p)

#	print 'energy = ',u_long_range

	for step in xrange(p.number_of_steps):
		print step,
		sys.stdout.flush()
		
		for i in xrange(p.natoms):
				
			surface_move(mol,p,i)
		mol.center(0)
		mol.write_dcd_step(dcdoutfile,0,step)

	mol.close_dcd_write(dcdoutfile)
	mol.write_pdb("final_coor.pdb",0,"w")

	return

if __name__ == "__main__":

	restartpdb = False
	restartpdb = 'initial_structure.pdb'
    import time¬
    start_time = time.time()¬
    mc_run(restartpdb)
    elapsed_time = time.time() - start_time¬
    print 'elapsed time = ', elapsed_time¬
   
    # hmm, replae this in the future; means nothing for now

    nframes = 1000¬
    print 'seconds per frame (1000) = ', elapsed_time/nframes¬





	print '\n\n>>> DONE <<<\n\n'
