import sassie.sasmol.sasmol as sasmol
import numpy,math,random
import sys

sys.path.append("./")

import make_ring as make_ring

class parameters():

        def __init__(self):
		self.number_of_steps = 20000
		self.temperature = 3.0
		self.rho = 0.65
		self.natoms = 7001
		self.sigma_1 = 1.0
		self.sigma_2 = 1.0
		self.radius = 30.0

		self.epsilon_ab_a = 1.0
		self.epsilon_ab_r = 1.0
		self.Rab_a = 1.0	
		self.Rab_r = 2.0	
		self.sigma_ab = 1.0

def calc_dist(coor1,coor2):

	dx = coor1[0] - coor2[0]
	dy = coor1[1] - coor2[1]
	dz = coor1[2] - coor2[2]

	r = math.sqrt(dx*dx + dy*dy + dz*dz)

	return r

def energy(mol,p):

	u_long_range = 0.0 

	for i in xrange(p.natoms-1):
#		print i+1,	
#		sys.stdout.flush()

		coor_i = mol.coor()[0][i][:]	
			
		for j in xrange(i+1,p.natoms):

			coor_j = mol.coor()[0][j][:]

			r = calc_dist(coor_i,coor_j)

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
	
	mol.coor()[0][i][0] = x
	mol.coor()[0][i][1] = y
	mol.coor()[0][i][2] = z

	return

def mc_run(restartpdb):

	p = parameters()
	
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

	restartpdb = 'final_coor.pdb'
	mc_run(restartpdb)

	print '\n\n>>> DONE <<<\n\n'
