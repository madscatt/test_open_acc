import string,os,math,numpy
import sasmol.sasmol as sasmol

PI = numpy.pi


def get_pdb_values(m1,natoms):

        atom=[] ; index=[] ; name=[] ; loc=[] ; resname=[] ; chain=[] ; resid=[] ; rescode=[]
        x=[] ; y=[] ; z=[]
        occupancy=[] ; beta=[] ; segname=[] ; element=[] ; charge=[] ; moltype=[]

        for i in xrange(natoms):
                atom.append("ATOM  ")
                index.append(i+1)
#                name.append("C")
                loc.append(" ")
                resname.append("LIP")
                chain.append(" ")
                resid.append(i+1)
                rescode.append(" ")
                occupancy.append("  0.00")
                beta.append("  0.00")
                segname.append("LIP")
                element.append("C")
                charge.append("  ")
                moltype.append("other")
        m1.setAtom(atom) ; m1.setIndex(index) ; m1.setName(name) ; m1.setLoc(loc) ; m1.setResname(resname)
        m1.setChain(chain) ; m1.setResid(resid) ; m1.setRescode(rescode) ; m1.setOccupancy(occupancy)
        m1.setBeta(beta) ; m1.setSegname(segname) ; m1.setElement(element) ; m1.setCharge(charge)
        m1.setMoltype(moltype)
        m1.setNatoms(natoms)

        return


def ring(ring_diameter,number_of_balls,ball_diameter):
	

	gr = (1.0 + math.sqrt(5))/2.0

	n = number_of_balls

	start = -(n-1)/2
	end = (n-1)/2

	#outfile = open('ring.xyz','w')
	#outfile.write('%i\n\n' % (number_of_balls*2))
	#outfile.write('%i\n\n' % (number_of_balls))

	offset = ball_diameter / 2.0

	ring_diameter2 = ring_diameter - ball_diameter - offset
	
	ring_diameter -= offset

	dum = 0

	m1 = sasmol.SasMol(0)

        coor = numpy.zeros((1,number_of_balls-1,3),numpy.float)

	name = []
	
	for k in xrange(start,end):		

		xp = ring_diameter * math.cos(math.asin(2.0*k/n))*math.cos(2.0*PI*k/gr)		
#		xp2 = ring_diameter2 * math.cos(math.asin(2.0*k/n))*math.cos(2.0*PI*k/gr)		
		yp = ring_diameter * math.cos(math.asin(2.0*k/n))*math.sin(2.0*PI*k/gr)	
#		yp2 = ring_diameter2 * math.cos(math.asin(2.0*k/n))*math.sin(2.0*PI*k/gr)	
	
		zp = 2.0*k*ring_diameter/n
#		zp2 = 2.0*k*ring_diameter2/n

		if((dum % 2) == 0):			
			#outfile.write("%s\t%f\t%f\t%f\n" % ("H",xp,yp,zp))
			name.append("H")
		else:
			#outfile.write("%s\t%f\t%f\t%f\n" % ("D",xp,yp,zp))
			name.append("D")
		coor[0][dum][0] = xp
		coor[0][dum][1] = yp
		coor[0][dum][2] = zp
		
#		outfile.write("%s\t%f\t%f\t%f\n" % ("N",xp2,yp2,zp2))

		dum += 1

	#outfile.close()
	
        get_pdb_values(m1,number_of_balls-1)
	m1.setName(name)
        m1.setCoor(coor)
	filename = 'initial_structure.pdb'
	print '\n>>> writing initial coordinates to: '+filename
        m1.write_pdb(filename,0,"w")

	return m1

if __name__ == "__main__":

	ring_diameter = 400.0 

	circumference = ring_diameter * PI

	membrane_thickness = 50.0

	ball_diameter = circumference / membrane_thickness

	print 'ball_diameter = ',ball_diameter
	
	number_of_balls = 1001

	m1 = ring(ring_diameter,number_of_balls,ball_diameter)



