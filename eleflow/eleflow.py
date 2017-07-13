'''                                  ELEFLOW.py
				Author: Subodh Joshi

This code computes analytical phi and psi for elementary flows. Superposition is carried out to simulate complex flows. 
Some instructions:
	1. use paraview or Mayavi for visualization. output file is in vtk format.
	2. You might want to adjust the strengths of elementary flows such that it is visible in post processing. Otherwise just center points of singularities will be visible and rest flow features will be lost.
	3. plot contours and streamlines for better visualization
	4. Use 'Rainbow Desaturated' color map for better visualization 
	5. Refer Fundamentals of Aerodynamics (J. Anderson) for theory
	6. use very small value for double / source /sink strenght (something like 0.1 ) and of the order of 10-100 for uniform flow to see the effects properly

To run the code-
type 

python eleflow.py 

in terminal. Follow the instructions

It is meant to be for academic exercise purpose. Feel free to change, delete, share, use the part of / all of code. Applicable only for incompressible, inviscid and irrotational flows

Thank you. 
'''
#__________________________________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________________________________


# numpy is a python module for computations (NUMerical PYthon ). Import it with pseudo-name np
import numpy as np
#______________________________________________________________________________________________

sin = np.sin;    # Its tedius writing np.sin every time so assign np.sin to a word sin here you can put your name also! (thats the beauty of python)
cos = np.cos;

#______________________________________________________________________________________________
# Datastructure
# google for what class is (OOP). refer- 'http://openbookproject.net/thinkcs/python/english3e/'
class pt:
	def __init__(self,x, y):
		# attributes
		self.x =x;
		self.y =y;
		self.phi = 0;
		self.psi = 0;
		self.u = 0;
		self.v = 0;
		pass;
	pass;

#______________________________________________________________________________________________
# Functions
#______________________________________________________________________________________________

# Google for VTK file format. This function takes a list of points (each member of which is an instance of class pt)
def write_vtk(points):
	# You need not worry about this function at this stage. Its just producing the output file...
	with open('eleflow' +'.vtk','w') as file:
		print >> file, '# vtk DataFile Version 5.8';
		print >> file, 'Elementary flows';
		print >> file, 'ASCII';
		print >> file, 'DATASET UNSTRUCTURED_GRID';
		#print >> file, 'DIMENSIONS', len(x), len(y), 1;
		print >> file, 'POINTS', len(x)*len(y), 'FLOAT';
		
		for i in range(len(x)):
			for j in range(len(y)):
				print >> file, x[i], y[j], 0;
				pass;
			pass;
		"""
		print >> file, 'DATASET STRUCTURED_POINTS';
		print >> file, 'DIMENSIONS', len(x), len(y), 1;
		print >> file, 'ASPECT_RATIO 1 1 1';
		print >> file, 'ORIGIN 0 0 0';
		"""
		print >> file, 'CELLS', len(points), 5*len(points);
		for i in range(len(x)-1):
			for j in range(len(y)-1):
				I = len(y)*i + j;
				print >> file, 4, I, I+1, I+len(y)+1, I+len(y);
				pass;
			pass
		print >> file, 'CELL_TYPES', len(points);
		for p in points:
			print >> file, 9;
			pass;

		print >> file, 'CELL_DATA', len(points) ;
		print >> file, 'SCALARS phi FLOAT';
		print >> file, 'LOOKUP_TABLE default';
		for p in points:
			print >> file, p.phi;
			pass;
		print >> file, 'SCALARS psi FLOAT';
		print >> file, 'LOOKUP_TABLE default';
		for p in points:
			print >> file, p.psi;
			pass;

		print >> file, 'VECTORS velocity FLOAT'
		for p in points:
			print >> file, p.u, p.v, 0.0;
			pass;
		pass;
	pass;
#______________________________________________________________________________________________
#______________________________________________________________________________________________


# an empty list 
points = [];

# n = number of points along each axis
n = 300;

# domain is (0,1)X(0,1) square
x = np.linspace(0,1,n); # what is linspace??   
y = x.copy();  # find out why I did not write y=x here. 


#______________________________________________________________________________________________
# Initialization
# google for for loop syntax of python. find out what len does.
for i in range(len(x)-1):
	for j in range(len(y)-1):
		# here an instance of class pt is created and immediately appended to list points
		points.append(pt((x[i]+x[i+1])/2., (y[j]+y[j+1])/2.));
		pass;
	pass;

# now points is not empty anymore. It has all instances of class pt. each instance represents a geometrical point in the space. 
#______________________________________________________________________________________________

ip = -1; # Initial value 
print 'Elementary flows:- \n 1. Uniform flow \n 2. Source \n 3. Sink \n 4. Doublet \n 5. Vortex \n 0. Thats enough...  \n';
while ip != 0:
	ip = input('Enter your choice\n' )
	i += 1;

	
	#______________________________________________________________________________________________
	if (ip == 1): # uniform flow
		theta = input('enter direction of uniform flow (theta in degrees) \n');
		V = input("Enter speed of the flow \n");
		theta = theta * np.pi/ 180; # convert into radians
		print 'Proccessing... Please wait... ';
		for p in points:
			# iterate over the list
			# use formulae for uniform flow. (in random direction)
			p.phi += p.x * V*cos(theta) + p.y * V*sin(theta);
			p.psi += p.y * V*cos(theta) - p.x * V*sin(theta);
			p.u += V*cos(theta);
			p.v += V*sin(theta);
			pass;
		pass;
	
	#______________________________________________________________________________________________
	elif (ip == 2): # Source
		Lambda = input('enter strength of the source ');
		print 'enter x,y coordinates of the singularity \n'
		X = input('X coordinate: ');
		Y = input('Y coordinate: ');
		print 'Proccessing... Please wait... ';
		for i in range (len(points)):
			p = points[i];

			if (p.x == X and p.y == Y): # singularity point. we don't know what happens at source center
				# so values of neighbor point taken here
				p.u += points[i-1].u;
				p.v += points[i-1].v;
				p.phi += points[i-1].phi; 
				p.psi += points[i-1].psi;
			else:
				# source computations
				r = np.sqrt((p.x - X)**2 + (p.y - Y)**2);
				theta = np.arctan2((p.y-Y), (p.x-X)); # why arctan2? find out.
				
				p.phi += Lambda / (2*np.pi) * np.log(r);
				p.psi += Lambda / (2*np.pi) * theta;

				vr = Lambda / (2*np.pi*r);
				vt = 0;
				
				p.u += vr*cos(theta);
				p.v += vr*sin(theta);
				pass;
			pass;
		pass;

	#______________________________________________________________________________________________
	elif (ip == 3): # Sink
		Lambda = input('enter strength of the sink (magnitude) ');
		Lambda = Lambda * -1;
		print 'enter x,y coordinates of the singularity \n'
		X = input('X coordinate: ');
		Y = input('Y coordinate: ');
		print 'Proccessing... Please wait... ';
		for p in points:
			if (p.x == X and p.y == Y):
				p.u += points[i-1].u;
				p.v += points[i-1].v;
				p.phi += points[i-1].phi; 
				p.psi += points[i-1].psi;
			else:
				r = np.sqrt((p.x - X)**2 + (p.y - Y)**2) + 10**(-12);
				theta = np.arctan2((p.y-Y), (p.x-X));
				
				p.phi += Lambda / (2*np.pi) * np.log(r);
				p.psi += Lambda / (2*np.pi) * theta;

				vr = Lambda / (2*np.pi*r);
				vt = 0;
				
				p.u += vr*cos(theta);
				p.v += vr*sin(theta);
				pass;
			pass;
		pass;

	#______________________________________________________________________________________________
	elif (ip == 4): # Doublet
		Lambda = input('enter strength of the doublet ');
		print 'enter x,y coordinates of the singularity \n'
		X = input('X coordinate: ');
		Y = input('Y coordinate: ');
		print 'Proccessing... Please wait... ';
		for p in points:
			if (p.x == X and p.y == Y):
				p.u += points[i-1].u;
				p.v += points[i-1].v;
				p.phi += points[i-1].phi; 
				p.psi += points[i-1].psi;
			else:
				r = np.sqrt((p.x - X)**2 + (p.y - Y)**2) + 10**(-12);
				theta = np.arctan2((p.y-Y), (p.x-X));
				
				p.phi += Lambda / (2*np.pi) * cos(theta) / r;
				p.psi += -Lambda / (2*np.pi) * sin(theta) / r;

				vr = -Lambda / (2*np.pi*r**2)*cos(theta);
				vt = -Lambda / (2*np.pi*r**2)*sin(theta);
				
				p.u += vr*cos(theta) - vt*sin(theta);
				p.v += vr*sin(theta) + vt*cos(theta);
				pass;
			pass;
		pass;
#______________________________________________________________________________________________
#______________________________________________________________________________________________
# now all phis and psis have been computed and superposition at points is performed. Only plotting remains.


print "\nwriting vtk file... ";

write_vtk(points); # function vtk called to write output. 

#______________________________________________________________________________________________
#______________________________________________________________________________________________
'''					END							'''


'''
If you want to learn better, nothing better than doing it yourself!
Following questions are for your own practice. You can choose to not to attempt! 

1. Although the list says vortex, I have not written a function for vortex. You can add that here
2. What condition would you give at the vortex center?
3. write a function for computations of Cp. call that at appropriate time to compute Cp in the domain. Validate with circular cylinder case
4. Exercise,
           Take two sources at 0.5,a and 0.5.-a (a<0.5)(tutorial question). Simulate using this code. See whether you can simulate wall effect
	   Using uniform flow and doublet, simulate flow past circular cylinder
	   using a source at (0.5,0.5), take a sink at (a,0.5) reduce 'a' to a very small distance. See if same effect as doublet is seen

5. Have fun! 

'''

# the code has been written in hurry. So mistakes can be there. please point out the mistakes and even better, try to correct.
# Thank you and cheers
# - Subodh 
