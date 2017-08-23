''' written by:  Subodh Joshi'''
import numpy as np
import grid as g

'''
#_____________________________________________________________________________________________________________________________________
Geometrical Functions:
#_____________________________________________________________________________________________________________________________________

'''
def length(one,two):
	# one and two are objects of class node. While finding length of any edge, enter  length(edge#.one, edge#.two) to get the length
	y2 = two.y
	y1 = one.y
	x2 = two.x
	x1 = one.x
	return np.sqrt((y2-y1)**2 + (x2 - x1)**2);

def area(nodes):
	#nodes should be a list of points each an instance of class node . Refer grid.py; 
	if len(nodes) == 3:
		#It is a triangle
		a = length(nodes[0],nodes[1])
		b = length(nodes[1],nodes[2])
		c = length(nodes[2],nodes[0])

		'''
		a = faces[0].length;
		b = faces[1].length;
		c = faces[2].length;

		'''
		s = (a + b + c)/2.
		A = np.sqrt(s*(s-a)*(s-b)*(s-c));

	elif len(nodes) == 4:
		# quadrilateral
		a = length(nodes[0],nodes[1])
		b = length(nodes[1],nodes[2])
		c = length(nodes[2],nodes[0])

		s = (a + b + c)/2.
		Aone = np.sqrt(s*(s-a)*(s-b)*(s-c));

		aa = length(nodes[2],nodes[3])
		bb = length(nodes[3],nodes[0])
		cc = length(nodes[0],nodes[2])

		ss = (aa + bb + cc)/2.
		Atwo = np.sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));
		'''

		a = faces[0].length;
		b = faces[1].length;
		c = length(faces[1].b, faces[0].a);
		s = (a + b + c)/2.
		Aone = np.sqrt(s*(s-a)*(s-b)*(s-c));

		aa = faces[2].length;
		bb = faces[3].length;
		cc = length(faces[3].b, faces[2].a);
		
		if (np.abs(cc - c) > 10**(-3)):
			print 'check the orientation of the faces in cell'
			pass;
		ss = (aa + bb + cc)/2.
		Atwo = np.sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));
		'''
		A = Aone + Atwo;

	return A;



'''
#_____________________________________________________________________________________________________________________________________
Vector Operations:
#_____________________________________________________________________________________________________________________________________
'''
class Vector:
	def __init__(self,x,y):
		self.x = x;
		self.y = y;
		pass;
	pass;

def dot(A,B):
	#A and B are objects of class vector (refer grid.py)
	return A.x * B.x + A.y * B.y

def norm(A):
	# A is an object of class vecotr
	return np.sqrt(A.x*A.x + A.y * A.y);


def findphi(face):
	# finds angle made by face normal with positive x axis;
	return np.arctan2(face.normal.y, face.normal.x);

def rotate(vect, angle):
	# vect is an object of class vector
	# angle is in radians 
	# returns rotated vector which is an object of class Vector
	un = vect.x *np.cos(angle) + vect.y*np.sin(angle);
	ut = vect.y *np.cos(angle) - vect.x*np.sin(angle);
	return Vector(un,ut);

def rotback(vect, angle):
	# similar like rotate, but in reverse direction (clockwise)
	# same effects can be achieved by calling rotate and giving angle argument as -angle
	u = vect.x * np.cos(angle) - vect.y * np.sin(angle);
	v = vect.x * np.sin(angle) + vect.y * np.cos(angle);
	return Vector(u,v);

'''
#_____________________________________________________________________________________________________________________________________
Euler equations related functions
#_____________________________________________________________________________________________________________________________________
'''

''' Flux computations '''
def find_flux(Q):
	# Q is a numpy array of shape (1,3)
	rho = Q[0]
	u = Q[1] / (float(rho) + 10**-6)
	v = Q[2] / (float(rho) + 10**-6)
	e = Q[3]
	p = (0.4) * (e - 1/2.*rho*(u**2+v**2));

	E = np.array([rho*u, rho*u**2+p, rho*u*v, u*(p + e)]);
	F = np.array([rho*v, rho*u*v, rho*v**2 + p, v*(p+e)]);

	return Vector(E,F);
#___________________________________________________________________________________________________
def directed_flux(flux, normal):
	# flux is an object of class vector with flux.x and flux.y as numpy arrays storing E and F fluxes. 
	# normal is the face normal. Again a numpy array storing normals at the face. 
		
	F = np.zeros(4)
	for i in range(len(flux.x)):
		F[i] = flux.x[i]*normal.x + flux.y[i]*normal.y;
		pass;
	return F;

#___________________________________________________________________________________________________

# Functions related to Eigenvalues and Eigenvectors defined here

def eival(U,normal):
	# Returns eigenvalues of the jacobian matrix dF / dU 
	# np.array returned: ([EI1,  EI2, EI3]);
	# U is a numpy array of the shape (1,3) it is the vector of conserved variables. 

	rho = float(U[0]) + 10**-6;
	u   = U[1] / rho;
	v   = U[2] / rho;
	e   = U[3];
	p   = (1.4-1) * (e - 1/2.*rho *(u**2+v**2));
	h   = (e + p) / rho;
	a   = np.sqrt((1.4 - 1) * (h - 1/2.*(u**2+v**2)));
	
	un = u*normal.x + v*normal.y;
	k1 = un - a;
	k2 = un;
	k3 = un;
	k4 = un + a;

	return np.array([k1,k2,k3,k4]);


def eivect(U, normal):
	# Returns eigenvectors of the jacobian matrix dF / dU 
	# np.matrix L and R returned;
	# U is a numpy array of the shape (1,3) it is the vector of conserved variables. 

	rho = float(U[0]) + 10**-6;
	u   = U[1] / rho;
	v   = U[2] / rho;
	e   = U[3];
	p   = (1.4-1) * (e - 1/2.*rho *(u**2+v**2));
	h   = (e + p) / rho;
	a   = np.sqrt((1.4 - 1) * (h - 1/2.*(u**2+v**2)));
	

	nx = normal.x;
	ny = normal.y;
	lx = -ny;
	ly = nx;

	un = u*normal.x + v*normal.y;
	ul = u*lx + v*ly
	
	# Matrix of right Eigenvectors
	R   = np.matrix([ [  1      ,     1     ,        0 ,    1       ],
			  [ u-a*nx  ,     u     ,        lx,    u+a*nx  ],
			  [ v-ny*a  ,     v     ,        ly,    v+ny*a  ],
			  [ h-un*a  ,  1/2.*(u**2+v**2), ul,    h+un*a  ]]);
	
	# Matrix of left Eigenvectors
	L   = R.I;
	
	return L, R;


#___________________________________________________________________________________________________

def ROE(U, V):
	# Eigenvalue computation based on Roe's averaging at the face
	# U and V are the conserved variable vectors on left and right of the face. 
	# Left state of the face
	rhol = float(U[0]) + 10**-6;
	ul   = U[1] / float(rhol);
	vl   = U[2] / float(rhol);
	el   = U[3];
	pl   = 0.4 * (el - 1/2.*rhol*(ul**2+ vl**2));
	hl   = (el + pl) / float(rhol);
	
	# Right state of the face
	rhor = float(V[0]) + 10 **-6;
	ur   = V[1] / float(rhor);
	vr   = V[2] / float(rhor);
	er   = V[3];
	pr   = 0.4 * (er - 1/2.*rhor*(ur**2+ vr**2));
	hr   = (er + pr) / float(rhor);
	
	# Roe averaging
	rho = np.sqrt(rhol * rhor)
	u   = ( (ul * np.sqrt(rhol) + ur * np.sqrt(rhor)) / (np.sqrt(rhol) + np.sqrt(rhor)));
	v   = ( (vl * np.sqrt(rhol) + vr * np.sqrt(rhor)) / (np.sqrt(rhol) + np.sqrt(rhor)));
	h   = ( (hl * np.sqrt(rhol) + hr * np.sqrt(rhor)) / (np.sqrt(rhol) + np.sqrt(rhor)));
	
	e = 1/1.4 * ( h*rho + 0.4/2. * rho *( u**2 + v**2));

	uavg = np.array([ rho, rho * u, rho*v,  e]);
	
	return uavg;
	pass;
'''
#_____________________________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________________________
'''

