# specific to triangle software
# Imports:

from pylab import *
import numpy as np
import sys
import funcs 
import vtkgen as v
#____________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________

# Reading 'Triangle' output files. To run it, type 		$python grid.py file_name.*  '

#____________________________________________________________________________________________________________________________

class data:
	'''
	Following files should be present in working directory,
		name.ele, name.edge, name.v.edge, name.node, name.v.node, name.neigh 
	
	class data reads all the files and creates useful python lists for all entities. This will be later read by class grid and data structure will be created.


	arguments:  option: a list containing name of all the required files; 

		    passed as sys.argv, it should contain ele, nodes, vnodes, vedges, edges files. 
		    To generate all these, use following flags while generating mesh with triangle

		    Triangle -pqcevnra name.poly
		    ra is optional. used for refinement. a number should follow a. refer triangle documentation.

	'''

	def __init__(self,option):
		
		#____________________________________________________________________________________________________________________________

		# Actual list Population takes place 
		self.elements = [];
		self.edges = [];
		self.nodes = [];
		self.vnodes = [];
		self.vedges = [];
		
		op = '';
		for o in option:
			op += o + ' ';
			pass;

		print '#--------------------------------------------------------------------------------------------------'
		print 'Checking Required Input Files... '
		assert '.ele' in op and '.v.edge' in op and '.edge' in op and '.node' in op and '.v.node' in op;
		print 'Done!'
		print '#--------------------------------------------------------------------------------------------------'

		for j in range(len(option)):
			if '.ele' in option[j]:
				self._populate_elements(option[j]);
			elif '.edge' in option[j]:
				if '.v.' in option[j]:
					self._populate_vedges(option[j]);
				else:
					self._populate_edges(option[j]);
			elif '.node' in option[j]:
				if '.v.' in option[j]:
					self._populate_vnodes(option[j]);
				else:
					self._populate_nodes(option[j]);
					pass;
				pass;
			pass;
		print '#--------------------------------------------------------------------------------------------------'
		print '#--------------------------------------------------------------------------------------------------'
		print '\nAll arrays created successfully, building datastructure...\n';
		print '#--------------------------------------------------------------------------------------------------'
		print '#--------------------------------------------------------------------------------------------------'
		pass;

	# Functions defined to generate the lists
	def _populate_edges(self,option):
		print '\nCreating an array of edges' , ' using ' , option, '...';
		cols = []
		with open(option) as file:
			for line in file:
				if '#' not in line:
					cols.append([int(z) for z in line.split()])
					pass;
				pass;
			pass;
		for i in range(len(cols)):
			if len(cols[i]) == 4:
				self.edges.append(cols[i]);
				pass;
			pass;
		print 'Done!\n'
		#print edges

	#____________________________________________________________________________________________________________________________

	def _populate_elements(self,option):
		print '\nCreating an array of elements' , ' using ' ,option,'...';
		cols = []
		with open(option) as file:
			for line in file:
				# following if condition ensures that the comments are not uppended in the list
				if '#' not in line:
					cols.append([int(z) for z in line.split()])
					pass;
				pass;
			pass;
		for i in range(len(cols)):
			# following if condition makes sure that only elements are upended and not other information
			if len(cols[i]) == 4:
				self.elements.append(cols[i]);
				pass;
			pass;
		print 'Done!\n'
		#print elements
		pass;

	# Note: Same if conditions are used in all the functions in this section. So not re-documented again. 
	#____________________________________________________________________________________________________________________________

	def _populate_vedges(self,option):
		print '\nCreating an array of voronoi-edges' , ' using ' , option, '...';
		cols = []
		with open(option) as file:
			for line in file:
				if '#' not in line:
					co = []
					for i,z in enumerate(line.split()):
						if i in range(3):
							co.append(int(z));
						else:
							co.append(float(z));
							pass;
						pass;
					cols.append(co)
					pass;
				pass;
			pass;
		for i in range(len(cols)):
			if len(cols[i]) != 2:
				self.vedges.append(cols[i]);
				pass;
			pass;
		print 'Done!\n'
		#print vedges
	#____________________________________________________________________________________________________________________________

	def _populate_nodes(self,option):
		print '\nCreating an array of nodes' , ' using ' , option, '...';
		cols = []
		with open(option) as file:
			for line in file:
				if '#' not in line:
					co = []
					for i,z in enumerate(line.split()):
						if i in range(1):
							co.append(int(z));
						elif i in range(1,3):
							co.append(float(z));
						else:
							co.append(int(z));
							pass;
						pass;
					cols.append(co)
					pass;
				pass;
			pass;
		for i in range(1,len(cols)):
			self.nodes.append(cols[i]);
			pass;
		print 'Done!\n'
		#print nodes 
	#____________________________________________________________________________________________________________________________

	def _populate_vnodes(self,option):
		print '\nCreating an array of voronoi-nodes' , ' using ' , option, '...';
		cols = []
		with open(option) as file:
			for line in file:
				if '#' not in line:
					co = []
					for i,z in enumerate(line.split()):
						if i in range(1):
							co.append(int(z));
						else:
							co.append(float(z));
							pass;
						pass;
					cols.append(co)
					pass;
				pass;
			pass;
		for i in range(1,len(cols)):
			if len(cols[i]) == 3:
				self.vnodes.append(cols[i]);
				pass;
			pass;
		print 'Done!\n'
		#print vnodes
	#____________________________________________________________________________________________________________________________


		pass;
	pass;
#____________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________

# Datastructure-elements begin
#____________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________
class Vector:
	def __init__(self,x,y):
		self.x = x;
		self.y = y;
		pass;
	pass;

class node:
	'''
	class node:
		inputs required: index, x, y, marker
		object represents a grid point in 2d 
	'''
	def __init__(self,index,x,y,marker):
		self.x = x
		self.y = y
		self.index = index;
		self.marker = marker;
		pass;

	def __str__(self):
		return 'Node index = %d,\n(x,y) = (%f,%f), \nmarker = %d'%(self.index, self.x, self.y, self.marker);
	pass;



class face:
	'''
	class face:
		inputs required: index, node1, node2, marker
		objects represent faces of the FVM grid

		attributes:
			
			* NODE 1                                     normal = (dely, -delx)
			 \    ny				     tangent= (delx, dely)	
			  \   | 
			   \  |
			    \ | 
			     \|_______ nx
			      \
			       \
			        \
				 \
				  \
				   * NODE 2
		
		Geometrical:
		0. index
		1. a,b =  node 1, node 2  : objects of class node
		2. left, right = pointer to cells on left and right : objects of class cell
		3. length = geometrical length : float
		4. normal, tangent = normal ,tangent (refer fig.) : objects of class Vector
		5. midpoint = xc,yc : object of class node with (index, xc,yc), 
				where index = f<faceindex>; e.g. f34 as index for midpoint of 34th face
				
		Variables:
		1. ql, qr = conserved variable vectors on left and right :floats
		2. Fl, Fr = F indicates dot prodcut (fnx+gny) l,r indicate left and right: floats
		3. flux = vector of numerical flux at the face :float
		4. self.un = normal velocity to face:float
	'''

	def __init__(self,index,a,b,marker):
		'''location details'''
		self.index = index;
		self.a = a;
		self.b = b;
		self.marker = marker;
		self.left = None;
		self.right = None;
		'''geometry and var details'''
		self._setup_face();
		pass;

	def _setup_face(self):
		'''geometry details'''

		self.length = funcs.length(self.a,self.b);
		
		delx = self.b.x - self.a.x;
		dely = self.b.y - self.a.y;

		norm = np.sqrt(delx**2 + dely**2);
		self.tangent = Vector(delx / norm ,dely / norm);
		self.normal = Vector(dely / norm, -delx / norm);
		
		'''
		# This might need in SV

		xc = (self.b.x + self.a.x)/2.
		yc = (self.b.y + self.a.y)/2.

		self.midpoint = node('f%d'%(self.index),xc,yc,'m');
		'''

		'''var details'''	

		self.ql = None; # q : conserved variable vector to be given as np.array.
		self.qr = None;
		self.Fl = None; # capital F : fnx + gny   where f and g are np.arrays
		self.Fr = None;
		self.flux = np.zeros(4,dtype='float'); # np.array
		pass;

	def __str__(self):
		a = '\nface index = %d'%(self.index);
		if self.marker != 0:
			b = '\nLeft Cell =  %d <--- FACE(%d) ---> Right Cell = %s'%(self.left.index,self.index,self.right);
		else:
			b = '\nLeft Cell =  %d <--- FACE(%d) ---> Right Cell = %d'%(self.left.index,self.index,self.right.index);
		c = '\nTangent [%f,%f], Normal [%f,%f]'%(self.tangent.x,self.tangent.y,self.normal.x,self.normal.y);
		d = '\nNode1 == %d , Node2 == %d'%(self.a.index,self.b.index);
		e = '\nmarker = %d'%(self.marker);
		return a+b+c+d+e;
	pass;



class cell:
	'''
	class cell:
		inputs required: index, list of nodes
		objects represent faces of the FVM grid

		attributes:
			
			                         *  NODE3
						/\
					       /  \
					      /    \
					     /      \
					    /        \
					   /          \
					  /            \
					 /              \
					/                \
				NODE1  *__________________* NODE2


		Geometrical:
		0. index
		1. nodes  =  [node 1, node 2,...]  : list containing objects of class node
		2. faces  =  [face1, face2, ... ]  : list containing objects of class face
		3. area   = geometrical area : float
		4. center = xc,yc : object of class node with (index, xc,yc), 
				where index = c<faceindex>; e.g. c34 as index for midpoint of 34th cell
				
		Variables:
		1. q = conserved variable vector  :float
		2. f,g = fluxes  : floats 
		3. flux = Vector(f,g)
		4. R = Residual = residuals at cell center : float
		5. u,v = x and y components of convective velocities : floats
		6. vel = Vector(u,v);
		others:
		1. marker
		2. methods:  a)  belongsto := to be used for subcell formation
			   

	'''
	def __init__(self,index,nodes):
		self.index = index;
		self.nodes = nodes;	 
		
		self.faces = [];
		self.area = funcs.area(self.nodes);
		
		x = []
		y = []
		for i,Node in enumerate(self.nodes):
			x.append(Node.x)
			y.append(Node.y)
			pass;
		self.center = node('c%d'%(self.index),sum(x)/float(len(x)), sum(y)/float(len(y)),0);

		self.q = np.zeros(4); # np.array of conserved variables
		self.Q = np.zeros(shape=(5,4));   # 5 runge time levels of 4 conserved variable each

		self.flux = Vector(None,None);

		self.vel = Vector(None,None);


		self.R = np.zeros(shape=(4,4));  # 

		self.marker = 0 # indicating the parent element

		'''Assigning vtk types to the cells '''
		if len(self.nodes) == 3:
			self.Type = 5;
		elif len(self.nodes) == 4:
			self.Type = 9;
		else:
			raise NameError('Polygon with %d nodes has not been considered yet!! '%(len(self.nodes)));
		


		pass;

	def __str__(self):
		a = 'Cell index = %d\n'%(self.index);
		'''
		b = 'Faces = [%s]'%','.join(map(str,self.faces));
		'''
		faces =[]
		for face in self.faces:
			faces.append('%d'%(face.index));
			pass;

		nodes = []

		for node in self.nodes:
			nodes.append('%d'%(node.index));
			pass;

		b = 'Faces = [%s]'%','.join(map(str,faces));
		c = '\nNodes = [%s]'%','.join(map(str,nodes));
		d = '\nCenter = [%f,%f]'%(self.center.x,self.center.y);
		
		e = None;
		if self.marker == 0:
			e = 'Parent'
		else:
			e = 'Subcell'
			pass;

		e = '\nType = ' + e
		f = '\nArea = %f'%(self.area);
		return a + b + c + d + f + e;
	pass;



class grid:
	'''
	Attributes,
	nodes, cells, faces
	#__________________________________________________________________________________
	arguments:  option: a list containing name of all the required files created by Triangle
	#__________________________________________________________________________________

	this class creates the datastructure. 
	example,
	#__________________________________________________________________________________

	create an object of  class. say G.

	then to see all the nodes (which is a list of objects of class node)

	for node in G.nodes:
		print node;
		pass;

	to access all information of any node, 
	node = G.nodes[i]
	node.x , node.y, node.marker, node.index will give all the information 
	Refer class node
	#__________________________________________________________________________________

	to see all the faces (which is a list of objects of class face), 
	(Refer class face)
	for face in G.faces:
		print face;
		pass;
	
	to access all the info,
	face = G.faces[i]
	print face.a, face.b, face.index, face.normal, face.tangent, face.flux, face.ql, face.qr, face.Fl, face.Fr, face.midpoint.x, face.midpoint.y, face.left (which will be an object of class cell), face.right (again an object of class cell);

	#__________________________________________________________________________________
	similarly G.cells(list containing objects of of class cell). Each cell in G.cells will have a list of faces (G.cells[i].faces), which will contain objects of class face making that cell. For each of those objects, left / right will indicate this cell. 
	#__________________________________________________________________________________
	
	
	'''
	
	def __init__(self,option):
		self.name = option[1].split('.',1)[0]; #interesting! splits the name on first occurance of '.' i.e. crates a list with two parts. eg. if the meshfile is name.msh, the above line will do ['name','msh'] and will take 'name' only.

		self.data = data(option); # Creating instance of class data (also named data). i.e. self.data in this case.
		self.order = 0
		# its attrebutes will be nodes,edges,vedges,vnodes,elements

		self.nodes = [];
		self.cells = [];
		self.faces = [];

		print '\nBuilding nodes...\n'
		self._build_nodes();
		print '\nBuilding faces...\n'
		self._build_faces();
		print '\nBuilding cells...\n'
		self._build_cells();
		print '\nConnecting faces to cells...\n'
		self._connect_f_to_c()

		self.nnodes = len(self.nodes);
		self.nfaces = len(self.faces);
		self.ncells = len(self.cells);
		self.nsubcells = 1;

		print '#--------------------------------------------------------------------------------------------------'
		print '\nDone!'
		print '\nDatastructure is created.'
		print '#--------------------------------------------------------------------------------------------------'
		print '#--------------------------------------------------------------------------------------------------'

		pass;

	def _build_nodes(self):
		for Node in self.data.nodes:
			# note: Node is used instead of node, because later is a class
			n = node(Node[0],Node[1],Node[2],Node[3]);
			self.nodes.append(n);
			pass;
		pass;

	def _build_faces(self):
		for edge in self.data.edges:
			f = face(edge[0],self.nodes[edge[1]],self.nodes[edge[2]],edge[3]);
			self.faces.append(f);
			pass;
		pass;

	def _build_cells(self):
		for ele in self.data.elements:
			Nodes = [self.nodes[ele[1]],self.nodes[ele[2]],self.nodes[ele[3]]];
			e = cell(ele[0],Nodes);
			self.cells.append(e);
			pass;
		pass;

	def _connect_f_to_c(self):
		'''
			1. check face index. Go to same index (row) in voronoi edge chart. The nodes of voronoi edge are actually indices of left and right element of original face. 
			2. refer: http://www.cs.cmu.edu/~quake/triangle.topo.html
			3. Triangle creates edges in a manner so that the None cell always lies on the right side of an edge. 
		'''
		for i,face in enumerate(self.faces):
			vedge = self.data.vedges[i];
			if len(vedge) == 3:
				# both neighbors exist;
				l = vedge[1];
				r = vedge[2]; # indices of voronoi nodes making vedge. these will be triangles neighbor to delaunay edge i

				face.left = self.cells[l];
				face.right = self.cells[r];

				# immediately, we append this face in list of faces of left and right cells. 
				face.left.faces.append(face);
				face.right.faces.append(face);
			else:
				l = vedge[1];
				face.left = self.cells[l];
				face.right = None;
				face.left.faces.append(face);
				pass;
			pass;
		pass;



class Initialize:
	def __init__(self,meshfile):
		
		self.grid = grid(meshfile);
		print 'Initializing... \n'
		self._initialize();
		self._write_vtk();
		pass;
	def _initialize(self):
		grid = self.grid;

		cells = grid.cells;
		faces = grid.faces;
		nodes = grid.nodes;
		for cell in cells:
			x = cell.center.x
			y = cell.center.y
			
			# Setting Gaussian initial conditions in pressure

			if x <= 0.1:
				u = 0;
			else:
				u = 0;
			
			v = 0.0;
			p = np.exp(-np.log(2)*((x-0.5)**2 + (y-0.5)**2)/0.2**2);
			#p = 1.0
			rho = 1.0;

			cell.Q[0,0] = rho;
			cell.Q[0,1] = rho*u;
			cell.Q[0,2] = rho*v;
			cell.Q[0,3] = p/(0.4) + 1/2.*rho*(u**2 + v**2);

			#cell.vel.x = 1.0;
			#cell.vel.y = 1.0;
			cell.q[:] = cell.Q[0,:]
			pass;
		self.grid = grid;

		pass;

	def _write_vtk(self):
		import vtkgen as v
		vtk = v.vtkgen(self.grid,'_0');
		pass;


		pass;

if __name__ == '__main__':
	I = grid(sys.argv);
