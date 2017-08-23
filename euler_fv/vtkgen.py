# Common to triangle and gmsh
import numpy as np
import datetime
import os

class vtkgen:
	def __init__(self,grid,arg):
		self.grid = grid;
		self.arg = arg;
		cells = self.grid.cells;
		cell = cells[0];
		print '\nWriting vtk file... \n'

		self.directory = self.grid.name + '_%d_'%(self.grid.ncells)+'%d'%(self.grid.order+1)+'_'+ str(datetime.date.today())

		if not os.path.exists(self.directory): # Checks wether the directory of grid.name already exists or not
			os.makedirs(self.directory); # if not, creates a directory
		if self.grid.nsubcells == 1:
			self._write_vtk_FV();
		else:
			self._write_vtk_SV();
			pass;
		pass;

	def _write_vtk_FV(self):
		''' 
			Writes the vtk file for first order grid. i.e. Finite Voluem Method (FVM).
		For more information on how to write vtk files, please refer,
			1) Writing vtk using python:     http://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python
			2) vtk formats:     http://www.vtk.org/VTK/img/file-formats.pdf
			3) A very good example to get you started with vtk:
					http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
			4) vtk user's guide:      http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html

		'''
		    

		grid = self.grid;
		cells = grid.cells;
		faces = grid.faces;
		nodes = grid.nodes;
		ncells = grid.ncells;
		nnodes = grid.nnodes;
		nfaces = grid.nfaces;

		name = grid.name.split('/',1)[0]
		with open(self.directory+'/'+name + self.arg +'.vtk','w') as file:
			#Header
			print >> file, '# vtk DataFile Version 5.8'
			print >> file, 'Field data visualization for %s.vtk'%(grid.name);
			print >> file, 'ASCII'
			#Unstructure grid
			print >> file, 'DATASET UNSTRUCTURED_GRID'
			#Pointvalues
			print >> file, 'POINTS', nnodes, 'FLOAT'
			for node in nodes:
				print >> file, node.x, node.y, 0.0;
				pass;

			#Celldata
			size = 0
			for cell in cells:
				size += len(cell.nodes) + 1;
				pass;

			print >> file, 'CELLS', ncells, size;
			for cell in cells:
				file.write('%d '%(len(cell.nodes)))
				for node in cell.nodes:
					file.write('%d '%(node.index));
					pass;
				file.write('\n');
				pass;
			
			print >> file, 'CELL_TYPES', ncells;
			for cell in cells:
				print >> file, cell.Type;
				pass;

			#Pointdatasets
			'''
			As we are coding FV method, we do not require point (node) data. So that part is not coded in this module.
			Refer documentation to know about PointData. 
			it looks like,

			POINT_DATA <#nnodes>
			SCALARS some_name_like_temperature FLOAT
			LOOKUP_TABLE default
			value_for_point_0_int
			value_for_point_1_int
			.
			.
			.
			value_for_last_point_int
			
			SCALARS some_other_variable_like_pressure FLOAT
			LOOKUP_TABLE default
			value_for_point_0_int
			value_for_point_1_int
			.
			.
			.
			value_for_last_point_int
		
			VECTORS some_vector_data_like_velocity FLOAT
			LOOKUP_TABLE default
			i_value_for_point_0_int		j_value		k_value
			value_for_point_1_int		j_value		k_value

			.
			.
			.
			value_for_last_point_int	j_value		k_value
			'''
			#Celldatasets:
			

			print >> file, 'CELL_DATA', ncells;
			# Density
			print >> file, 'SCALARS rho FLOAT'
			print >> file, 'LOOKUP_TABLE default'
			'''
			for i in range(ncells):
				print >> file, np.sin(2*np.pi*i*0.01)
				pass;
			'''
			for cell in cells:
				print >> file, cell.q[0];
				pass;
			
			# Pressure
			print >> file, 'SCALARS P FLOAT'
			print >> file, 'LOOKUP_TABLE default'
			'''
			for i in range(ncells):
				print >> file, np.sin(2*np.pi*i*0.01)
				pass;
			'''
			for cell in cells:
				e = cell.q[3]
				rho = cell.q[0]
				u = cell.q[1]/rho
				v = cell.q[2]/rho
				P = (0.4) * (e - 1/2.*rho*(u**2 + v**2));
				print >> file, P;
				pass;

			# velocity
			print >> file, 'VECTORS velocity FLOAT'
			for cell in cells:
				u = cell.q[1]/cell.q[0];
				v = cell.q[2]/cell.q[0]
				print >> file, u, v, 0.0;
				pass;
			pass;
		pass;

	def _write_vtk_SV(self):
		''' 
			Writes the vtk file for first order grid. i.e. Finite Voluem Method (FVM).
		For more information on how to write vtk files, please refer,
			1) Writing vtk using python:     http://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python
			2) vtk formats:     http://www.vtk.org/VTK/img/file-formats.pdf
			3) A very good example to get you started with vtk:
					http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
			4) vtk user's guide:      http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html

		'''
		grid = self.grid;
		cells = grid.cells;
		faces = grid.faces;
		nodes = grid.nodes;
		ncells = grid.ncells;
		nnodes = grid.nnodes;
		nfaces = grid.nfaces;

		name = grid.name.split('/',1)[0]
		with open(self.directory+'/'+name +self.arg +'.vtk','w') as file:
			#Header
			print >> file, '# vtk DataFile Version 5.8'
			print >> file, 'Field data visualization for %s.vtk'%(grid.name);
			print >> file, 'ASCII'
			#Unstructure grid
			print >> file, 'DATASET UNSTRUCTURED_GRID'
			#Pointvalues
			print >> file, 'POINTS', ncells * len(cells[0].subcells)*len(cells[0].subcells[0].nodes), 'FLOAT'

			for cell in cells:
				for subcell in cell.subcells:
					for node in subcell.nodes:
						print >> file, node.x, node.y, 0.0;
						pass;
					pass;
				pass;

			#Celldata
			size = 0
			for cell in cells:
				for subcell in cell.subcells:
					size += len(subcell.nodes) + 1;
				pass;

			print >> file, 'CELLS', ncells*len(cells[0].subcells), size;
			i = 0
			for cell in cells:
				for subcell in cell.subcells:

					file.write('%d '%(len(subcell.nodes)))
					for node in subcell.nodes:
						file.write('%d '%(i));
						i += 1
						pass;
					file.write('\n');
					pass;
			
			print >> file, 'CELL_TYPES', ncells*len(cells[0].subcells);
			for cell in cells:
				for subcell in cell.subcells:
					print >> file, subcell.Type;
				pass;

			#Pointdatasets
			'''
			As we are coding FV method, we do not require point (node) data. So that part is not coded in this module.
			Refer documentation to know about PointData. 
			it looks like,

			POINT_DATA <#nnodes>
			SCALARS some_name_like_temperature FLOAT
			LOOKUP_TABLE default
			value_for_point_0_int
			value_for_point_1_int
			.
			.
			.
			value_for_last_point_int
			
			SCALARS some_other_variable_like_pressure FLOAT
			LOOKUP_TABLE default
			value_for_point_0_int
			value_for_point_1_int
			.
			.
			.
			value_for_last_point_int
		
			VECTORS some_vector_data_like_velocity FLOAT
			LOOKUP_TABLE default
			i_value_for_point_0_int		j_value		k_value
			value_for_point_1_int		j_value		k_value

			.
			.
			.
			value_for_last_point_int	j_value		k_value
			'''
			#Celldatasets:
			

			print >> file, 'CELL_DATA', ncells*len(cells[0].subcells);
			print >> file, 'SCALARS u FLOAT'
			print >> file, 'LOOKUP_TABLE default'
			'''
			for i in range(ncells):
				print >> file, np.sin(2*np.pi*i*0.01)
				pass;
			'''
			for cell in cells:
				for subcell in cell.subcells:
					print >> file, subcell.q;
				pass;

			pass;
		pass;
	pass;



