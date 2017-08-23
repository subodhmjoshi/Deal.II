# Common to triangle and gmsh
import grid as g
import numpy as np
import funcs as f
import copy


class reconstruct:
	''' 
		This class reconstructs a required order polynomial in the sv domain. Cell averaged variables are used for polynomial reconstruction.
		Inputs: 
			1. grid: object of class grid. (Refer module grid.py). 
				the sequence goes like,    Initialize - Reconstruct - Evolve. So grid object has arrived after initialization. 
			2. order: required order of accuracy. 1 will produce standard FV results
			3. options: to be decided later 

		Outputs: no explicit outputs. THe reconstruct class will have a method reconstruct(order) which will set left and right variable as well as flux vectors for all the faces. So that we can use approximate Riemann solver for flux computation. From another module 'solver.py' we will import object grid which is done with reconstruction. 

		References: 
			1. spectral volume method papers, series by z.j.wang
			2. Riemann solvers and numerical methods for fluid dynamics TORO, Leveque
	'''


	def __init__(self,grid,tl,order=1,option=0):
		self.grid = grid;
		self.order = order;
		self.option = option;

		for cell in self.grid.cells:
			cell.R[tl,:] = 0.0;
			pass;

		
		if self.order == 1:
			self._reconstruct_FO(); # First Order reconstruction done
		else:
			self._reconstruct_HO(); # Higher order reconstruction done
			pass;
		
		pass;

	def _reconstruct_FO(self):
		grid = self.grid;
		faces = grid.faces;
		cells = grid.cells;

		for face in faces:
			# Boundary treatment
			if face.marker == 0:
				# Internal cells
				face.ql = face.left.q[:];
				face.qr = face.right.q[:];

				face.Fl = f.find_flux(face.ql[:]);
				face.Fr = f.find_flux(face.qr[:]);
				
			elif face.marker == 1:
				# Boundary cells with outgoing boundary condition
				face.ql = face.left.q[:];
				face.qr = copy.deepcopy(face.ql);

				face.Fl = f.find_flux(face.ql[:]);
				face.Fr = f.find_flux(face.qr[:]);
				pass;

			elif face.marker == 2:
				# Solid wall boundary condition
				face.ql = face.left.q[:];
				
				ul = float(face.ql[1] / face.ql[0]);
				vl = float(face.ql[2] / face.ql[0]);

				vel_left = f.Vector(ul,vl);
				phi = f.findphi(face);
				
				''' velocity vector is rotated by angle phi, this gives unormal and utangent, then unormal is reversed and 
				everything is rotated back by angle phi. Refer notes for details
				'''
				vel_rot = f.rotate(vel_left,phi);
				vel_rot.x = -vel_rot.x # face-normal component reversed
				vel_new = f.rotback(vel_rot,phi);
				
				ur = vel_new.x;
				vr = vel_new.y;
				rhor = face.ql[0];
				er = face.ql[3];

				face.qr = np.array([rhor, rhor*ur, rhor*vr, er]);

				face.Fl = f.find_flux(face.ql[:]);
				face.Fr = f.find_flux(face.qr[:]);
				

			Uavg = f.ROE(face.ql, face.qr)
			L, R = f.eivect(Uavg, face.normal);
			
			eivaluel = f.eival(face.ql,face.normal);
			eivaluer = f.eival(face.qr,face.normal);
			eivaluem = f.eival(Uavg,face.normal);


			ql = face.ql[:];
			qr = face.qr[:];
			fl = f.directed_flux(face.Fl, face.normal);
			fr = f.directed_flux(face.Fr, face.normal);

			ql = ql.reshape(4,1);
			qr = qr.reshape(4,1);
			fl = fl.reshape(4,1);
			fr = fr.reshape(4,1);

			vl = np.dot(L,ql);
			vr = np.dot(L,qr);
			fl = np.dot(L,fl);
			fr = np.dot(L,fr);

			fchar = np.zeros(4,dtype='float');

			for k in range(4):
				alpha = 1.3 * np.max(np.abs([eivaluel[k], eivaluem[k], eivaluer[k]]));
				fchar[k] = 1/2. * (fl[k,0] + alpha * vl[k,0] + fr[k,0] - alpha * vr[k,0]);
				pass;

			fchar = fchar.reshape(4,1);

			flux = np.dot(R,fchar);

			for k in range(4):
				face.flux[k] = flux[k,0];
				pass;
			
		self.grid.faces = faces;

		pass;

	def _reconstruct_HO(self):
		# Higher-order reconstruction pending
		pass;
	pass;



				



























