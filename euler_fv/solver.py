# Common to gmsh and triangle except the last part < if __name__ == '__main__' > 

import reconstruction as r
import sys
import grid as g
import funcs as f
import numpy as np
import vtkgen as v

def solve(grid,time,CFL,TI,op):
	if TI == 'FE':
		_forward_euler(grid,time,CFL,op)
	elif TI == 'RK2':
		_RK2(grid,time,CFL,op);
	elif TI == 'RK3':
		_RK3(grid,time,CFL,op);

def _forward_euler(grid,time,CFL,op):
	t = 0
	i = 1;
	while t <= time:
		for cell in grid.cells:
			cell.q[:] = cell.Q[0,:];
			pass;

		R = r.reconstruct(grid,0);
		#flux(grid);
		dt = compute_dt(grid,CFL);

		for face in grid.faces:
			if face.marker == 0:
				face.right.R[0,:] += face.flux * face.length;
				face.left.R[0,:] -= face.flux * face.length;
			else:
				face.left.R[0,:] -= face.flux * face.length;
				pass;
		
		for cell in grid.cells:
			cell.R[0,:] = cell.R[0,:]/ cell.area;

			cell.Q[1,:] = cell.Q[0,:] + dt*cell.R[0,:];
			cell.Q[0,:] = cell.Q[1,:];
			cell.q[:] = cell.Q[0,:];
			pass;

		print t/time*100, dt;
		#print '_____________________________________',i,t/time*100,t/time*100-i*100/5.0,50*dt/time;
		grid = R.grid;
		
		if np.abs(t/time*100 - i*100/float(op)) <= 50*dt/time:
			#print '_____________________________________',t/time*100,t/time*100-i*100/5.0,50*dt/time;
			v.vtkgen(grid,'_%d'%(int(i*100/float(op))));
			i += 1;
			pass;
		t += dt;
		pass;

	v.vtkgen(grid,'_done');
	pass;

def _RK2(grid,time,CFL,op):
	t = 0;
	i = 1;
	
	while t <= time:

		for tl in range(2):

			for cell in grid.cells:
				cell.q[:] = cell.Q[tl,:];
				pass;

			R = r.reconstruct(grid,tl);
			#flux(grid); # flux is computed inside reconstruction module for eulers equations
			dt = compute_dt(grid,CFL);

			for face in grid.faces:
				if face.marker == 0:
					face.right.R[tl,:] += face.flux * face.length;
					face.left.R[tl,:] -= face.flux * face.length;
				else:
					face.left.R[tl,:] -= face.flux * face.length;
					pass;
			
			for cell in grid.cells:

				cell.R[tl,:] = cell.R[tl,:]/ cell.area;
				'''
				if tl == 0:
					cell.Q[1,:] = cell.Q[0,:] + dt*cell.R[0,:];
				else:
					cell.Q[2,:] = 0.5*(cell.Q[0,:] + cell.Q[1,:] + dt*cell.R[1,:]);
					pass;
				'''
				cell.Q[1,:] = cell.Q[0,:] + dt*cell.R[0,:];
				cell.Q[2,:] = 0.5*(cell.Q[0,:] + cell.Q[1,:] + dt*cell.R[1,:]);
				pass;
			pass;
		
		for cell in grid.cells:
			cell.Q[0,:] = cell.Q[2,:];
			cell.q[:] = cell.Q[0,:];
			pass;

		print t/time*100;
		#print '_____________________________________',i,t/time*100,t/time*100-i*100/5.0,50*dt/time;
		grid = R.grid;
		
		if np.abs(t/time*100 - i*100/float(op)) <= 50*dt/time:
			#print '_____________________________________',t/time*100,t/time*100-i*100/5.0,50*dt/time;
			v.vtkgen(grid,'_%d'%(int(i*100/float(op))));
			i += 1;
			pass;
		t += dt;
		pass;

	v.vtkgen(grid,'_done');
	pass;


def _RK3(grid,time,CFL,op):
	t = 0;
	i = 1;
	
	while t <= time:
		for tl in range(3):
			for cell in grid.cells:
				cell.q[:] = cell.Q[tl,:];
				pass;
			R = r.reconstruct(grid,tl);
			#flux(grid);
			dt = compute_dt(grid,CFL);

			for face in grid.faces:
				if face.marker == 0:
					face.right.R[tl,:] += face.flux * face.length;
					face.left.R[tl,:] -= face.flux * face.length;
				else:
					face.left.R[tl,:] -= face.flux * face.length;
					pass;
			
			for cell in grid.cells:
				cell.R[tl,:] = cell.R[tl,:]/ cell.area;
				'''
				if tl == 0:
					cell.Q[1,:] = cell.Q[0,:] + dt*cell.R[0,:];
					pass;
				elif tl == 1:
					cell.Q[2,:] = 3./4. * cell.Q[0,:] + 1/4.*cell.Q[1,:] + 1/4.*dt*cell.R[1,:];
					pass;

				else:
					cell.Q[3,:] = 1/3.*cell.Q[0,:] + 2/3.*cell.Q[2,:] + 2/3.*dt*cell.R[tl,:];
					pass;
				'''
				cell.Q[1,:] = cell.Q[0,:] + dt*cell.R[0,:];
				cell.Q[2,:] = 3./4. * cell.Q[0,:] + 1/4.*cell.Q[1,:] + 1/4.*dt*cell.R[1,:];
				cell.Q[3,:] = 1/3.*cell.Q[0,:] + 2/3.*cell.Q[2,:] + 2/3.*dt*cell.R[2,:];
				
				pass;
			pass;
		
		for cell in grid.cells:
			cell.Q[0,:] = cell.Q[3,:];
			cell.q[:] = cell.Q[0,:];
			pass;

		print t/time*100;
		#print '_____________________________________',i,t/time*100,t/time*100-i*100/5.0,50*dt/time;
		grid = R.grid;
		
		if np.abs(t/time*100 - i*100/float(op)) <= 50*dt/time:
			#print '_____________________________________',t/time*100,t/time*100-i*100/5.0,50*dt/time;
			v.vtkgen(grid,'_%d'%(int(i*100/float(op))));
			i += 1;
			pass;
		t += dt;
		pass;

	v.vtkgen(grid,'_done');
	pass;

'''
Following flux function is never called from here. It is computed in reconstruction module
'''
def flux(grid):
	for face in grid.faces:
		face.flux = 1/2.*(face.Fl + face.Fr - np.abs(face.un)*(face.qr - face.ql));
		pass;
	pass;

def compute_dt(grid, CFL):
	dt = 1000;

	for cell in grid.cells:
		denom = 0;
		for face in cell.faces:

			U = f.ROE(face.ql, face.qr);

			rho = float(U[0]);
			u   = U[1] / rho;
			v   = U[2] / rho;
			e   = U[3];
			p   = (1.4-1) * (e - 1/2.*rho *(u**2+v**2));
			h   = (e + p) / rho;
			a   = np.sqrt((1.4 - 1) * (h - 1/2.*(u**2+v**2)));

			un = u*face.normal.x + v*face.normal.y;	
			denom += face.length*(un+a);
			pass;
		tstep = CFL * cell.area / denom;
		if dt >= np.abs(tstep):
			dt = np.abs(tstep);
			pass;
		pass;
	return dt;

if __name__ == '__main__':

	option = []
	k = 0
	for arg in sys.argv:
		k+=1
		option.append(arg)
		pass;
	#d = grid(option);
	I = g.Initialize(option);
	# for gmsh, instead of this, do      I = g.Initialize(sys.argv[1]);
	solve(I.grid,0.1,0.8,'FE',25);

