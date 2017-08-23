# Common to gmsh and triangle except the last part < if __name__ == '__main__' > 
import os
import reconstruction as r
import sys
import grid as g
import funcs as f
import numpy as np
import vtkgen as v
from solver import solve

input_file ='test.poly';
time = 2;
CFL = 0.8;
TI = 'FE'
outfiles = 150 # no. of output files required
minimumarea = 0.0007

# Refer Triangle.py for understanding how the triangle files are created from poly file.
os.system('python Triangle.py %s %f'%(input_file,minimumarea));
			

#d = grid(option);
name = input_file.split('.',1)[0];
# Arguments list is manually created. 
option = [name+'/'+name+'.2.node',name+'/'+name+'.2.v.node',name+'/'+name+'.2.edge',name+'/'+name+'.2.v.edge',name+'/'+name+'.2.neigh',name+'/'+name+'.2.ele']
# Refer grid.py 
I = g.Initialize(option);
# for gmsh, instead of this, do      I = g.Initialize(sys.argv[1]);
# Plug in your solver here
solve(I.grid,time,CFL,TI,outfiles);


