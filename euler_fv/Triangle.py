import sys
import os

#Refer Triangle documentation for understanding following commands:
#https://www.cs.cmu.edu/~quake/triangle.html

filename = sys.argv[1].split('.',1)[0];

if not os.path.exists(filename):
	os.makedirs(filename);

	pass;

os.chdir(filename);
os.system('cp ../%s .'%(sys.argv[1]));
os.system('triangle -pqcven %s'%(sys.argv[1]));
os.system('triangle -pqcvenra%s %s'%(sys.argv[2],filename+'.1'));
os.system('cd ..');
