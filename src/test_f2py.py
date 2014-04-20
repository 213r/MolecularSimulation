from Potential_MM_f2py import potential_mm_f2py, potential_mm_f2py_mc, potential_qmmm_f2py
from Molecule import Molecule 
import  numpy as np 
natom = 3
ndims = 3
atomnumbers = [18]*natom
#positions=np.array([[5.,0.,0.],[-5.,0.,0.],[0., 5., 0.],[0., -5., 0.]])  
positions=np.array([[105.,0.,0.],[-105.,0.,0.],[0., 105., 0.]])  
#positions = list(positions)

mol = Molecule(['C', 'Ar','H'],positions= [[0.,5.,0.],[5.,0.,0.],[-5.,0.,0.]]) 
positions = mol.get_positions() 
atomnumbers = mol.get_atomnumbers() 
ind = 0
rlimit = 30 

#print  potential_mm_f2py.__doc__ 
print  potential_qmmm_f2py.__doc__ 
sys.exit() 
print  potential_mm_f2py_mc.__doc__ 
energy = potential_mm_f2py_mc(ind, atomnumbers, positions.T, rlimit) 
print energy 
ind = 1
energy = potential_mm_f2py_mc(ind, atomnumbers, positions.T, rlimit) 
#energy, forces = potential_mm_f2py(atomnumbers, positions.T, rlimit) 
#energy, forces1, forces2 = potential_qmmm_f2py(atomnumbers1,atomnumbers2, positions1.T, positions2.T,rlimit) 
print energy 
