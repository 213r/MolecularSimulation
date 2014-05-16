import Potential_MM_f2py  
import Potential_MM_f2py_backup  
from Constants import ang2bohr
import numpy as np

for i in np.arange(0.2, 3.0, 0.1):
    #print i,get_lj_potential_arar(i * ang2bohr)
    print i,Potential_MM_f2py.get_lj_force_arar(i * ang2bohr)
    print i,Potential_MM_f2py_backup.get_lj_force_arar(i * ang2bohr)
