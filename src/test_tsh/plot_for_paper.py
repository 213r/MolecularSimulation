import  matplotlib.pyplot as plt
import re
import numpy as np
from Constants import hartree2kj, hartree2ev
from AnalyzeTrj import get_data 
from Molecule import Molecule
from math import pi
import sys

gen = get_data("md_tsh_qmmm.trj", start=0.0)
mol = Molecule(["C","H","O","O","H"]) 
t_ary = []
l_ary = []
l2_ary = []
l3_ary = []

lim = 1200 

for a in gen: 
    t, x, v = a
    if t > lim: break
    mol.set_positions(x,unit="bohr")
    #print t, mol.get_bond_length(1,2,unit="ang") 
    t_ary.append(t)
    l_ary.append(mol.get_bond_length(3,4,unit="ang"))
    l2_ary.append(mol.get_bond_length(0,1,unit="ang"))
    l3_ary.append(mol.get_bond_length(1,4,unit="ang"))
    #l_ary.append(mol.get_dihedral(1,0,2,4)*180/pi)

f, (ax1, ax2) = plt.subplots(2, sharex=True)

ax2.plot(t_ary,l_ary,label="O-H") 
ax2.plot(t_ary,l2_ary,label="C-H") 
ax2.plot(t_ary,l3_ary,label="H-H") 
#ax2.legend(loc="best")
ax2.legend()

def convert_num_ary(fname):
    ary = []  
    with open(fname) as f:
        line = f.readline()  
        if re.match("\s*#",line): line = f.readline() 
        while line: 
            ary.append(map(float, line.split()))
            line = f.readline() 
    return np.array(ary).T 

data = convert_num_ary("md_tsh_qmmm_energy.dat") 
#data = convert_num_ary("length_min_arh2.dat") 
t, ene = data[0], data[1:] 
ind = t < lim
t = t[ind] 
ene = ene[:,ind]
ene0 = ene[2][0]
toUnit = hartree2ev
for i in xrange(len(ene)): 
#    if i == 5: continue 
#    if 2 <= i <= 6: plt.plot(t, ene[i])
    if 2 <= i <= 4: ax1.plot(t, (ene[i]-ene0)*toUnit,label="S{}".format(i-2))
#    if i == 1: ax1.plot(t, ene[i]*toUnit,label="kinetic")
#    if  i == 7: plt.plot(t, ene[i])
ene6 = (ene[6]-ene0)*toUnit
freq = 400
ax1.plot(t[::freq], ene6[::freq],"k.",linewidth=1)
ax1.legend()
#plt.xlim(0,1200)
ax1.set_ylim(-1,16)
#ax1.set_yticks((0,2,4,6,8,10,12))
ax1.set_yticks((0,4,8,12,16))
ax2.set_xlabel("Time [fs]")
ax2.set_yticks((0,2,4,6,8,10))
ax1.set_ylabel("Energy [eV]")
ax2.set_ylabel("Distance [ang]")
f.subplots_adjust(hspace=0)
plt.show()
