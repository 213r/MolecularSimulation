import numpy as np 
from Molecule import Molecule 
from Constants import bohr2ang

def read_xyz(fname):
    f = open(fname, "r") 
    line = f.readline() 
    while line: 
        natom = int(line)
        doc = line  
        for _ in xrange(natom+1): doc += f.readline()   
        yield doc  
        line = f.readline() 

def get_positions(fname):
    f = open(fname, "r") 
    line = f.readline() 
    while line: 
        natom = int(line)
        positions = []  
        f.readline() 
        for _ in xrange(natom): positions.append(map(float,f.readline().split()[1:]))
        yield np.array(positions)  
        line = f.readline() 

if __name__ == '__main__':
    tot = ''  
    for i,doc in enumerate(read_xyz('md.xyz')):
        if i%10 == 0: tot += doc 
    print tot 
#
#mol = Molecule(["C", "H", "O", "O", "H"])
#for i, pos in enumerate(get_positions("md.xyz")):
#    mol.set_positions(pos, unit="ang")
#    print i, mol.get_bond_length(0,2)*bohr2ang  

