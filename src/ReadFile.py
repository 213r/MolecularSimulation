
def read_xyz(fname, freq = 1):
    with open(fname) as f:
        line = f.readline()
        sym, coord = [], []  
        while line:
            if "#Coordinates [" in line:
                aline = line.split() 
                natom = int(aline[0]) 
                unit = aline[2][1:-1]
                line = f.readline()
                for i in xrange(natom): 
                    aline = f.readline().split() 
                    sym.append(aline[0])
                    coord.append(map(float, aline[1:]))
                mol = Molecule(atomnames = sym)
                mol.set_positions(coord,unit=unit) 
                return mol 
            line = f.readline()

def read_positions_formated(fname): 
    with open(fname) as f:
        line = f.readline()
        sym, coord = [], []  
        while line:
            if "#Coordinates [" in line:
                aline = line.split() 
                natom = int(aline[0]) 
                unit = aline[2][1:-1]
                line = f.readline()
                for i in xrange(natom): 
                    aline = f.readline().split() 
                    sym.append(aline[0])
                    coord.append(map(float, aline[1:]))
                mol = Molecule(atomnames = sym)
                mol.set_positions(coord,unit=unit) 
                return mol 
            line = f.readline()

def read_velocities_formated(fname): 
    with open(fname) as f:
        line = f.readline()
        sym, coord = [], []  
        while line:
            if "#Velocities [" in line:
                aline = line.split() 
                natom = int(aline[0]) 
                unit = aline[2][1:-1]
                line = f.readline()
                for i in xrange(natom): 
                    aline = f.readline().split() 
                    sym.append(aline[0])
                    coord.append(map(float, aline[1:]))
                mol = Molecule(atomnames = sym)
                mol.set_velocities(coord,unit=unit) 
                return mol 
            line = f.readline()


