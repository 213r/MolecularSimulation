from Molecule import Molecule  
def ReadRestartMC(fname): 
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
                return Molecule(atomnames = sym, positions = coord)  
            line = f.readline()

print ReadRestartMC("mc_restart.dat").get_positions_formated()  
