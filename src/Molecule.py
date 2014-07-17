import numpy as np
import sys
from Atom import Atom 
from math import sqrt, acos, pi
from Constants import bohr2ang, ang2bohr, fs2tau 

"""
    Molecule Module:
    
    -include-
    class: Molecule(self, atomnames) 
    method: bind_molecule(mol1, mol2):
     
"""

class Molecule:
    """ 
        Molecule class defines  
    
    """
    def __init__(self, atomnames = [], atomnumbers = [], positions = None, velocities = None,\
            forces = None, ndims = 3):
        assert atomnames != [] or atomnumbers != [], "define atomic symbols or atomic numbers"
        
        self.atomnames,self.atomnumbers = [], []
        self.masses = []
        self.ndims = ndims
        
        if atomnames != []:
            self.atomnames = [ atom.lower().capitalize() for atom in atomnames ]
            for atom in self.atomnames: 
                self.masses.append(Atom(symbol=atom).GetMass())
                self.atomnumbers.append(Atom(symbol=atom).GetAtomicNumber())

        if atomnumbers != []:
            self.atomnumbers = atomnumbers
            for natom in self.atomnumbers: 
                self.masses.append(Atom(Z=natom).GetMass())
                self.atomnames.append(Atom(Z=natom).GetChemicalSymbol())

        self.atomnames = np.array(self.atomnames)
        self.atomnumbers = np.array(self.atomnumbers)
        self.masses = np.array(self.masses) 
        self.natoms = len(self.atomnames) 

        if positions == None: self.positions = np.zeros((self.natoms,self.ndims))
        else:
            assert np.shape(positions) == (self.natoms,self.ndims), "Molecule class size error"
            self.positions = np.array(positions)
        
        if velocities == None: 
            self.velocities = np.zeros((self.natoms,self.ndims))
        else: 
            assert np.shape(velocities) == (self.natoms,self.ndims), "Molecule class size error" 
            self.velocities = np.array(velocities)

        if forces == None: 
            self.forces = np.zeros((self.natoms,self.ndims))
        else: 
            assert np.shape(forces) == (self.natoms,self.ndims), "Molecule class size error" 
            self.forces = np.array(forces)

    def __len__(self):
        return self.natoms

    def __iter__(self):
        for i in self.atomnames:
            yield i
    
    def get_shape(self):
        return (self.natoms,self.ndims)

    def get_ndims(self):
        return self.ndims
    
    def copy(self):
        return Molecule(atomnames = self.atomes,\
                   atomnumbers = self.atomnumbers,\
                   positions = self.positions,\
                   velocities = self.velocities,\
                   forces = self.forces,\
                   ndims = self.ndims)

    def set_addinfo(self,addinfo):
        self.addinfo = addinfo
    
    def get_addinfo(self):
        return self.addinfo 
    
    def get_atomnames(self):
        return self.atomnames 

    def get_atomnumbers(self):
        return self.atomnumbers

    def get_atomnumbers(self):
        return self.atomnumbers
    
    def get_positions(self):
        return self.positions
    
    def get_positions_index(self, i):
        return self.positions[i]
    
    def set_positions(self, positions, unit='bohr'): 
        p = np.array(positions)
        if self.natoms != len(p): 
            print "the length of the positions array is invalid!"
            sys.exit() 
        if unit == "ang": p *= ang2bohr 
        self.positions = p 

    def set_positions_index(self, i, position,  unit='bohr'): 
        x = np.array(position)  
        if unit == "ang": x *= ang2bohr 
        self.positions[i] = x 
        
    def get_velocities(self):
        return self.velocities
    
    def set_velocities(self, velocities, unit='bohr/tau'):
        v = np.array(velocities)
        if self.natoms != len(v): 
            print "the length of the velocities array is invalid!"
            sys.exit() 
        if unit == 'ang/fs': v *= bohr2ang / tau2fs 
        self.velocities = v 
    
    def set_velocities_index(self, i, velocities, unit='bohr/tau'):
        v = np.array(velocities)
        if unit == 'ang/fs': v *= bohr2ang / tau2fs 
        self.velocities[i] = v 
    
    def get_masses(self):
        return self.masses  
    
    def set_masses(self, masses):
        self.masses = np.array(masses)
    
    def get_momentum(self):
        return self.velocities * self.masses[:,np.newaxis]
    
    def set_momentum(self, mlist):
        self.velocities = np.array(mlist) / self.masses[:,np.newaxis]
    
    def get_forces(self):
        return self.forces 
  
    def set_forces(self, forces):
        if forces.shape != (self.natoms,self.ndims):
            print "the shape of forces is incorrect"
            sys.exit() 
        self.forces = np.array(forces) 
    
    def get_accelerations(self):
        return self.forces / self.masses[:, np.newaxis] 

    def get_potential_energy(self):
        return self.energy 
    
    def get_potential_energy_multi(self):
        return self.energy_multi 
    
    def set_potential_energy(self, energy):
        self.energy = energy

    def set_potential_energy_multi(self, energy_multi):
        self.energy_multi = energy_multi
    
    def get_kinetic_energy(self):
        return np.sum(self.masses[:, np.newaxis] * \
        self.velocities * self.velocities) * 0.5
    
    def get_kinetic_energy_per_atoms(self):
        return np.sum(self.masses[:, np.newaxis] * \
        self.velocities * self.velocities, axis=1) * 0.5

    def get_bond_length(self,i,j, unit="bohr"):
        a = self.positions[i] - self.positions[j]
        if unit == "ang": return sqrt(np.dot(a, a))*bohr2ang
        else: return sqrt(np.dot(a, a))

    def get_bond_angle(self,i,j,k):
        a = self.positions[i] - self.positions[j]
        b = self.positions[k] - self.positions[j]
        x = np.dot(a,b) / sqrt(np.dot(a, a) * np.dot(b, b)) 
        if abs(x) > 1.0:
            if x < 0.0: x += 1.e-15 
            else: x -= 1.e-15
        return acos(x)   

#    def get_dihedral(self,i,j,k,l):
#        vji = self.positions[i] - self.positions[j]
#        vjk = self.positions[k] - self.positions[j]
#        vkl = self.positions[l] - self.positions[k]
#        norm_vjk = np.dot(vjk, vjk)
#        vi_pro = vji - np.dot(vjk,vji) / norm_vjk * vjk   
#        vl_pro = vkl + np.dot(vjk,vkl) / norm_vjk * vjk   
#        x = np.dot(vi_pro,vl_pro) / sqrt(np.dot(vi_pro, vi_pro) * np.dot(vl_pro, vl_pro)) 
#        if abs(x) > 1.0:
#            if x < 0.0: x += 1.e-15 
#            else: x -= 1.e-15
#        return acos(x)   

    def get_dihedral(self,i,j,k,l):
        v_ji = self.positions[i] - self.positions[j]
        v_jk = self.positions[k] - self.positions[j]
        v_kl = self.positions[l] - self.positions[k]
        v_kj = -v_jk 
        v_ijk = np.cross(v_ji,v_jk) 
        v_jkl = np.cross(v_kj,v_kl) 
        #print v_ijk, v_jkl
        x = np.dot(v_ijk,v_jkl) / sqrt(np.dot(v_ijk, v_ijk) * np.dot(v_jkl, v_jkl)) 
        if abs(x) > 1.0:
            if x < 0.0: x += 1.e-15 
            else: x -= 1.e-15
        return acos(x)

    def get_positions_formated(self, unit='bohr', label = True, message = None):
        st = "    {}\n".format(self.natoms)
        if label: 
            if message == None: st += "   #Coordinates [{}] \n".format(unit) 
            else: st += "   #Coordinates [{}] {} \n".format(unit,message) 
        else: 
            if message == None: st += "\n"  
            else: st += "    {} \n".format(message) 
            
        for i in xrange(self.natoms):
            sym = self.atomnames[i]
            if unit == 'bohr': coord = self.positions[i]
            elif unit == 'ang': coord = self.positions[i]*bohr2ang
            else: 
                print 'Error: the unit of positions is wrong'
                sys.exit() 
            st += "{:s}   {:15.8f}   {:15.8f}   {:15.8f}\n".format(sym,*coord) 
        return st

    def get_velocities_formated(self, unit='bohr/tau', label = True, message = None):
        str = "    {}\n".format(self.natoms)
        if label: 
            if message == None: str += "   #Velocities [{}] \n".format(unit) 
            else: str += "   #Velocities [{}] {}\n".format(unit,message) 
        else: 
            if message == None: str += "\n"  
            else: str += "    {} \n".format(message) 
 
        for i in xrange(self.natoms):
            sym = self.atomnames[i]
            if unit == 'bohr/tau': vel = self.velocities[i]
            elif unit == 'ang/fs': vel = self.velocities[i] * bohr2ang / tau2fs
            else: 
                print 'Error: the unit of positions is wrong'
                sys.exit() 
            str += "{:s}   {:15.8f}   {:15.8f}   {:15.8f}\n".format(sym,*vel) 
        return str 

    def sort_by_positions(self, origin = np.zeros(3)):
        origin = np.array(origin)
        positions_ref = self.positions - origin[np.newaxis,:] 
        length2_ref = np.sum(positions_ref * positions_ref, axis=1) 
        index =  np.argsort(length2_ref)
        self.atomnames = self.atomnames[index] 
        self.masses = self.atomnumbers[index] 
        self.atomnumbers = self.atomnumbers[index] 
        self.positions = self.positions[index] 
        self.velocities = self.velocities[index] 
        self.forces = self.forces[index] 

    def substract_molecule(self, *which_atoms):
        # which_atoms assign the symbol of atoms which you want to substrace.  
        index = np.array([False] * self.natoms ,dtype=bool)
        for i, name in enumerate(self.atomnames):
            if name in which_atoms: index[i] = True 
        return Molecule(atomnames = self.atomnames[index],\
                positions = self.positions[index], \
                velocities = self.velocities[index],
                forces = self.forces[index], ndims = self.ndims) 

    def substract_molecule_index(self, index):
        index  = np.array(index) 
        return Molecule(atomnames = self.atomnames[index],\
                positions = self.positions[index], \
                velocities = self.velocities[index],
                forces = self.forces[index], ndims = self.ndims) 

    def read_coord_from_file(self, file, unit="bohr"):
        f = open(file)
        positions = np.array([map(float,i.split()) for i in f.readlines()])
        self.set_positions(positions, unit)

    def read_velocity_from_file(self, file, unit="bohr/tau"):
        f = open(file)
        velocities = np.array([map(float,i.split()) for i in f.readlines()])
        self.set_velocities(velocities, unit)

    def move2centermass(self):
        m = self.masses
        centermass = np.sum(m[:,np.newaxis] * self.positions,axis=0) / np.sum(m) 
        self.positions -= centermass
        return centermass 

def atom2mol(atom):
    return Molecule(atomnames=[atom.GetChemicalSymbol()],positions=[atom.GetCartesianPosition()], \
    velocities=[atom.GetCartesianVelocity()],forces=[atom.GetCartesianForce()]) 

def bind_molecule(mol1, mol2):
    if isinstance(mol1,Atom): mol1 = atom2mol(mol1) 
    if isinstance(mol2,Atom): mol2 = atom2mol(mol2) 
     
    atomnumbers = np.r_[mol1.get_atomnumbers(), mol2.get_atomnumbers()]
    positions = np.r_[mol1.get_positions(), mol2.get_positions()]
    velocities = np.r_[mol1.get_velocities(), mol2.get_velocities()]
    forces = np.r_[mol1.get_forces(), mol2.get_forces()]
    ndims  = mol1.get_ndims()  
    return Molecule(atomnumbers = atomnumbers, positions = positions,\
           velocities = velocities,forces = forces, ndims = ndims)


