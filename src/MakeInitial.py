from Constants import Kboltz, ang2bohr
import numpy as np 
from Atom import Atom 
from Molecule import Molecule
from math import ceil, sqrt, cos, sin 
import copy
"""
    MakeInitial Module 
    provides the useful tool for setting up the initial coordinates or velocities.  
     
    - include - 
     : SetLattice(self, atomname, estimate_num, dens) 
     : SetMaxwell(self,mol, temp)
     : superpositioned_atoms_delete(mol_deleted, mol_super)
     : euler_rotation(alpha, beta, gamma) 

"""


class SetLattice:
    """
        SetLattice provides the coordinates of lattice structure. 
    """

    def __init__(self, atomname, estimate_num, dens):
        self.atomname = atomname 
        self.n = estimate_num 
        mass = Atom(symbol=self.atomname).GetMass() 
        self.numdens = dens/9.10939e-4/mass/ang2bohr**3  
        self.vlength = (self.n/self.numdens)**(1.0/3.0)
        self.make_lattice() 

    def get_molecule(self):
        n = len(self.positions)
        return Molecule([self.atomname] * n, positions = \
                self.positions) 
    
    def get_positions(self):
        return self.positions  

    def get_positions_formated(self):
        str = ""
        for coord in self.positions: 
            str += "{:15.8f}   {:15.8f}   {:15.8f}\n".format(*coord) 
        return str 

    def get_lattice_length(self):
        return self.vlength

    def make_lattice(self):
        iq = int(ceil((self.n/4.0)**(1.0/3)))
        vlength_unit = (4.0/self.numdens)**(1.0/3.0) 
        base = np.array([[0.5,0.5,0],[0.0,0.5,0.5],\
                [0.5,0.0,0.5],[0.0,0.0,0.0]]) 
        ary = np.copy(base)   
        for i in xrange(iq):
            for j in xrange(iq):
                for k in xrange(iq):
                    if i == 0 and j == 0 and k == 0: continue 
                    x = np.array([i,j,k])
                    ary = np.vstack((ary, base + x)) 
        self.positions = ary * vlength_unit - 0.5 * self.vlength * np.ones(3) 

class SetMaxwell:
    """
        SetMaxwell produces the velocities from the Maxwell Boltzman distribution.  
    """

    def __init__(self,mol, temp):
        self.mol = mol
        self.temp = temp
    
    def set_velocities(self):
        self.make_velocities()
        self.mol.set_velocities(self.vel)

    def get_velocities(self):
        self.make_velocities()
        return self.vel
    
    def get_velocities_formated(self):
        self.make_velocities()
        str = ""
        for vel in self.vel: 
            str += "{:15.8f}   {:15.8f}   {:15.8f}\n".format(*vel) 
        return str 

    def make_velocities(self):
        vel_mw = sqrt(Kboltz * self.temp) * np.random.randn(*self.mol.get_shape())
        self.vel = vel_mw / np.sqrt(self.mol.get_masses()[:,np.newaxis])

def superpositioned_atoms_delete(mol_deleted, mol_super):
    """
        When two atoms are too close, the one of them is deleted.
    """
    index = [True] * len(mol_deleted)
    for i, ni in enumerate(mol_deleted.get_atomnumbers()):
        for j, nj in enumerate(mol_super.get_atomnumbers()):
            # Currently, the specific atoms are used only.  
            if nj == 1 and ni == 18: rlim = 3.62896 * 0.9  
            if nj == 6 and ni == 18: rlim = 3.70352 * 0.9 
            if nj == 8 and ni == 18: rlim = 3.83503 * 0.9 
            vec_ij = mol_deleted.get_positions()[i] - mol_super.get_positions()[j] 
            if np.sum(vec_ij * vec_ij) < rlim * rlim:
                index[i] = False
                break 
    return mol_deleted.substract_molecule_index(index)   

def euler_rotation(alpha, beta, gamma):
    """
    xyz type euler angle
    alpha, gamma : [-pi:pi] 
    beta : [-pi/2: pi/2] 
    R_rot = RX(alpha) * RY(beta) * RZ(gamma) 
    """

    def matrix_rot_x(r):
        return np.matrix([
            [1., 0., 0.],
            [0., cos(r), sin(r)],
            [0.,-sin(r), cos(r)]
        ]) 
    
    def matrix_rot_y(r):
        return np.matrix([
            [cos(r), 0., -sin(r)],
            [0., 1., 0.],
            [sin(r), 0.,  cos(r)]
        ]) 

    def matrix_rot_z(r):
        return np.matrix([
            [cos(r), sin(r), 0.],
            [-sin(r), cos(r), 0.],
            [0., 0., 1.]
        ]) 

    return matrix_rot_x(alpha) * matrix_rot_y(beta) \
            * matrix_rot_z(gamma) 

