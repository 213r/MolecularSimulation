from Molecule import Molecule 
import numpy as np
from Constants import fs2tau, bohr2ang, ang2bohr
from VelocityVerlet import  VelocityVerlet

class Potential_TEST:

    def __init__(self,mol,pot_func,nrange = 1):
        self.mol = mol
        self.pot_func = pot_func 
        self.nrange = nrange

    def calc(self):
        self.mol.set_potential_energy_multi(self.get_potential_energy_multi())
        self.mol.set_potential_energy(self.get_potential_energy())
        self.mol.set_forces(-self.get_gradient())

    def get_potential_energy_multi(self):
        self.x = self.mol.get_positions()
        self.ene = self.pot_func.get_pot(self.x)
        return self.ene

    def get_gradient(self):
        return self.pot_func.get_grad(self.x)

    def get_potential_energy(self):
        return self.ene[0]

    def get_check_pbc(self): return False
    
    def get_pbc_adjusted(self, position):
        return position 

    def get_nrange(self): return self.nrange
         

class Test_Pot():
    def __init__(self):
        self.k = 0.0001

    def get_pot(self,xm):
        x = xm[0][0] 
        return np.array([0.5 * self.k * x * x])

    def get_grad(self,xm):
        x = xm[0] 
        return np.array([self.k * x])


x = Molecule(["X"],ndims=1)
x.set_positions([[0.0]])
x.set_velocities([[0.1]])

func = Test_Pot()
pot = Potential_TEST(x,func)
md = VelocityVerlet(x,pot,dt=0.1*fs2tau,nstep=100) 
md.run() 
