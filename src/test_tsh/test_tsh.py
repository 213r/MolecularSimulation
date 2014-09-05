from Molecule import Molecule 
import numpy as np
from Constants import fs2tau, bohr2ang, ang2bohr
from VelocityVerlet import  VelocityVerlet
import sys
from math import exp, sqrt
import matplotlib.pyplot as plt
from TullySurfaceHopping import TullySurfaceHopping 

class Potential_TEST:

    def __init__(self,mol,pot_func,now_state = 0, nrange = 1):
        self.mol = mol
        self.pot_func = pot_func 
        self.nrange = nrange
        self.now_state = now_state

    def set_now_state(self,now_state):
        self.now_state = now_state
    
    def calc(self):
        self.mol.set_potential_energy_multi(self.get_potential_energy_multi())
        self.mol.set_potential_energy(self.get_potential_energy())
        self.mol.set_forces(-self.get_gradient())
        self.x = self.mol.get_positions()
        self.d_multi = self.pot_func.get_nacme(self.x) 

    def get_potential_energy_multi(self):
        self.x = self.mol.get_positions()
        self.ene = self.pot_func.get_pot(self.x)
        return self.ene

    def get_gradient(self):
        return self.pot_func.get_grad(self.x)[self.now_state]

    def get_potential_energy(self):
        return self.ene[self.now_state]

    def get_check_pbc(self): return False
    
    def get_pbc_adjusted(self, position):
        return position 

    def get_nrange(self): return self.nrange
    
    def get_now_state(self): return self.now_state

    def get_inner_v_d(self):
        #return self.mol.get_velocities()*self.d_multi
        return np.einsum("pq,ijpq -> ij",self.mol.get_velocities(), self.d_multi)

    def get_velocity_adjustment_vecotr(self): 
        return self.d_multi  

    def check_nacme_sign(self):
        pass
    
    def get_nacme_multi(self):
        return self.d_multi 

class Test_Pot():
    def __init__(self):
        self.k = 0.0001

    def get_pot(self,xm):
        x = xm[0][0] 
        return np.array([0.5 * self.k * x * x])

    def get_grad(self,xm):
        x = xm[0] 
        return np.array([self.k * x])

#def nac(x):
#    tan2 = 2.0 * v12(x) / (v11(x) - v22(x))  
#    dadr = 2.0 * ( dv12(x) * (v11(x) - v22(x)) -  v12(x) * \
#            (dv11(x) - dv22(x))) / (v11(x) -v22(x))**2
#    return -dadr / (2 * (1 + tan2 * tan2))

class Test_Pot_Tully():

    def __init__(self):
        self.a = 0.01
        self.b = 1.6
        self.c = 0.005
        self.d = 1.0

    def get_pot(self,xm):
        x = xm[0][0]
        v11 = self.v11(x)
        v22 = -v11
        v12 = self.v12(x)
        v21 = v12
        a = 0.5*(v11+v22)
        b = 0.5*sqrt((v11 + v22)**2 - 4.0*(v11*v22-v12*v21))
        return np.array([a-b,a+b])

    def get_grad(self,xm):
        x = xm[0][0]
        v11 = self.v11(x)
        v22 = -v11
        v12 = self.v12(x)
        v21 = v12
        dv11 = self.dv11(x)
        dv22 = -dv11
        dv12 = self.dv12(x)
        dv21 = dv12
        a = 0.5*(dv11+dv22)
        b = 0.5*((v11 + v22)*(dv11+dv22) - 2.0 * (dv11*v22 + v11*dv22 \
          - dv12*v21 - v12*dv21)) / sqrt((v11 + v22)**2 - 4.0*(v11*v22-v12*v21))
        return np.array([[[a-b,0,0]],[[a+b,0,0]]])

#    def diff_adiabatic(x):
#        a = 0.5*(dv11(x) + dv22(x))
#        b = 0.5*((v11(x) - v22(x)) * (dv11(x) - dv22(x)) + \
#            4 * v12(x) * dv12(x)) / \
#            np.sqrt((v11(x) - v22(x))**2 + 4*v12(x)**2)
#        return a - b, a + b

    def get_nacme(self,xm):
        x = xm[0][0]
        v11 = self.v11(x)
        v22 = -v11
        v12 = self.v12(x)
        v21 = v12
        dv11 = self.dv11(x)
        dv22 = -dv11
        dv12 = self.dv12(x)
        dv21 = dv12
        tan2 = 2.0 * v12 / (v11 - v22)  
        dadr = 2.0 * (dv12 * (v11 - v22) -  v12 * \
                     (dv11 - dv22)) / (v11 -v22)**2
        x = -dadr / (2 * (1 + tan2 * tan2))
        z = np.zeros((2,2,1,3))
        z[0,1] = np.array([x,0.,0.]) 
        z[1,0] = -z[0,1]
        return z 

    def v11(self, x):
        if x <= 0:
            return -self.a*(1.0 - exp(self.b*x))
        else:
            return self.a*(1.0 - exp(-self.b*x))

    def dv11(self, x):
        if x <= 0:
            return self.a * self.b * exp(self.b*x)
        else:
            return self.a * self.b * exp(-self.b*x)

    def v12(self, x):
        return self.c * exp(-self.d*x*x)

    def dv12(self,x):
        return -2 * x * self.c * self.d * exp(-self.d*x*x)


#pot = Test_Pot_Tully()
#x = np.linspace(-10,10,100)
#p_ary = np.array([ pot.get_nacme([[i]]) for i in x]) / 50
#plt.plot(x,p_ary[:])
#p_ary = np.array([ pot.v12(i) for i in x])
#plt.plot(x,p_ary[:])
#p_ary = np.array([ pot.get_pot([[i]]) for i in x])
#plt.plot(x,p_ary[:,0])
#plt.plot(x,p_ary[:,1])
#plt.show()
#sys.exit()
#

x = Molecule(["X"])
x.set_masses([2000])
x.set_positions([[-5.0,0,0]])
x.set_velocities([[0.0180,0,0]])
#print x.get_masses()
#print x.get_kinetic_energy()
print x.get_momentum()
#sys.exit()
func = Test_Pot_Tully()
pot = Potential_TEST(x,func,0,2)
#md = VelocityVerlet(x,pot,dt=0.2*fs2tau,nstep=5000) 
md = TullySurfaceHopping(x,pot,dt=0.2*fs2tau,nstep=1500) 
md.run() 
