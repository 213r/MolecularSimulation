import sys
import numpy as np
from math import sqrt

class LennardJonesPotential:

    def __init__(self,  a1,  a2): 
        if a1 < a2: self.atom1, self.atom2 = a1, a2
        else: self.atom1, self.atom2 = a2, a1 
        self.set_parameter_hcooh()

    def set_parameter_hcooh(self):
        # these parameters are deived from the optimization of the energy 
        # along the distance between Ar and each atom of HCOOH
        # the ab initio level is ccsd(t) and basis set is aug-cc-pVTZ  

        if self.atom1 == "ar" and self.atom2 == "ar": 
            # this LJ potential is  between two Ars
            self.sigma, self.eps = 6.899046647,0.000375861
        elif self.atom1 == "ar" and self.atom2 == "h1": 
            # this LJ potential is  between H and  Ar
            self.sigma, self.eps = 5.006593287, 0.000845957
            # the parameter for disociated H  
            # self.sigma, self.eps = 6.079175558, 0.00015622
        elif self.atom1 == "ar" and self.atom2 == "c2": 
            # this is the LJ potential between C2 and Ar
            self.sigma, self.eps = 5.927657824, 0.0010312
        elif self.atom1 == "ar" and self.atom2 == "o3": 
            # this is the LJ potential between O and Ar
            self.sigma, self.eps = 5.775431915, 0.000751319
        elif self.atom1 == "ar" and self.atom2 == "o4":
            # this is the LJ potential between O and Ar
            self.sigma, self.eps = 5.877809228, 0.000727761
        elif self.atom1 == "ar" and self.atom2 == "h5":
            # this is the LJ potential between O and Ar
            self.sigma, self.eps = 4.519668972, 0.001304569
        else: 
            print "no data of LJ potential between atom {} and {}"\
                .format(self.atom1, self.atom2)
            sys.exit()

    def get_energy(self,  r): 
        num6 = (self.sigma/ r)**6 
        return 4.0 * self.eps * num6 * (- 1.0 + num6)

    def get_force(self,  r):
        num6 = (self.sigma / r)**6 
        return - 24.0 * self.eps * num6 * (1.0 - 2 * num6) / r

def potential_mm(positions,atomnumbers, rlimit = 50.0):

    n = len(positions)
    dim = len(positions[0])
    energy = 0.0 
    force = np.zeros((n,dim), dtype=np.float64)

    for i in xrange(n-1):
        for j in xrange(i+1,n):
            vec_ji =  positions[i] - positions[j]
            lj_pot = LennardJonesPotential(atomnumbers[i],atomnumbers[j]) 
            rji = sqrt(np.sum(vec_ji * vec_ji)) 
            if rji > rlimit: next  
            energy += lj_pot.get_energy(rji) 
            force_ji =   vec_ji * (lj_pot.get_force(rji) / rji )
            force[j] -= force_ji 
            force[i] += force_ji
    return energy, force 

def potential_qmmm(positions1,atomnumbers1,positions2,atomnumbers2, \
                          rlimit = 50.0):
    n1 = len(positions1)
    n2 = len(positions2)
    dim = len(positions1[0])
    energy = 0.0 
    force1 = np.zeros((n1,dim), dtype=np.float64)
    force2 = np.zeros((n2,dim), dtype=np.float64)

    from Constants import bohr2ang

    for i in xrange(n1):
        for j in xrange(n2):
            vec_ji =  positions1[i] - positions2[j]
            lj_pot = LennardJonesPotential(atomnumbers1[i],atomnumbers2[j]) 
            rji = sqrt(np.sum(vec_ji * vec_ji)) 
            if rji > rlimit: next  
            energy += lj_pot.get_energy(rji) 
            print  i,j,rji*bohr2ang, lj_pot.get_energy(rji) 
            force_ji = vec_ji * (lj_pot.get_force(rji) / rji )
            force1[i] += force_ji
            force2[j] -= force_ji 
    return energy, force1, force2

