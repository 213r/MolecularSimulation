import numpy as np 
from copy import copy 
from math import exp, sqrt, pi
from WriteOutput import WriteOutputMC
from Constants import kelvin2hartree
from random import random, randint 
from MakeInitial import euler_rotation
import sys


"""
MonteCarlo (MC) simulation program 
"""

class MonteCarlo:
    
    def __init__(self, mol, pot, delta, nstep, temp, restart = False, frozen_atom_number = None ):
       
        self.mol, self.pot = mol, pot
        self.nstep = nstep 
        self.temp, self.restart = temp, restart
        self.interaction_save = np.zeros((len(self.mol), len(self.mol)))
        self.which_frozen = [False] * len(self.mol)  
        if not frozen_atom_number == None: 
            for i in frozen_atom_number: self.which_frozen[i] = True  

        if self.pot.get_check_pbc: 
            self.mol.set_positions([self.pot.get_pbc_adjusted(pos) for pos in self.mol.get_positions()])  

        if isinstance(delta, float): 
            self.delta = [delta] * len(self.mol) 
        elif isinstance(delta, list) or isinstance(delta, np.ndarray):
            if len(delta) == len(mol):
                self.delta = np.array(delta)
            else: 
                print 'check the array size of delta' 
                sys.exit() 
        else: 
            print  'delta type error'    
            sys.exit()  
        
        self.beta = 1.0 / (temp * kelvin2hartree)  
        self.count, self.count_accept = 0, 0  

        self.setup_output() 

    def access_writeoutput(self):
        return self.woutp

    def setup_output(self):
        self.woutp = WriteOutputMC()
        if self.restart: self.woutp.restart(self)
        else: self.woutp.start()

    def run(self):
        self.prepare() 
        # log at the initial point 
        if not self.restart: self.woutp.logging(self)  

        # loop of sampling  
        while self.count < self.nstep:
            self.count += 1 
            self.step() 
            self.woutp.logging(self)

        self.woutp.finalize(self)

    def prepare(self):
        for i in xrange(len(self.mol)):
            self.interaction_save[i,:] = self.pot.get_interaction_energy(i)
        self.mol.set_potential_energy(0.5 * np.sum(self.interaction_save)) 

    def step(self):
        # select the transfered atom  
        ind = randint(0,len(self.mol)-1)
        if self.which_frozen[ind]: return 
        x_save = copy(self.mol.get_positions_index(ind))  
        if self.pot.get_check_pbc(): 
            self.mol.set_positions_index(ind, self.pot.get_pbc_adjusted(x_save \
            + self.delta[ind] * random_v(self.mol.get_ndims())))
        else:
            self.mol.set_positions_index(ind, x_save \
            + self.delta[ind] * random_v(self.mol.get_ndims()))
        int_trial = self.pot.get_interaction_energy(ind) 
        delta_e = np.sum(int_trial) - np.sum(self.interaction_save[ind])  
        judge = exp(- self.beta * delta_e) > random()  
        # the trial movement is judged 
        if judge: 
            self.count_accept += 1 
            self.interaction_save[ind,:] = int_trial 
            self.interaction_save[:,ind] = int_trial 
            self.mol.set_potential_energy(self.mol.get_potential_energy() + delta_e)  
        else: self.mol.set_positions_index(ind, x_save)

class MonteCarlo_Ex(MonteCarlo):
    
    def __init__(self, mol, pot, delta, nstep, temp, restart = False, frozen_atom_number = None,\
            treated_as_molecule= None, delta_mol = None):
        
        MonteCarlo. __init__(self, mol, pot, delta, nstep, temp, restart, frozen_atom_number)
        self.treated_as_molecule = np.array(treated_as_molecule)
        if treated_as_molecule != None: 
            if delta_mol != None:
                self.alpha_save, self.beta_save, self.gamma_save = 0., 0., 0.   
                self.count_tot_mol, self.count_accept_mol = 0, 0 
                self.delta_mol = delta_mol  
                self.mol_save = self.mol.substract_molecule_index(self.treated_as_molecule)  
            else: 
                print 'provide the step size of targeted molecule'
                sys.exit() 
        self.count_tot = 0

    def step(self):
        # select the transfered atom  
        ind = randint(0,len(self.mol)-1)
        if ind in self.treated_as_molecule: 
            self.count_tot_mol += 1
            x_mol_save = copy(self.mol.get_positions_index(self.treated_as_molecule))  
            #print euler_rotation(*self.delta_mol) * x_mol_save.T
            alpha, beta, gamma = self.delta_mol[1] * random_v(3)
            alpha *= pi; beta *= pi * 0.5; gamma *= pi 
            x_mol_try = (euler_rotation(alpha, beta, gamma) * x_mol_save.T).T
            x_mol_try += self.delta_mol[0]*random_v(self.mol.get_ndims()) 
            self.mol.set_positions_index(self.treated_as_molecule, x_mol_try) 
            int_trial_matrix = np.array([self.pot.get_interaction_energy(i) for i in self.treated_as_molecule])
            int_save_matrix = np.array([self.interaction_save[i] for i in self.treated_as_molecule])
            delta_e = np.sum(int_trial_matrix) - np.sum(int_save_matrix)  
            judge = exp(- self.beta * delta_e) > random()  
            # the trial movement is judged 
            if judge: 
                self.count_accept_mol += 1 
                for i in self.treated_as_molecule: self.interaction_save[i,:] = int_trial_matrix[i]  
                for i in self.treated_as_molecule: self.interaction_save[:,i] = int_trial_matrix[i]  
                self.mol.set_potential_energy(self.mol.get_potential_energy() + delta_e)  
                self.alpha_save += alpha; self.beta_save += beta; self.gamma_save += gamma 
            else: self.mol.set_positions_index(self.treated_as_molecule, x_mol_save)
        else:
            self.count_tot += 1
            if self.which_frozen[ind]: return 
            x_save = copy(self.mol.get_positions_index(ind))  
            if self.pot.get_check_pbc(): 
                self.mol.set_positions_index(ind, self.pot.get_pbc_adjusted(x_save \
                + self.delta[ind] * random_v(self.mol.get_ndims())))
            else:
                self.mol.set_positions_index(ind, x_save \
                + self.delta[ind] * random_v(self.mol.get_ndims()))
            int_trial = self.pot.get_interaction_energy(ind) 
            delta_e = np.sum(int_trial) - np.sum(self.interaction_save[ind])  
            judge = exp(- self.beta * delta_e) > random()  
            # the trial movement is judged 
            if judge: 
                self.count_accept += 1 
                self.interaction_save[ind,:] = int_trial 
                self.interaction_save[:,ind] = int_trial 
                self.mol.set_potential_energy(self.mol.get_potential_energy() + delta_e)  
            else: self.mol.set_positions_index(ind, x_save)


def random_v(size):
    rand = 2.0 * np.random.random(size) - 1.0
    return rand / sqrt(np.sum(rand**2)) 

class PeriodicBoundaryCondition:

    def __init__(vlength): 
        self.box = np.diag(np.ones(3) * vlength)
        self.invbox = np.linalg.inv(self.box) 

    def periodic_boundary_adjustment(vlength, r):
        s = np.dot(self.invbox, r)  
        return np.dot(self.box, s - np.round(s))  

if __name__ == '__main__':
    print np.sum(random_v(3)**2) 

