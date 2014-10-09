#from scipy import linalg
from numpy import linalg
from random import random  
from Constants import fs2tau, tau2fs 
from math import sqrt 
import os, sys
import numpy as np
from MolecularDynamics import MolecularDynamics
from WriteOutput import WriteOutputMD_TSH, WriteOutputMD_TSH_QMMM 

class TullySurfaceHopping(MolecularDynamics):
    """ 
    Tully Surface Hopping method 
        is for dealng with non adiabatic transition process during classical MD simulation. 

         J. C. Tully, J. Chem. Phys. 93, 1061 (1990).
    
    """
    
    def __init__(self, mol, pot, dt, nstep, tsh_times = 5, restart = False,\
            tsh_ediff_thresh=0.1, tsh_ediff_factor=0.25, check_mdstop_dispersion = False,\
            limit_dispersion = 10.0, tlim = float('inf')):
      
        MolecularDynamics.__init__(self, mol, pot, dt, nstep, restart,\
                check_mdstop_dispersion, limit_dispersion, tlim)

        self.now_state = self.pot.get_now_state()
        self.nrange = self.pot.get_nrange()
        self.m_root = np.sqrt(self.mol.get_masses()) 

        self.tsh_times = tsh_times  
        self.delta_st = 0.5 * self.dt / self.tsh_times 
        self.dt_bak = self.dt 
    
        self.tsh_ediff_thresh = tsh_ediff_thresh
        self.tsh_ediff_factor = tsh_ediff_factor

        self.c = np.zeros(self.nrange)
        self.c[self.now_state] = 1.0 
   
        self.count = 0
        self.count_tsh = 0
        self.setup_output() 

    def setup_output(self):
        self.woutp = WriteOutputMD_TSH()
        if self.restart: self.woutp.restart(self)
        else: self.woutp.start(self)

    def run(self):
        self.prepare()
        self.woutp.logging(self)
        while self.count < self.nstep:
            
            self.count += 1

            self.judge_tsh()

            self.step()
            self.pre_judge_tsh()
            self.check_ediff()
            self.judge_tsh()
         
            self.count_tsh = 0
            self.woutp.logging(self)
        
            if self.mdstop_time_trans(): break
            if self.check_mdstop_dispersion and self.mdstop_atom_dispersion(): break

        self.woutp.finalize(self)

    def set_coefficients(self, coefficients):
        self.c = np.array(coefficients)

    def prepare(self):
        #self.pot.calc(self.count)
        self.pot.calc()
        self.pre_judge_tsh()
        self.check_ediff()
    
    def step(self):
        dt = self.dt
        get_accelerations_save =  self.mol.get_accelerations()
        self.mol.set_positions(self.mol.get_positions() + \
                self.mol.get_velocities() *  dt + 0.5 * \
                get_accelerations_save * dt * dt)
        #self.pot.calc(self.count) 
        self.pot.calc() 
        self.mol.set_velocities(self.mol.get_velocities() + \
                0.5 * (self.mol.get_accelerations() + \
                get_accelerations_save) * dt)
        self.elaptime += dt 

    def pre_judge_tsh(self):
        self.pot_multi = self.pot.get_potential_energy_multi()
        #print self.pot.get_inner_v_d()
        #self.nac =  np.diag(self.pot_multi) -  1.j*self.vd
        self.nac =  np.diag(self.pot_multi) -  1.j*self.pot.get_inner_v_d() 
        self.eig, self.trans = linalg.eigh(self.nac)
        self.woutp.write_transition_prob("Time: {}, Count: {}, Current: {}\n".format(self.elaptime * tau2fs, self.count, self.now_state))

    def judge_tsh(self):
        #print self.c
        while self.count_tsh < self.tsh_times:

            # solving eq.8 in the reference papaer 
            # for getting the expansion coordinate
            # exponential expansion 
            c = np.dot(np.conjugate(self.trans).T, self.c)
            c = np.dot(np.diag(np.exp(-1j*self.eig*self.delta_st)),c)
            self.c = np.dot(self.trans,c)
            print self.count_tsh, self.c 
            #print self.delta_st
            # calculation of the swithing probability 
            #   bkl  =-2Re(akl*VDkl)
            #   gkj  = deltat*bjk/akk
            a = np.outer(self.c,np.conjugate(self.c))
            b = -np.real(np.conjugate(a)*self.nac*2.j)
            #b = -2.0 * np.real(np.conjugate(a)*self.vd)
            #print b
            #print b2
            # judgement  
            rand = random()
            tot = 0.0
            self.woutp.write_transition_prob("step{}: ".format(self.count_tsh))
            for i in xrange(self.nrange): self.woutp.write_transition_prob("a{}:{:5.3f} ".format(i,a[i,i].real)) 
            self.woutp.write_transition_prob("| rand:{:8.5f}\n       Prob ".format(rand))
            
            for i in xrange(self.nrange):
                if i < self.now_state - 1 or i == self.now_state or i > self.now_state + 1: continue 
                g = self.delta_st*b[i,self.now_state] / \
                    a[self.now_state,self.now_state].real
                # log the hopping probability 
                self.woutp.write_transition_prob(" {}->{}:{:5.3f}".format(self.now_state,i,g))
                if tot < rand < tot + g:
                    self.woutp.write_transition_prob("\nHopping Starts\n")
                    self.woutp.write_trj_message("Hopping Starts [Probability Ocuuracne] \n")
                    if self.velocity_adjustment(i): return
                tot += g
            self.count_tsh += 1
            self.woutp.write_transition_prob("\n")

    def velocity_adjustment(self,cand_state):
            
        # unit vector alog the mass-weighted coordinate
        d = self.pot.get_velocity_adjustment_vecotr()[self.now_state, cand_state]
        d = d / self.m_root[:,np.newaxis] 
        d = d / linalg.norm(d) 
       
        velo = self.mol.get_velocities()

        in_prd = np.sum(self.m_root[:,np.newaxis] * velo * d)
        diff_pot12 = self.pot_multi[self.now_state] - \
                self.pot_multi[cand_state] 
        in_root = in_prd ** 2 + 2 * diff_pot12

        if in_root < 0: 
            # hopping failer because of the insuffiency of kinetic energy 
            self.woutp.write_trj_message("Hopping Failur [Kinetic Enerty Insuffiency]\n")
            return False
        else:
            self.woutp.write_trj_message("Hopping Success!! State has changed from {0} to {1}\n"\
                        .format(self.now_state+1,cand_state+1))
            self.woutp.write_trj_message("Velocitiies are adjusted\n")
            
            # renewal for energy and force
            # time integration restarts at the last configuration  
            n = sqrt(in_root) - in_prd       
            self.mol.set_velocities(velo + n * d / self.m_root[:,np.newaxis])
            self.now_state = cand_state 
            self.pot.set_now_state(self.now_state) 
            self.pot.calc() 
            return True 
    
    def check_ediff(self):
        #if np.any(np.diff(self.pot_multi) < self.tsh_ediff_thresh):
        #    self.dt = self.tsh_ediff_factor * self.dt_bak
        #else:
        #    self.dt = self.dt_bak
        if self.now_state < self.nrange -1 and self.pot_multi[self.now_state+1] - self.pot_multi[self.now_state] < self.tsh_ediff_thresh:
            self.dt = self.tsh_ediff_factor * self.dt_bak
        elif self.now_state > 0 and self.pot_multi[self.now_state] - self.pot_multi[self.now_state-1] < self.tsh_ediff_thresh:
            self.dt = self.tsh_ediff_factor * self.dt_bak
        else:
            self.dt = self.dt_bak
        self.delta_st = 0.5 * self.dt / self.tsh_times 

class TullySurfaceHopping_QMMM(TullySurfaceHopping):
    
    def __init__(self, mol_qm, mol_mm, pot_qm, pot_mm, pot_qmmm, dt, nstep, tsh_times= 10, restart = False, \
            tsh_ediff_thresh=0.05, tsh_ediff_factor=0.25, check_mdstop_dispersion = False, \
            limit_dispersion = 10.0, tlim = float('inf')):

        self.mol_mm = mol_mm
        self.pot_mm = pot_mm
        self.pot_qmmm = pot_qmmm 
        
        if self.pot_mm.get_check_pbc: 
            self.mol_mm.set_positions(self.pot_mm.get_pbc_adjusted(self.mol_mm.get_positions()))  
        
        TullySurfaceHopping.__init__(self,mol_qm ,pot_qm, dt, nstep, tsh_times, restart,\
            tsh_ediff_thresh, tsh_ediff_factor, check_mdstop_dispersion, limit_dispersion, tlim)

    def setup_output(self):
        self.woutp = WriteOutputMD_TSH_QMMM()
        if self.restart: self.woutp.restart(self)
        else: self.woutp.start(self)

    def prepare(self):
        #self.pot.calc(self.count); self.pot_mm.calc(); self.pot_qmmm.calc()   
        self.pot.calc(); self.pot_mm.calc(); self.pot_qmmm.calc()   
        self.pre_judge_tsh()
        self.check_ediff()

    def step(self):
        dt = self.dt
        get_accelerations_save =  self.mol.get_accelerations()
        self.mol.set_positions(self.mol.get_positions() + \
                self.mol.get_velocities() *  dt + 0.5 * \
                get_accelerations_save * dt * dt)
        get_accelerations_save_mm =  self.mol_mm.get_accelerations()
        self.mol_mm.set_positions(self.mol_mm.get_positions() + \
                self.mol_mm.get_velocities() *  dt + 0.5 * \
                get_accelerations_save_mm * dt * dt)
        if self.pot_mm.get_check_pbc: 
            self.mol_mm.set_positions(self.pot_mm.get_pbc_adjusted(self.mol_mm.get_positions()))  
        #self.pot.calc(self.count); self.pot_mm.calc(); self.pot_qmmm.calc()   
        self.pot.calc(); self.pot_mm.calc(); self.pot_qmmm.calc()   
        self.mol.set_velocities(self.mol.get_velocities() + \
                0.5 * (self.mol.get_accelerations() + \
                get_accelerations_save) * dt)
        self.mol_mm.set_velocities(self.mol_mm.get_velocities() + \
                0.5 * (self.mol_mm.get_accelerations() + \
                get_accelerations_save_mm) * dt)
        self.elaptime += dt 

class TullySurfaceHopping_QMMM_Rigid(TullySurfaceHopping):
    
    def __init__(self, mol_qm, mol_mm, pot_qm, pot_mm, pot_qmmm, dt, nstep, tsh_times, restart = False, \
            tsh_ediff_thresh=0.05, tsh_ediff_factor=0.25, check_mdstop_dispersion = False, limit_dispersion = 10.0, \
            tlim = float('inf')):
        self.mol_mm = mol_mm
        self.pot_mm = pot_mm
        self.pot_qmmm = pot_qmmm 
        
        TullySurfaceHopping.__init__(self,mol_qm ,pot_qm, dt, nstep, tsh_times, restart,\
            tsh_ediff_thresh, tsh_ediff_factor, check_mdstop_dispersion, limit_dispersion, tlim)

    def setup_output(self):
        self.woutp = WriteOutputMD_TSH_QMMM()
        if self.restart: self.woutp.restart(self)
        else: self.woutp.start(self)

    def prepare(self):
        self.pot.calc(); self.pot_mm.calc(); self.pot_qmmm.calc()   
        self.pre_judge_tsh()
        self.check_ediff()

    def step(self):
        dt = self.dt
        get_accelerations_save =  self.mol.get_accelerations()
        self.mol.set_positions(self.mol.get_positions() + \
                self.mol.get_velocities() *  dt + 0.5 * \
                get_accelerations_save * dt * dt)
        self.pot.calc(); self.pot_qmmm.calc()   
        self.mol.set_velocities(self.mol.get_velocities() + \
                0.5 * (self.mol.get_accelerations() + \
                get_accelerations_save) * dt)
        self.elaptime += dt 

