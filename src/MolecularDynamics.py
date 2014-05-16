import numpy as np 
from datetime import datetime
from Constants import ang2bohr 
from WriteOutput import WriteOutputMD
import sys

"""
Molecular Dynamics (MD) simulation program 

Author: Keisuke Niimi

This class must be overrided by any subclasse according to the supecific method.

"""

class MolecularDynamics:
    
    def __init__(self, mol, pot, dt, nstep, restart, check_mdstop_dispersion, tlim):
        
        self.mol, self.pot = mol, pot
        self.dt, self.nstep = dt, nstep 
        self.restart = restart
        self.check_mdstop_dispersion = check_mdstop_dispersion 
        self.elaptime = 0.0 
        self.tlim = tlim

    def access_writeoutput(self):
        return self.woutp

    def run(self):
        # this method should be override depending on your chonsen method 
        pass 

    def set_elaptime(self, elaptime):
        self.elaptime = elaptime 

    def step(self):
        # this method should be override depending on your chonsen method 
        pass 

    def mdstop_time_trans(self):
        if self.tlim < self.elaptime:
            self.woutp.finalize(self)
            self.woutp.write_log_message("MD is terminated because time limitation is transcended.")
            return True
        else: return False

    def mdstop_atom_dispersion(self,rlim = 10*ang2bohr): 
        coord =  self.mol.get_positions()
        for i in xrange(len(coord)-1):
            for j in xrange(i+1,len(coord)):
                x = coord[i] - coord[j] 
                if np.sum(x*x) > rlim * rlim:
                    self.woutp.finalize(self)
                    self.woutp.write_log_message("MD is terminated because of the large dispersion of atoms.")
                    return True
        return False

