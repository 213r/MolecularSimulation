""" 
    Verlcity Verlet method is used in Molcular Dynaics(MD).  
        This class inherits from Molecular Dynamics class. 
        All used variables are in atomic units   
"""

from MolecularDynamics import MolecularDynamics
from WriteOutput import WriteOutputMD, WriteOutputMD_QMMM

class VelocityVerlet(MolecularDynamics):
    
    """
    VelocityVerlet(self, mol, pot,  deltat, nstep, restart, check_mdstop_dispersion, tlim)
        ** This is performing QM or MM MD-Simulation.  
            
            mol => Molecular instance. Refer to the Molecule.py 
            pot => Potential instance. Refer to the Potential.py 
            deltat => time step 
            nstep => step number  
            restart => whether md starts in the middle. default is False 
                       if restart is True, the initial condition is loaded from the restart file 
                       and the output files are written continuously. 
            check_mdstop_dispersion => when atoms are widely dispersed, md is terminated.
                                        default is False 
            tlim => set the time when md is finished. 
                    default is inf. Usually, the run of md is controlled by nstep.  
    """
    
    def __init__(self, mol, pot, dt, nstep, restart=False, check_mdstop_dispersion = False,\
                 limit_dispersion = 20.0, tlim = float('inf')):
        MolecularDynamics.__init__(self, mol, pot, dt, nstep, restart, check_mdstop_dispersion, limit_dispersion, tlim)
        if self.pot.get_check_pbc: self.mol.set_positions(self.pot.get_pbc_adjusted(self.mol.get_positions())) 
        self.count = 0 
        self.nrange = self.pot.get_nrange()
        self.setup_output() 
            
    def setup_output(self):
        self.woutp = WriteOutputMD()
        self.woutp.start(self)
        if self.restart: self.woutp.restart(self)

    def run(self):
        #this is the potential energy caluculation at the initial point 
        self.pot.calc() 
        
        #this is for logging at the initial point 
        self.woutp.logging(self)  

        # loop for MD 
        for _ in xrange(self.nstep):
            self.step() 
            self.count += 1
            self.woutp.logging(self)
            self.mdstop_time_trans()
            if self.mdstop_time_trans(): break
            if self.check_mdstop_dispersion and self.mdstop_atom_dispersion(): break
        self.woutp.finalize(self)
 
    def step(self):
        accelerations_save =  self.mol.get_accelerations()
        print self.elaptime, self.mol.get_positions(), self.mol.get_velocities(),accelerations_save 
        positions = self.mol.get_positions() + \
                self.mol.get_velocities() *  self.dt + 0.5 * \
                accelerations_save * self.dt * self.dt
        if self.pot.get_check_pbc(): self.mol.set_positions(self.pot.get_pbc_adjusted(positions)) 
        else:  self.mol.set_positions(positions) 
        self.pot.calc() 
        self.mol.set_velocities(self.mol.get_velocities() + \
                0.5 * (self.mol.get_accelerations() + \
                accelerations_save) * self.dt)
        self.elaptime += self.dt 

class VelocityVerlet_QMMM(MolecularDynamics):
    """ 
        VelocityVerlet_QMMM(self, mol_qm, mol_mm, pot_qm, pot_mm, pot_qmmm, dt, nstep, \
                            restart, check_mdstop_dispersion, tlim)
        ** This is performing QM/MM MD-Simulation.  
           The necessary arguments are alomost same with the above. 
           However, when you put the molecule and potential instance, you must separate 
           QM and MM parts. And, you must add the pot_qmmm instnace, which handles the QM/MM
           interface. 
    """ 
    
    def __init__(self, mol, mol_mm, pot, pot_mm, pot_qmmm, dt, nstep, \
            restart=False, check_mdstop_dispersion = False, limit_dispersion = 20.0,tlim = float('inf')):
        MolecularDynamics.__init__(self, mol, pot, dt, nstep, restart, check_mdstop_dispersion, limit_dispersion, tlim)
        self.mol_mm = mol_mm
        self.pot_mm = pot_mm
        self.pot_qmmm = pot_qmmm 
        #print "position-3" 
        #print self.mol_mm.get_positions() 
        #print self.pot_mm.get_check_pbc 
        if self.pot_mm.get_check_pbc(): 
            self.mol_mm.set_positions(self.pot_mm.get_pbc_adjusted(self.mol_mm.get_positions()))  
        self.count = 0
        self.nrange = self.pot.get_nrange()
        self.setup_output() 
        #print "position-2" 
        #print self.mol_mm.get_positions() 

    def run(self):
        self.pot.calc(); self.pot_mm.calc(); self.pot_qmmm.calc()
        #print "position-1" 
        #print self.mol_mm.get_positions() 
        #this is for logging at initial point 
        
        self.woutp.logging(self)  

        # loop for MD 
        while self.count < self.step:
            self.count += 1 
            self.step()
            self.woutp.logging(self)
            self.mdstop_time_trans()
            if self.mdstop_time_trans(): break
            if self.check_mdstop_dispersion and self.mdstop_atom_dispersion(): break

        self.woutp.finalize(self)
 
    def setup_output(self):
        self.woutp = WriteOutputMD_QMMM()
        self.woutp.start(self)
        if self.restart: self.woutp.restart(self)

    def step(self):
        #print "position0" 
        #print self.mol_mm.get_positions() 
        dt = self.dt
        get_accelerations_save = self.mol.get_accelerations()
        self.mol.set_positions(self.mol.get_positions() + \
                self.mol.get_velocities() *  dt + 0.5 * \
                get_accelerations_save * dt * dt)
        get_accelerations_save_mm =  self.mol_mm.get_accelerations()
        #print "position1" 
        #print self.mol_mm.get_positions() 
        positions_mm = self.mol_mm.get_positions() + \
                self.mol_mm.get_velocities() *  dt + 0.5 * \
                get_accelerations_save_mm * dt * dt
        #print "position2" 
        #print self.mol_mm.get_positions() 
        if self.pot_mm.get_check_pbc(): 
            self.mol_mm.set_positions(self.pot_mm.get_pbc_adjusted(positions_mm)) 
        else:  self.mol_mm.set_positions(positions_mm) 
        self.pot.calc(); self.pot_mm.calc(); self.pot_qmmm.calc()
        self.mol.set_velocities(self.mol.get_velocities() + \
                0.5 * (self.mol.get_accelerations() + \
                get_accelerations_save) * dt)
        self.mol_mm.set_velocities(self.mol_mm.get_velocities() + \
                0.5 * (self.mol_mm.get_accelerations() + \
                get_accelerations_save_mm) * dt)
        self.elaptime += dt 

