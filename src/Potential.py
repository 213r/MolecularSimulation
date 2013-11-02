import numpy as np
import pwd, os
from copy import deepcopy
from math import sqrt 
from IO_MOLPRO import InputMOLPRO, OutputMOLPRO
#from Potential_MM_Cython import potential_mm_cython, potential_qmmm_cython  
from Potential_MM_f2py import potential_mm_f2py, potential_mm_f2py_mc, potential_qmmm_f2py

"""
Potential Module 
    Wrapper for obtaining the potential energy and the gradinet  

"""

class Potential_MM:

    def __init__(self,mol, rlimit = 50.0, check_pbc = False, vlength = None ):
        self.mol = mol
        self.rlimit = rlimit 
        self.check_pbc = check_pbc 
        self.box = np.zeros((3,3)) 
        self.invbox = np.zeros((3,3))  
        if self.check_pbc and vlength != None: self.set_pbc(vlength)  

    def set_pbc(self, vlength):
        self.rlimit = vlength / 2 
        self.box = np.diag(np.ones(3) * vlength)
        self.invbox = np.linalg.inv(self.box) 

    def get_check_pbc(self):
        return self.check_pbc

    def get_pbc_adjusted(self, r): 
        s = np.dot(self.invbox, r.T)  
        return np.dot(self.box, s - np.round(s)).T  

    def calc(self):
        positions = self.mol.get_positions() 
        atomnumbers = self.mol.get_atomnumbers()
        #self.energy, self.forces = potential_mm_cython(positions, atomnumbers, self.rlimit) 
        self.energy, forces = potential_mm_f2py(atomnumbers, positions.T, self.rlimit, self.check_pbc, self.box, self.invbox) 
        self.forces = forces.T 
        self.mol.set_potential_energy(self.energy)
        self.mol.set_forces(self.forces) 

    def get_interaction_energy(self, ind):
        #This method is only for MonteCarlo 
        positions = self.mol.get_positions() 
        atomnumbers = self.mol.get_atomnumbers()
        return potential_mm_f2py_mc(ind, atomnumbers, positions.T, self.rlimit, self.check_pbc, self.box, self.invbox) 

    def get_potential_energy(self):
        return self.energy 
    
    def get_forces(self):
        return self.forces

class Potential_QMMM(Potential_MM):

    def __init__(self, mol_qm, mol_mm, rlimit = 50.0, check_pbc = False, vlength = None ):
        self.mol_qm = mol_qm
        self.mol_mm = mol_mm
        self.rlimit = rlimit 
        self.check_pbc = check_pbc 
        self.box = np.zeros((3,3)) 
        self.invbox = np.zeros((3,3))  
        if self.check_pbc and vlength != None: self.set_pbc(vlength)  

    def calc(self):
        positions_qm = self.mol_qm.get_positions() 
        atomnumbers_qm = self.mol_qm.get_atomnumbers()
        positions_mm = self.mol_mm.get_positions() 
        atomnumbers_mm = self.mol_mm.get_atomnumbers()
        self.energy, forces_qm, forces_mm = potential_qmmm_f2py(atomnumbers_qm, \
                atomnumbers_mm, positions_qm.T, positions_mm.T, self.rlimit, \
                self.check_pbc, self.box, self.invbox)
        self.mol_qm.set_forces(self.mol_qm.get_forces() + forces_qm.T) 
        self.mol_mm.set_forces(self.mol_mm.get_forces() + forces_mm.T) 

    def get_interaction_energy(self):
        return self.energy 

class Potential_QM:
    
    def __init__(self, mol, inp):
        self.mol = mol
        self.inp = inp
        if str(self.inp) == "InputMOLPRO":
            root, ext = os.path.splitext(self.inp.get_inputname())
            self.outp = OutputMOLPRO(root + ".out",self.mol)

    def override_output(self, wfile,mol=None):
        if mol is None: mol = self.mol 
        self.outp = OutputMOLPRO(wfile,self.mol)

    def calc(self):
        self.inp.set_molecule(self.mol) 
        self.inp.make_input()
        self.inp.write()
        os.system(self.inp.get_command())            
        self.mol.set_potential_energy(self.get_potential_energy())
        self.mol.set_forces(-self.get_gradient())

    def get_potential_energy(self):
        pass 

    def get_gradient(self):
        pass


class Potential_QM_RHF(Potential_QM):

    def __init__(self, mol, inp):
        Potential_QM.__init__(self, mol, inp)

    def get_potential_energy(self):
        return self.outp.get_potential_energy_rhf()

    def get_gradient(self):
        return self.outp.get_gradient_multi(0) 

class Potential_QM_MP2(Potential_QM):

    def __init__(self, mol, inp):
        Potential_QM.__init__(self, mol, inp)

    def get_potential_energy(self):
        return self.outp.get_potential_energy_mp2()

    def get_gradient(self):
        return self.outp.get_gradient_multi(0) 


class Potential_QM_CASSCF(Potential_QM):

    def __init__(self, mol, inp, now_state, nrange):
        Potential_QM.__init__(self, mol, inp)

        self.now_state = now_state
        self.nrange = nrange
        self.inp.set_cpmcscf_istate(self.now_state,self.nrange) 
    
    def get_potential_energy(self):
        return self.outp.get_potential_energy_mcscf_multi(self.nrange)[self.now_state]

    def get_gradient(self):
        return self.outp.get_gradient_multi(self.now_state) 


class Potential_QM_CASPT2(Potential_QM):
    
    def __init__(self, mol, inp, now_state, nrange):
        Potential_QM.__init__(self, mol, inp)

        self.now_state = now_state
        self.nrange = nrange
        self.inp.set_caspt2_istate(self.now_state,self.nrange) 
    
    def get_potential_energy(self):
        ene = self.outp.get_data_caspt2_multi(self.nrange)[0]
        return ene[self.now_state]

    def get_gradient(self):
        return self.outp.get_gradient_multi(self.now_state) 

class Potential_TSH(Potential_QM):

    def __init__(self, mol, inp, now_state, nrange):
        #super(Potential_TSH, self).__init__(mol, inp)
        Potential_QM.__init__(self, mol, inp)
        self.now_state = now_state
        self.nrange = nrange
        self.d_multi_old = None
    
        self.set_now_state(self.now_state) 

    def calc(self):
        self.inp.set_molecule(self.mol) 
        self.inp.make_input()
        self.inp.write()
        os.system(self.inp.get_command())
        self.check_nacme_sign() 
        self.mol.set_potential_energy_multi(self.get_potential_energy_multi())
        self.mol.set_forces(-self.get_gradient())

    def set_now_state(self,now_state):
        pass 

    def get_potential_energy_multi(self):
        pass  

    def check_nacme_sign(self):
        self.d_multi = self.outp.get_nacme(self.nrange)
        if self.d_multi_old == None:
            self.d_multi_old = self.d_multi
            return 
          
        for i in xrange(self.nrange-1):
            for j in xrange(i+1,self.nrange):
                sum_d = np.sum(self.d_multi[i,j] * \
                        self.d_multi_old[i,j])
                if sum_d < 0.0:
                    self.d_multi[i,j] = - self.d_multi[i,j]
                    self.d_multi[j,i] = - self.d_multi[j,i]
        self.d_multi_old = self.d_multi

    def get_gradient(self):
        return self.outp.get_gradient_multi(self.now_state) 
    
    def get_nrange(self):
        return self.nrange

    def get_now_state(self):
        return self.now_state

    def get_nacme_multi(self):
        return self.d_multi 

    def get_inner_v_d(self):
        pass

class Potential_TSH_CASSCF(Potential_TSH):

    def __init__(self, mol, inp, now_state, nrange):
        #super(Potential_TSH_CASSCF, self).__init__(mol, inp, now_state, nrange)
        Potential_TSH.__init__(self, mol, inp, now_state, nrange)

    def set_now_state(self,now_state):
        self.now_state = now_state
        self.inp.set_cpmcscf_istate(self.now_state,self.nrange) 

    def get_potential_energy_multi(self):
        return self.outp.get_potential_energy_mcscf_multi(self.nrange)

    def get_inner_v_d(self):
        # the sign of nacv must be checked befor it is used 
        return np.einsum("pq,ijpq -> ij",self.mol.get_velocities(), self.d_multi)

    def get_velocity_adjustment_vecotr(self): 
        return self.d_multi  

class Potential_TSH_CASPT2(Potential_TSH):

    def __init__(self, mol, inp, now_state, nrange):
        Potential_TSH.__init__(self, mol, inp, now_state, nrange)
        #super(Potential_TSH_CASPT2, self).__init__(mol, inp, now_state, nrange)
        self.setup()

    def set_now_state(self,now_state):
        self.now_state = now_state
        self.inp.set_caspt2_istate(self.now_state,self.nrange) 

    def setup(self):
        self.eps = 0.01
        
        dinp = deepcopy(self.inp)
        dinp.set_inputname("tmp2.com")
        s = dinp.get_method_save().replace("<force>\n","")
        s = s.replace("<cpmcscf>\n","").replace("force\n","") 
        dinp.set_method_save(s)
        dinp.set_caspt2_istate(self.now_state, self.nrange) 
        self.dinp = dinp
        if str(self.dinp) == "InputMOLPRO":
            root, ext = os.path.splitext(self.dinp.get_inputname())
            self.doutp = OutputMOLPRO(root + ".out",self.mol)

    def get_potential_energy_multi(self):
        self.ene, self.co, self.s = self.outp.get_data_caspt2_multi(self.nrange)   
        return  self.ene
    
    def get_inner_v_d(self):
        # the sign of nacv must be checked befor it is used 
        # refert to eq.13
        pos_save,vel_save = self.mol.get_positions(), self.mol.get_velocities()
        len_vel =  sqrt(np.sum(vel_save)**2)
        self.mol.set_positions(pos_save + self.eps * vel_save / len_vel) 
        self.dinp.set_molecule(self.mol)  
        self.dinp.make_input()
        self.dinp.write()
        os.system(self.dinp.get_command())
        self.mol.set_positions(pos_save) 
        ene2, co2, s2 = self.doutp.get_data_caspt2_multi(self.nrange)   
        dco = (co2 - self.co) / self.eps 
        # refert to eq.9
        v_d1 = np.einsum("ip, jq, pq -> ij",self.co,dco,self.s)
        v_d2 = np.einsum("ip,pqlm,jq -> ijlm",self.co,self.d_multi,self.co)
        return len_vel * v_d1 + np.einsum("pq,ijpq -> ij",vel_save,v_d2) 

    def get_velocity_adjustment_vecotr(self): 
        return self.mol.get_forces() 

#class Potential_for_tully:
#    from TestPotential_TSH_type1 import nac, adiabatic, diff_adiabatic 
#
#    def __init__(self, mol,now_state, nrange):
#        self.mol = mol
#        self.now_state = now_state
#        self.nrange = nrange
#
#    def calc(self):
#        self.mol.set_potential_energy(self.get_potential_energy())
#        self.mol.set_forces(-self.get_gradient())
#
#    def get_potential_energy_multi(self):
#        return adiabatic(self.mol.get_positions())      
#
#    def get_potential_energy(self):
#        return adiabatic(self.mol.get_positions())[self.now_state]      
#
#    def get_gradient(self):
#        return diff_adiabatic(self.mol.get_positions())[self.now_state]      
#
#    def get_nacme_multi(self):
#        d12 = nac(self.mol.get_positions())     
#        return  np.array([[0.0,d12],[-d12,0.0]])     
#
#    def get_nrange(self):
#        return self.nrange
#
#    def get_now_state(self):
#        return self.now_state
#
#
#class Potential:
#
#    """    
#    To obtain the potential energy and forces of 
#    target molecules 
#    First, this class submits job file 
#    ,read the output class
#    and then, get the necessary information. 
#
#    In present, this covers only Molpro program  
#
#    """
#
#    def __init__(self,mol):
#        # set the target molecule as molecule instance 
#        self.mol = mol
#        self.user = pwd.getpwuid(os.getuid())[0]   
#        self.cwd = os.getcwd() 
#        #self.root, self.ext = os.path.splitext(self.inp.get_inputname())
#
#    def calc(self):
#        self.inp.set_molecule(self.mol) 
#        self.write_input() 
#        os.system(self.get_command())            
#        self.energy = self.outp.get_potentialenergy()
#        self.forces = -self.outp.get_gradient()
#        self.mol.set_potential_energy(self.energy[0])
#        self.mol.set_forces(self.forces[0])
#
#    def calc_tsh(self,nrange,now_state):
#        self.inp.set_molecule(self.mol) 
#        self.write_input() 
#        os.system(self.get_command())            
#        self.mol.set_potential_energy(self.get_potential_energy()[now_state])
#        self.mol.set_forces(-self.get_gradient(now_state))
#
#    def get_nacme_multi(self, nrange):
#        return self.outp.get_nacme(nrange)
#
#    def get_potential_energy(self):
#        return self.outp.get_potential_energy()
#    
#    def get_gradient(self, now_state):
#        return self.outp.get_gradient(now_state)
#
#    def write_input(self):
#        self.inp.make_input()
#        self.inp.write()
#
#    def get_command(self):
#        # This method must be overrided. 
#        pass
#
#def test():
#    from Atom import Atom
#    from Molecule import Molecule 
#    from IO_MOLPRO import InputMOLPRO 
#    at1 = Atom(symbol='H', position=[2,0,0])
#    at2 = Atom(symbol='H', position=[0,0,0])
#    mol = Molecule([at1,at2])
#    inp = InputMOLPRO()
#    pot = Potential_QM(mol,inp)
#    pot.calc()
#
#def test1():
#    from Atom import Atom
#    from Molecule import Molecule 
#    from TullySurfaceHopping2 import TullySurfaceHopping 
#    from Constants import fs2tau, tau2fs 
#
#    at1 = Atom(symbol='X',mass=2000)
#    mol1 = Molecule([at1])
#    mol1.set_velocities([4.0]) 
#    mol1.set_positions([-10]) 
#    pot = Potential_for_tully(mol1,now_state=0,nrange=2)
#    tsh = TullySurfaceHopping(mol1,pot,dt=0.5*fs2tau,nstep=100,\
#            tsh_times=5)  # These are special fort TSH
#
#    tsh.judge_tsh()
#
#def test2():
#    from Atom import Atom
#    from Molecule import Molecule 
#    from TullySurfaceHopping import TullySurfaceHopping 
#    from Constants import fs2tau, tau2fs 
#    from IO_MOLPRO import molpro_input_parser
#    from VelocityVerlet import set_random_momentum  
#
#    mol, inp = molpro_input_parser("mlpinp1.com") 
#    set_random_momentum(mol,0.1)  
#    pot = Potential_TSH(mol,inp,now_state=0,nrange=3)
#    pot.override_output("mlpinp0.out_01") 
#    read_velocity(mol,"velocity5_1") 
#    read_coord(mol,"coord5_1") 
#    tsh = TullySurfaceHopping(mol,pot,dt=0.5*fs2tau,nstep=100,\
#            tsh_times=5)  # These are special fort TSH
#
#    tsh.run()
#    #tsh.check_nacme_sign()
#    #tsh.judge_tsh()
#
#def energy_arar(r): 
#    eps = 379386396928e-15    
#    sigma = 6.43640672
#    num = sigma / r
#    num = num**6 
#    ene = - num 
#    ene += num * num
#    return 4.0 * eps * ene
#
#def diff_energy_arar(r):
#    eps = 379386396928e-15    
#    sigma = 6.43640672
#    num = sigma / r 
#    return -24.0 * eps * (2 * num**13 - num ** 7) / sigma 
#
#
#def test3():
#    from Atom import Atom
#    from Molecule import Molecule
#    from VelocityVerlet import VelocityVerlet
#    from Constants import fs2tau
#
#    at1 = Atom(symbol='Ar', position=[5,0,0])
#    at2 = Atom(symbol='Ar', position=[0,0,0])
#    at3 = Atom(symbol='Ar', position=[0,5,0])
#    mol = Molecule([at1,at2,at3])
#    pot = Potential_MM(mol)
#    md = VelocityVerlet(mol,pot,deltat=0.5*fs2tau ,nstep=200)
#    md.run() 
#
#def potential(calc_type, mol, inp = None, tsh = False, \
#        now_state = None, nrange = None, mm_type = None):
#    if calc_type == "rhf": 
#        return Potential_QM_RHF(mol, inp)
#    elif calc_type == "casscf":
#        if tsh: return Potential_TSH_CASSCF(mol, inp, now_state, nrange)
#        else: return Potential_QM_CASSCF(mol, inp, now_state, nrange)
#    elif calc_type == "caspt2": 
#        if tsh: return Potential_TSH_CASPT2(mol, inp, now_state, nrange)
#        else: return Potential_QM_CASPT2(mol, inp, now_state, nrange)
#    elif calc_type == "mm":   
#        return Potential_MM(mol, mm_type) 
#    else: 
#        assert "specify the range of states or the current state"
#        sys.exit()
#
