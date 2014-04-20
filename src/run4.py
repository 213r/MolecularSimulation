from IO_MOLPRO import molpro_input_parser
from TullySurfaceHopping import TullySurfaceHopping
from Molecule import Molecule
from Potential import Potential_MM
from Constants import fs2tau, ang2bohr, bohr2ang
from VelocityVerlet import VelocityVerlet
from IO_MOLPRO import OutputMOLPRO
from MakeInitial import SetMaxwell, SetLattice 
import numpy as np
import sys

"""
n = 3
mol_mm = Molecule(["Ar"]*n,positions=np.array([[5.,0.,0.],[-5.,0.,0.],[0., 5., 0.]]))  
pot = Potential_MM(mol_mm)
pot.calc()
print pot.get_potential_energy() 
print pot.get_forces() 
sys.exit()
mol_mm = SetLattice("Ar",n,1.77).set_molecule() 
SetMaxwell(mol_mm,30).set_velocities()
#mol_mm = Molecule(["Ar"]*3,positions=np.array([[0.,0.,0.],[0.,0.,5.2*ang2bohr], [0., 2.5 * ang2bohr, 2.5*ang2bohr]]))  
pot_mm = Potential_MM(mol_mm) 
vel = VelocityVerlet(mol_mm, pot_mm, nstep=20000, dt=0.5*fs2tau)
vel.run() 
sys.exit()
"""

def do():
    n = 50 
    mol_mm = SetLattice("Ar",n,1.77).set_molecule() 
    SetMaxwell(mol_mm,30).set_velocities()
    pot_mm = Potential_MM(mol_mm) 
    vel = VelocityVerlet(mol_mm, pot_mm, nstep=200, dt=0.5*fs2tau)
    vel.run() 

import pstats, cProfile
import pyximport 
pyximport.install()
cProfile.runctx("do()", globals(), locals(), "Profile.prof")
s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()

sys.exit() 
for i in xrange(500):
    x = 1.0 + 0.1 * i 
    mol_mm.set_positions([[0.,0.,0.],[0.,0.,x*ang2bohr]]) 
    pot_mm.calc() 
    #print x, pot_mm.get_potential_energy(), pot_mm.get_forces()[0][2]
    print x,  pot_mm.get_forces()
#print pot_mm.get_forces() 
sys.exit() 
n = 256
mol_mm = Molecule(["Ar"]*n) 
SetLattice(mol_mm,1.77).set_positions() 
SetMaxwellVelocities(mol_mm,30).set_velocities()
pot_mm = Potential_MM(mol_mm) 
vel = VelocityVerlet(mol_mm, pot_mm, nstep=2000, dt=0.5*fs2tau)
vel.run() 
sys.exit() 
mol, inp = molpro_input_parser("template.com") 
read_coord(mol,"coord1") 
read_velocity(mol,"velocity1")
pot = Potential_TSH_CASSCF(mol, inp, now_state=1, nrange=2) 
tsh = TullySurfaceHopping(mol,pot,dt=0.5*fs2tau,\
       nstep=500,tsh_times=5) 
tsh.run()
