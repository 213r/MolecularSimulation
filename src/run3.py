from IO_MOLPRO import molpro_input_parser
from Potential import Potential_TSH_CASSCF, Potential_MM, Potential_QMMM
from TullySurfaceHopping import TullySurfaceHopping_QMMM
from Constants import fs2tau, bohr2ang
import sys
from IO_MOLPRO import OutputMOLPRO
from MakeInitial import SetLattice, SetMaxwell, superpositioned_atoms_delete 

mol_qm, inp = molpro_input_parser("template.com") 
mol_qm.read_coord_from_file("coord1") 
mol_qm.read_velocity_from_file("velocity1")
pot_qm = Potential_TSH_CASSCF(mol_qm, inp, now_state=2, nrange=4) 

n = 500 
lattice = SetLattice("Ar",n,1.77)
vlength = lattice.get_lattice_length() 
mol_mm = lattice.set_molecule() 
SetMaxwell(mol_mm, 300).set_velocities()
mol_mm = superpositioned_atoms_delete(mol_mm, mol_qm)
pot_mm = Potential_MM(mol_mm, rlimit = vlength / 2) 
pot_qmmm = Potential_QMMM(mol_qm, mol_mm, rlimit = vlength / 2) 
tsh = TullySurfaceHopping_QMMM(mol_qm, mol_mm, pot_qm, pot_mm,\
        pot_qmmm, dt=0.5*fs2tau, nstep=5000,tsh_times=5)
tsh.run() 
