from IO_MOLPRO import molpro_input_parser
from Potential import Potential_QM_MP2, Potential_QM_RHF, Potential_MM, Potential_QMMM 
from Constants import fs2tau, bohr2ang
import sys
from IO_MOLPRO import OutputMOLPRO
from MakeInitial import SetLattice, SetMaxwell, superpositioned_atoms_delete 
from Molecule import  read_positions_formated, read_velocities_formated
from VelocityVerlet import  VelocityVerlet, VelocityVerlet_QMMM

mol_qm, inp = molpro_input_parser("template_mp2.com") 
mol_qm.read_coord_from_file("coord1") 
mol_qm.read_velocity_from_file("velocity1")
#pot_qm = Potential_TSH_CASSCF(mol_qm, inp, now_state=2, nrange=4) 
pot_qm = Potential_QM_MP2(mol_qm, inp) 

n = 500 
lattice = SetLattice("Ar",n,1.77)
vlength = lattice.get_lattice_length() 
mol_mm = lattice.get_molecule() 
mol_mm = read_positions_formated("md_save.dat") 
mol_mm.set_velocities(read_velocities_formated("md_save.dat").get_velocities())
#print mol_mm.get_velocities_formated() 
#print mol_mm.get_positions_formated() 
#sys.exit() 
#SetMaxwell(mol_mm, 30).set_velocities()
mol_mm = superpositioned_atoms_delete(mol_mm, mol_qm)
pot_mm = Potential_MM(mol_mm, check_pbc = True, vlength = vlength) 
pot_qmmm = Potential_QMMM(mol_qm, mol_mm, check_pbc = True, vlength = vlength) 

vel =  VelocityVerlet_QMMM(mol_qm, mol_mm, pot_qm, pot_mm, pot_qmmm, 0.5*fs2tau, nstep=50000, restart=False)
#vel =  VelocityVerlet(mol_mm, pot_mm, 0.5*fs2tau, nstep=50000, restart=False)
vel.access_writeoutput().set_freq_xyz(10)
vel.access_writeoutput().set_freq_energy(10) 
vel.access_writeoutput().set_freq_trajectory(10) 
vel.run() 
#pot_qmmm = Potential_QMMM(mol_qm, mol_mm, rlimit = vlength / 2) 
#tsh = TullySurfaceHopping_QMMM(mol_qm, mol_mm, pot_qm, pot_mm,\
#        pot_qmmm, dt=0.5*fs2tau, nstep=5000,tsh_times=5)
#tsh.run() 
