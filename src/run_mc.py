from Molecule import bind_molecule 
from IO_MOLPRO import molpro_input_parser
from MakeInitial import SetLattice, SetMaxwell, superpositioned_atoms_delete 
from Potential import  Potential_MM 
from MonteCarlo import MonteCarlo, MonteCarlo_Ex 
import numpy as np 
import sys 

mol_qm, inp = molpro_input_parser("template.com") 
mol_qm.read_coord_from_file("coord1") 

n = 500 
lattice = SetLattice("Ar",n,1.77)
vlength = lattice.get_lattice_length() 
mol_mm = lattice.get_molecule() 
mol_mm = superpositioned_atoms_delete(mol_mm, mol_qm)
mol_tot = bind_molecule(mol_qm,mol_mm) 
pot_tot = Potential_MM(mol_tot, check_pbc=True, vlength=vlength) 
#print mol_tot.get_positions_formated() 
#print mol_tot.get_atomnames()  
#sys.exit() 
print len(mol_tot) 

delta = [0.2] * (len(mol_tot))  
#mc = MonteCarlo(mol_tot, pot_tot, delta, 1000000, 30.0, restart=False, frozen_atom_number = [0,1,2,3,4])  
mc = MonteCarlo_Ex(mol_tot, pot_tot, delta, 1000000, 30.0, restart=False, treated_as_molecule = [0,1,2,3,4], delta_mol = [0.02,0.02])  
#mc = MonteCarlo_Ex(mol_tot, pot_tot, delta, 100000, 30.0, restart=False,frozen_atom_number = [0,1,2,3,4])  
mc.access_writeoutput().set_freq_xyz(1000)
mc.access_writeoutput().set_freq_trajectory(1000) 
mc.access_writeoutput().set_freq_energy(1000) 
mc.run() 

