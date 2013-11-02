from MakeInitial import SetLattice 
from Potential import  Potential_MM 
from MonteCarlo import MonteCarlo_Ex 

# Set the initial coordinate
# lattice structure of 500 Ar atoms
n = 500 
lattice = SetLattice("Ar",n,1.77)
vlength = lattice.get_lattice_length() 
mol = lattice.get_molecule() 

# define the potential  
pot = Potential_MM(mol, check_pbc=True, vlength=vlength) 

# make the instance
delta = [0.2] * (len(mol)) 
mc = MonteCarlo_Ex(mol, pot, delta, 100000, 30.0, restart=False)  

# adjust the frequency of the output reloading 
mc.access_writeoutput().set_freq_xyz(10000)
mc.access_writeoutput().set_freq_trajectory(10000) 
mc.access_writeoutput().set_freq_energy(10000) 

# run montecarlo 
mc.run() 

