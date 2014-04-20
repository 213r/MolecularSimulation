import numpy as np 
import math
from Molecule import Molecule, bind_molecule
from datetime import datetime
from Restart import restart_mc, restart_md, restart_md_qmmm
from Constants import ang2bohr,bohr2ang, fs2tau, tau2fs
import sys,os

class WriteOutput:
    # This class must be overrided. 
    # This class handles the frequency of the output reaload and xyz output printing . 

    def __init__(self):
 
        self.time_start = datetime.now()   # log the time when MD starts
        self.freq_xyz, self.freq_energy, self.freq_trj, self.freq_vel_xyz = 0, 0, 0, 0 
        self.count_xyz, self.count_energy, self.count_trj, self.count_vel_xyz = 0, 0, 0, 0 
        self.count = 0

    def set_freq_xyz(self, freq):
        self.freq_xyz = freq - 1
        if self.not_restart: self.count_xyz = self.freq_xyz
 
    def set_freq_vel_xyz(self, freq):
        self.freq_vel_xyz = freq - 1
        if self.not_restart: self.count_vel_xyz = self.freq_vel_xyz

    def set_freq_energy(self, freq):
        self.freq_energy = freq - 1
        if self.not_restart: self.count_energy = self.freq_energy

    def set_freq_trajectory(self, freq):
        self.freq_trj = freq - 1
        if self.not_restart: self.count_trj = self.freq_trj
    
    def start(self):
        if os.path.isfile(self.file_xyz): os.remove(self.file_xyz)  
        if os.path.isfile(self.file_ene): os.remove(self.file_ene)  
        if os.path.isfile(self.file_log): os.remove(self.file_log)  
        if os.path.isfile(self.file_trj): os.remove(self.file_trj)  
        self.not_restart = True 

    def restart(self):
        with open(self.file_ene,'a') as f: 
            f.write("#Restart Here!!\n")
        with open(self.file_log,'a') as f: 
            f.write("Restart Here!!\n")
        with open(self.file_trj,'a') as f: 
            f.write("Restart Here!!\n")
            f.write("\n"+"-"*40+"\n")
        self.not_restart = False

    def logging(self, self_simu):
        pass 

    def write_xyz(self,self_simu):
        if self.count_xyz == self.freq_xyz:
            self.count_xyz = 0
            if hasattr(self_simu, "elaptime"): mess = "No.{}, Time: {} fs".format(self_simu.count,self_simu.elaptime * tau2fs) 
            else: mess = None  
            with open(self.file_xyz,'a') as f: 
                f.write(self_simu.mol.get_positions_formated(unit="ang",label=False, message=mess))   
        else: self.count_xyz += 1 
    
    def write_trj(self, self_simu):
        if self.count_trj == self.freq_trj:
            self.count_trj = 0
            with open(self.file_trj,'a') as f: 
                f.write('No.{0}: Time = {1} [fs]\n'.format(str(self_simu.count),str(self_simu.elaptime*tau2fs)))
                f.write(self_simu.mol.get_positions_formated(unit='bohr'))  
                f.write(self_simu.mol.get_velocities_formated(unit='bohr/tau'))  
                f.write("\n"+"-"*40+"\n")
        else: self.count_trj += 1 

    def write_ene(self, self_simu):
        pass 

    def finalize(self, self_simu):
        pass 

    def write_trj_message(self,message):
        with open(self.file_trj,'a') as f: f.write(message) 
    
    def write_ene_message(self,message):
        with open(self.file_ene,'a') as f: f.write(message) 

    def write_log_message(self,message):
        with open(self.file_log,'a') as f: f.write(message)  
            
class WriteOutputMC(WriteOutput):

    def __init__(self): 

        WriteOutput.__init__(self) 
        self.file_xyz, self.file_ene = "mc.xyz", "mc_energy.dat" 
        self.file_log, self.file_trj = "mc.log", "mc.trj"
        self.file_restart = "mc_restart.dat"
    
    def start(self):
        WriteOutput.start(self) 
        self.write_log_message("MC starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        self.write_ene_message("#  Step Potential[H] Acceptance ratio\n")

    def restart(self, self_mc):
        WriteOutput.restart(self) 
        self.write_log_message("MC starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        positions = restart_mc(self.file_restart) 
        self_mc.mol.set_positions(positions)

    def logging(self, self_mc):
        self.write_xyz(self_mc)
        self.write_energy(self_mc)

    def write_energy(self,self_mc):
        if self.count_energy == self.freq_energy:
            self.count_energy = 0
            pot = self_mc.mol.get_potential_energy()
            if self_mc.count == 0: ratio = 0 
            else: 
                ratio = float(self_mc.count_accept) / self_mc.count_tot 
                if hasattr(self_mc, "delta_mol"):
                    ratio_mol = float(self_mc.count_accept_mol) / self_mc.count_tot_mol 
                    self.write_ene_message("{0:7d}   {1: 8.6f}   {2: 8.6f}   {3: 8.6f}\n".\
                                format(self_mc.count, pot, ratio, ratio_mol))
                else: 
                    self.write_ene_message("{0:7d}   {1: 8.6f}   {2: 8.6f}\n".format(self_mc.count, pot, ratio))
        else: self.count_energy += 1 

    def finalize(self, self_mc):
        with open(self.file_restart,'w') as f: 
            f.write('Step = {0} \n'.format(str(self_mc.count)))
            f.write(self_mc.mol.get_positions_formated(unit='bohr'))  
            if hasattr(self_mc, "delta_mol"):
                f.write("\nOriginal Coordinates of fixed molecule [borh]\n")
                f.write(self_mc.mol_save.get_positions_formated(unit='bohr'))  
                f.write("\nalpha beta gamma\n")
                f.write("{0: 10.8f} {1: 10.8f}  {2: 10.8f}\n".format(self_mc.alpha_save, self_mc.beta_save, self_mc.gamma_save))

        time_start = self.time_start ; time_end = datetime.now()
        with open(self.file_log,'a') as f: 
            f.write("MC ends at {0:%Y-%m-%d %H:%M:%S}\n".format(time_end))
            f.write("Total MC run takes {0}\n".format(str(time_end - time_start).split('.')[0]))

class WriteOutputMD(WriteOutput):

    def __init__(self):
  
        WriteOutput.__init__(self) 
        self.file_xyz, self.file_ene = "md.xyz", "md_energy.dat" 
        self.file_log, self.file_trj = "md.log", "md.trj"
        self.file_restart = "md_restart.dat"

    def start(self, self_md):
        WriteOutput.start(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        #self.write_ene_message("#Time[fs] Kinetic[H] Potential[H] Total[H]\n")
        self.write_ene_message("#Time[fs] Kin(QM)[H] Pot(QM)(1~{0}) Total\n"\
                .format(self_md.nrange))

    def restart(self, self_md):
        WriteOutput.restart(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        at_start, positions, velocities = restart_md(self.file_restart)  
        self_md.set_elaptime(at_start) 
        self_md.mol.set_positions(positions); self_md.mol.set_velocities(velocities) 

    def logging(self, self_md):
        self.write_xyz(self_md)
        self.write_trj(self_md)
        self.write_energy(self_md)

    def write_energy(self,self_md):
        if self.count_energy == self.freq_energy:
            self.count_energy = 0
            kin = self_md.mol.get_kinetic_energy()
            pot_multi = self_md.mol.get_potential_energy_multi()
            tot = kin + pot_multi[0]
            with open(self.file_ene,'a') as f: 
                f.write("{0: 4.2f} ".format(self_md.elaptime*tau2fs))
                f.write("{0: 8.6f} ".format(kin))
                for i in pot_multi: f.write(" {0: 8.6f} ".format(i))
                f.write("{0: 8.6f}\n".format(tot))
        else: self.count_energy += 1 
            
    def finalize(self, self_md):
        with open(self.file_restart,'w') as f: 
            f.write('Time = {0} [tau]\n'.format(str(self_md.elaptime)))
            f.write(self_md.mol.get_positions_formated(unit='bohr'))  
            f.write(self_md.mol.get_velocities_formated(unit='bohr/tau'))  

        time_start = self.time_start ; time_end = datetime.now()
        with open(self.file_log,'a') as f: 
            f.write("MD ends at {0:%Y-%m-%d %H:%M:%S}\n".format(time_end))
            f.write("Total MD run takes {0}\n".format(str(datetime.now() - self.time_start).split('.')[0]))

class WriteOutputMD_QMMM(WriteOutput):

    def __init__(self):
        WriteOutput.__init__(self) 
        self.file_xyz, self.file_ene = "md_qmmm.xyz", "md_qmmm_energy.dat" 
        self.file_log, self.file_trj = "md_qmmm.log", "md_qmmm.trj"
        self.file_restart = "md_qmmm_restart.dat"
        self.file_vel_xyz = "md_qmmm_vel.xyz"

    def start(self, self_md):
        WriteOutput.start(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        #self.write_ene_message("#Time[fs] Kin(QM)[H] Pot(QM) Kin(MM) Pot(MM) Pot(QMMM) Total\n")
        self.write_ene_message("#Time[fs] Kin(QM)[H] Pot(QM)(1~{0}) Kin(MM) Pot(MM) Pot(QMMM) Total\n"\
                .format(self_md.nrange))

    def restart(self, self_md):
        WriteOutput.restart(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        at_start, positions_qm, velocities_qm, positions_mm, velocities_mm = \
                restart_md_qmmm(self.file_restart)  
        self_md.set_elaptime(at_start) 
        self_md.mol.set_positions(positions); self_md.mol.set_velocities(velocities) 
        self_md.mol_mm.set_positions(positions_mm); self_md.mol_mm.set_velocities(velocities_mm) 

    def logging(self, self_md):
        self.write_xyz(self_md)
        self.write_vel_xyz(self_md)
        self.write_trj(self_md)
        self.write_energy(self_md)
 
    def write_xyz(self,self_md):
        if self.count_xyz == self.freq_xyz:
            self.count_xyz = 0
            mol_tot = bind_molecule(self_md.mol,self_md.mol_mm) 
            mess = "No.{} Time: {} fs [ang]".format(self_md.count, self_md.elaptime * tau2fs) 
            with open(self.file_xyz,'a') as f: 
                f.write(mol_tot.get_positions_formated(unit="ang",label=False,message=mess))   
        else: self.count_xyz += 1 
  
    def write_vel_xyz(self,self_md):
        if self.count_vel_xyz == self.freq_vel_xyz:
            self.count_vel_xyz = 0
            mol_tot = bind_molecule(self_md.mol,self_md.mol_mm) 
            mess = "Time: {} fs  [bohr/tau]".format(self_md.elaptime * tau2fs) 
            with open(self.file_vel_xyz,'a') as f: 
                f.write(mol_tot.get_velocities_formated(unit="bohr/tau",label=False,message=mess))   
        else: self.count_vel_xyz += 1 
      
    def write_energy(self, self_md):
        if self.count_energy == self.freq_energy:
            self.count_energy = 0
            kin_qm = self_md.mol.get_kinetic_energy()
            pot_qm_multi = self_md.mol.get_potential_energy_multi()
            kin_mm = self_md.mol_mm.get_kinetic_energy()
            pot_mm = self_md.mol_mm.get_potential_energy()
            pot_qmmm = self_md.pot_qmmm.get_interaction_energy() 
            tot = kin_qm + pot_qm_multi[0] + kin_mm + pot_mm + pot_qmmm 
            with open(self.file_ene,'a') as f: 
                f.write("{0: 4.2f} ".format(self_md.elaptime*tau2fs))
                f.write("{0: 8.6f} ".format(kin_qm))
                for i in pot_qm_multi: f.write(" {0: 8.6f} ".format(i))
                f.write("{0: 8.6f} {1: 8.6f} ".format(kin_mm, pot_mm))
                f.write("{0: 8.6f} {1: 8.6f}\n".format(pot_qmmm, tot))
        else: self.count_energy += 1

    def finalize(self,self_md):
        with open(self.file_restart,'w') as f: 
            f.write('Time = {0} [tau]\n'.format(str(self_md.elaptime)))
            f.write('QM Part\n')
            f.write(self_md.mol.get_positions_formated(unit='bohr'))  
            f.write(self_md.mol.get_velocities_formated(unit='bohr/tau'))  
            f.write('\nMM Part\n')
            f.write(self_md.mol_mm.get_positions_formated(unit='bohr'))  
            f.write(self_md.mol_mm.get_velocities_formated(unit='bohr/tau'))  

        time_start = self.time_start
        time_end = datetime.now()
        with open(self.file_log,'a') as f: 
            f.write("MD ends at {0:%Y-%m-%d %H:%M:%S}\n".format(time_end))
            f.write("Total MD run takes {0}\n".format(str(time_end - time_start).split('.')[0]))

class WriteOutputMD_TSH(WriteOutputMD):
    "inherit write_xyz method only from WriteOutputMD" 
    
    def __init__(self):
    
        WriteOutput.__init__(self) 
        self.file_xyz, self.file_ene = "md_tsh.xyz", "md_tsh_energy.dat" 
        self.file_log, self.file_trj = "md_tsh.log", "md_tsh.trj"
        self.file_restart = "md_tsh_restart.dat"
        self.file_prob = "tsh_prob.dat" 

    def start(self,self_tsh):
        WriteOutput.start(self) 
        if os.path.isfile(self.file_prob): os.remove(self.file_prob)  

        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        self.write_ene_message("#Time[fs] State Kinetic[H] Potential[H](1~{0},current) Total[H]\n"\
                .format(self_tsh.nrange))

    def restart(self, self_tsh):
        WriteOutput.restart(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
            
        #with open(self.file_restart,'r') as f: 
        #    natom = len(self_tsh.mol)
        #    positions, velocities, coefficients = [], [], [] 
        #    line = f.readline() 
        #    while line:
        #        if "Time" in line: at_start = float(line.split()[2])
        #        if "Current state" in line: now_state = int(line.split()[2]) - 1
        #        if "Coordinates" in line:
        #            for _ in xrange(natom): positions.append(map(float, f.readline().split()[1:]))
        #        if "Velocities" in line:
        #            for _ in xrange(natom): velocities.append(map(float, f.readline().split()[1:]))
        #        if "Coefficients" in line:
        #            for _ in xrange(self_tsh.nrange):
        #                num = f.readline().split() 
        #                coefficients.append(float(num[3]) + float(num[5]) * 1.j)
        #        line = f.readline()
        #    self_tsh.set_elaptime(at_start); self_tsh.pot.set_now_state(now_state) 
        #    self_tsh.mol.set_positions(positions); self_tsh.mol.set_velocities(velocities) 
        #    self_tsh.set_coefficients(coefficients)

    def write_energy(self, self_tsh):
        if self.count_energy == self.freq_energy:
            self.count_energy = 0
            kin = self_tsh.mol.get_kinetic_energy()
            pot = self_tsh.mol.get_potential_energy_multi()
            with open(self.file_ene,'a') as f: 
                f.write("{0: 4.2f} {1:d} {2: 8.6f} ".format(self_tsh.elaptime*tau2fs,\
                    self_tsh.now_state+1, kin))
                for i in pot: 
                    f.write(" {0: 8.6f} ".format(i))
                f.write(" {0: 8.6f} ".format(pot[self_tsh.now_state]))
                f.write(" {0: 8.6f}\n".format(kin + pot[self_tsh.now_state]))
        else: self.count_energy += 1

    def write_trj(self, self_tsh):
        if self.count_trj == self.freq_trj:
            self.count_trj = 0
            with open(self.file_trj,'a') as f: 
                f.write('No.{0}: Time = {1} [fs]\n'.format(str(self_tsh.count),str(self_tsh.elaptime*tau2fs)))
                f.write(self_tsh.mol.get_positions_formated(unit='bohr'))  
                f.write(self_tsh.mol.get_velocities_formated(unit='bohr/tau'))  
                f.write("\nCoefficeients of each state\n")   
                for i, ic in enumerate(self_tsh.c): 
                    f.write("St{}:  Re {}  Img {}\n".format(i+1, ic.real, ic.imag))   
                f.write("\nThe length of nacme\n")   
                for i in xrange(self_tsh.nrange-1): 
                    for j in xrange(i+1,self_tsh.nrange): 
                        nacme = self_tsh.pot.get_nacme_multi()[i,j]
                        f.write("St{} to St{}: {}\n".format(i+1, j+1, np.sum(nacme ** 2)))
                f.write("\n" + "-"*40 + "\n\n")
        else: self.count_trj += 1

    def write_prob_message(self,message):
        # log the hopping probability 
        with open(self.file_prob,'a') as f: 
            f.write(message)

    def finalize(self,self_tsh):
        with open(self.file_restart,'w') as f: 
            f.write('Time = {0} [tau]\n'.format(str(self_tsh.elaptime)))
            f.write('Current state: {0}\n'.format(str(self_tsh.now_state+1)))
            f.write(self_tsh.mol.get_positions_formated(unit='bohr'))  
            f.write(self_tsh.mol.get_velocities_formated(unit='bohr/tau'))  
            f.write('Coefficients\n')
            for i, ic in enumerate(self_tsh.c): 
                f.write("St: {}  Re: {}  Img: {}\n".format(i, ic.real, ic.imag))   

        time_start = self.time_start ; time_end = datetime.now()
        with open(self.file_log,'a') as f: 
            f.write("MD ends at {0:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()))
            f.write("Total MD run takes {0}\n".format(str(datetime.now() - self.time_start).split('.')[0]))

class WriteOutputMD_TSH_QMMM(WriteOutputMD_TSH):

    def __init__(self):
        WriteOutput.__init__(self) 
        self.file_xyz, self.file_ene = "md_tsh_qmmm.xyz", "md_tsh_qmmm_energy.dat" 
        self.file_log, self.file_trj = "md_tsh_qmmm.log", "md_tsh_qmmm.trj"
        self.file_restart = "md_tsh_qmmm_restart.dat"
        self.file_prob = "tsh_prob.dat" 

    def start(self,self_tsh):
        WriteOutput.start(self) 
        if os.path.isfile(self.file_prob): os.remove(self.file_prob)  

        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
        self.write_ene_message("#Time[fs] State Kin(QM)[H] Pot(QM)(1~{0},current) Kin(MM) Pot(MM) Pot(QMMM) Total\n"\
                .format(self_tsh.nrange))

    def restart(self, self_tsh):
        WriteOutput.restart(self) 
        self.write_log_message("MD starts at {0:%Y-%m-%d %H:%M:%S}\n".format(self.time_start))
                            
#        with open(self.file_restart,'r') as f: 
#            natom = len(self_tsh.mol)
#            positions, velocities, coefficients = [], [], [] 
#            natom_mm = len(self_tsh.mol_mm)
#            positions_mm, velocities_mm = [], [] 
#            line = f.readline() 
#            while line:
#                if "Time" in line: at_start = float(line.split()[2])
#                if "Current state" in line: now_state = int(line.split()[2]) - 1
#                if "Coordinates (QM) " in line:
#                    for _ in xrange(natom): positions.append(map(float, f.readline().split()[1:]))
#                if "Coordinates (MM) " in line:
#                    for _ in xrange(natom_mm): positions_mm.append(map(float, f.readline().split()[1:]))
#                if "Velocities (QM) " in line:
#                    for _ in xrange(natom): velocities.append(map(float, f.readline().split()[1:]))
#                if "Velocities (MM) " in line:
#                    for _ in xrange(natom_mm): velocities_mm.append(map(float, f.readline().split()[1:]))
#                if "Coefficients" in line:
#                    for _ in xrange(self_tsh.nrange):
#                        num = f.readline().split() 
#                        coefficients.append(float(num[3]) + float(num[5]) * 1.j)
#                line = f.readline()
#            self_tsh.set_elaptime(at_start); self_tsh.pot.set_now_state(now_state) 
#            self_tsh.mol.set_positions(positions); self_tsh.mol.set_velocities(velocities) 
#            self_tsh.mol_mm.set_positions(positions_mm); self_tsh.mol_mm.set_velocities(velocities_mm) 
#            self_tsh.set_coefficients(coefficients)

    def write_xyz(self,self_tsh):
        if self.count_xyz == self.freq_xyz:
            self.count_xyz = 0
            mol_tot = bind_molecule(self_tsh.mol,self_tsh.mol_mm) 
            mess = "No.{} Time: {} fs [ang]".format(self_tsh.count, self_tsh.elaptime * tau2fs) 
            with open(self.file_xyz,'a') as f: 
                f.write(mol_tot.get_positions_formated(unit="ang",label=False,message=mess))   
        else: self.count_xyz += 1 

    def write_energy(self, self_tsh):
        if self.count_energy == self.freq_energy:
            self.count_energy = 0
            kin = self_tsh.mol.get_kinetic_energy()
            pot = self_tsh.mol.get_potential_energy_multi()
            kin_mm = self_tsh.mol_mm.get_kinetic_energy() 
            pot_mm = self_tsh.mol_mm.get_potential_energy()
            pot_qmmm = self_tsh.pot_qmmm.get_potential_energy() 
            with open(self.file_ene,'a') as f: 
                f.write("{0: 4.2f} {1:d} {2: 8.6f} ".format(self_tsh.elaptime*tau2fs,\
                    self_tsh.now_state+1, kin))
                for i in pot:  f.write(" {0: 8.6f} ".format(i))
                f.write(" {0: 8.6f} ".format(pot[self_tsh.now_state]))
                f.write(" {0: 8.6f} {1: 8.6f} {2: 8.6f}".format(kin_mm, pot_mm, pot_qmmm))
                f.write(" {0: 8.6f}\n".format(kin + pot[self_tsh.now_state] + kin_mm + pot_mm + pot_qmmm))
        else:
            self.count_energy += 1

    def finalize(self,self_tsh):
        with open(self.file_restart,'w') as f: 
            f.write('Time = {0} [tau]\n'.format(str(self_tsh.elaptime)))
            f.write('QM Part\n')
            f.write('Current state: {0}\n'.format(str(self_tsh.now_state+1)))
            f.write('Coordinates (QM) [bohr]\n')
            f.write(self_tsh.mol.get_positions_formated(unit='bohr'))  
            f.write(self_tsh.mol.get_velocities_formated(unit='bohr/tau'))  
            f.write('Coefficients\n')
            for i, ic in enumerate(self_tsh.c): 
                f.write("St: {}  Re: {}  Img: {}\n".format(i, ic.real, ic.imag))   
            f.write('\nMM Part\n')
            f.write(self_tsh.mol_mm.get_positions_formated(unit='bohr'))  
            f.write(self_tsh.mol_mm.get_velocities_formated(unit='bohr/tau'))  

        time_start = self.time_start
        time_end = datetime.now()
        with open(self.file_log,'a') as f: 
            f.write("MD ends at {0:%Y-%m-%d %H:%M:%S}\n".format(time_end))
            f.write("Total MD run takes {0}\n".format(str(time_end - time_start).split('.')[0]))

