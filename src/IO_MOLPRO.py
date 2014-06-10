from Atom import Atom
from Molecule import Molecule
import numpy as np
import os, sys, re, pwd
from Constants import bohr2ang, ang2bohr

""" Making the molpro input 
    for simulation """ 

class InputMOLPRO:

    def __init__(self):
        # set the default parameter 
        self.memory = 400
        self.option = '' 
        self.title = 'test' 
        self.geomtype = 'xyz'
        self.basis_keyword = 'sto-3g'
        self.method = 'hf'
        self.mem_unit = 'm'
        #self.initialize()
        self.addmethod = "\n"
        self.symmetry = 'nosym'

        self.wfile = "tmp.com"

        self.user = pwd.getpwuid(os.getuid())[0]   
        self.cwd = os.getcwd() 

    def __str__(self):
        return "InputMOLPRO"

    def initialize(self):
        # delete these resulting output before starting MD 
        output = os.path.splitext(self.wfile)[0] + ".out" 
        if os.path.isfile(output):
            command = "rm {0}* ".format(output)
            os.system(command) 

    def set_molecule(self, mol):
        # set atoms as numpy array 
        self.mol = mol 

    def set_memory(self, memory, mem_unit = 'm'):
        # set memory value as integer 
        self.memory = memory
        self.mem_unit = mem_unit
    
    def set_option(self, option):
        # set option message as string 
        self.option = option

    def set_title(self, title):
        # set title message as string 
        self.title = title 

    def set_geomtype(self, geomtype):
        # set geometry type as string 
        self.geomtype = geomtype

    def set_basis_keyword(self, basis_keyword):
        self.basis_keyword = basis_keyword
    
    def set_basis_manual(self, basis_manual):
        self.basis_manual = basis_manual 

    def set_method(self, method):
        self.method = method 
    
    def add_method(self, method):
        self.addmethod += method
    
    def set_method_save(self, method):
        self.method_save = method 
    
    def set_symmetry(self, symmetry):
        self.symmetry = symmetry

    def set_inputname(self, wfile):
         self.wfile = wfile
    
    def get_inputname(self):
        return self.wfile
    
    def get_method(self):
        return self.method

    def get_method_save(self):
        return self.method_save
    
    def set_cpmcscf_istate(self,now_state,nrange, q_nacm):
        
        str_method = self.method_save
        str_method = str_method.replace("<nrange>", str(nrange)) 
        
        num = 5100.1 
        keyword_cpmcscf = "cpmcscf,grad,{0:d}.1,acc=1.d-8,record={1:4.1f}\n" \
                .format(now_state+1,num)
        keyword_force = "force;samc,{0:4.1f}\n".format(num) 
        
        if q_nacm:  
            for i in xrange(2):
                
                if i == 0 and now_state == 0: continue
                if i == 1 and now_state == nrange-1: continue
            
                if i == 0: st_i, st_j = now_state, now_state + 1 
                if i == 1: st_i, st_j = now_state+1, now_state + 2 
            
                num += 100
                keyword_cpmcscf += "cpmcscf,nacm,{0:d}.1,{1:d}.1,acc=1.d-8,record={2:4.1f}\n" \
                    .format(st_i, st_j,num)
                keyword_force += "force;samc,{0:4.1f}\n".format(num) 

        str_method = str_method.replace("<cpmcscf>", keyword_cpmcscf) 
        self.method = str_method.replace("<force>", keyword_force) 

    def set_caspt2_istate(self,now_state,nrange, q_nacm):

        str_method = self.method_save
        str_method = str_method.replace("<nrange>", str(nrange)) 
        str_method = str_method.replace("<istate>", str(now_state+1)) 
        keyword_cpmcscf, keyword_force = "", ""
        num = 5100.1 

        if q_nacm:
            for i in xrange(nrange-1):
                for j in xrange(i+1,nrange):
                    num += 100
                    keyword_cpmcscf += "cpmcscf,nacm,{0:d}.1,{1:d}.1,acc=1.d-8,record={2:4.1f}\n" \
                        .format(i+1,j+1,num)
                    keyword_force += "force;samc,{0:4.1f}\n".format(num) 
            str_method = str_method.replace("<cpmcscf>", keyword_cpmcscf) 
            str_method = str_method.replace("<force>", keyword_force) 
        self.method = str_method

    def make_input(self):
        input_txt = '' 
        input_txt += 'memory,' + str(self.memory) + \
                ',' + self.mem_unit + '\n'
        input_txt += self.option + '\n'
        #input_txt += self.symmetry + '\n'
        input_txt += 'geomtyp=' + self.geomtype + '\n'
        input_txt += 'geometry={\n' 
        #input_txt += self.symmetry + ';\n'
        input_txt += str(len(self.mol)) + '\n' 
        input_txt += self.title + '\n' 
        names = self.mol.get_atomnames() 
        positions = self.mol.get_positions() * bohr2ang
        for i in xrange(len(self.mol)):
            input_txt += "{:5s} {:18.8f} {:18.8f} {:18.8f}\n".\
                    format(names[i],*positions[i])
        input_txt += '}\n\n'

        if hasattr(self, 'basis_manual'):
            input_txt += "basis={\n" 
            input_txt += self.basis_manual
            input_txt += "}\n" 
        else:
            input_txt += 'basis=' + self.basis_keyword 

        input_txt += '\n\n'
        input_txt += self.method
        input_txt += self.addmethod
        input_txt += '\n---'

        self.addmethod = "\n"
        self.input_txt = input_txt

    def write(self):
        f = open(self.wfile, 'w')
        f.write(self.input_txt)
        f.close()

    def get_command(self):
        #submit = "/share/apps/opt/molpro/2010.1/bin/molprop_2010_1_Linux_x86_64_i8"
        submit = "/share/apps/opt/molpro/2008.1/bin/molpros_2008_1_Linux_x86_64_i8" 
        scrdir = "/scr/" + self.user 
        #scrdir = os.getcwd() 
        return "{0} --no-xml-output -d {1} -I {1} -W {1} {2}/{3}".format(submit, \
        scrdir, self.cwd, self.get_inputname())


class OutputMOLPRO:

    def __init__(self,wfile,mol=None):
        self.wfile = wfile 
        self.mol = mol 
        #if mol is None: 
        #    mol = self.get_mol_info() 
         
        #self.natom = len(mol) 
        #if mol is None: self.mol = self.get_mol_info() 
        #else: self.mol = mol 
    
    def set_outputfile(wfile):
        self.wfile = wfile

#    def get_natom(self):
#        f = open(self.wfile)
#        alist = [] 
#        cnt = 0
#        l = f.readline()
#        ele = re.compile("(\D+)")
#        while l:
#            l = f.readline()
#            if l.find("ATOMIC COORDINATES") > -1:
#                for _ in xrange(4): l = f.readline()
#                while l != '\n': 
#                    cnt += 1 
#                    l = f.readline()
#                return cnt 
#
    def get_mol_info(self):
        f = open(self.wfile)
        symbols, coords = [], []  
        l = f.readline()
        ele = re.compile("(\D+)")
        while l:
            if l.find("ATOMIC COORDINATES") > -1:
                for _ in xrange(4): 
                    l = f.readline()
                while True: 
                    atom = l.split()  
                    if atom == []: break 
                    sym, coord = ele.match(atom[1]).group() ,map(float, atom[3:6])
                    symbols.append(sym)  
                    coords.append(coord)  
                    #alist.append(Atom(symbol=sym, position=coord)) 
                    l = f.readline()
                return Molecule(atomnames = symbols, positions = coords)
            l = f.readline()

    def get_potential_energy_rhf(self):
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 

        ene = 0.0 
        l = f.readline()
        while l:
            if l.find('!RHF STATE 1.1 Energy') > -1:
                ene = float(l.split()[4])
            l = f.readline()

        try:
            return ene
        except:
            print 'no convergence'
            sys.exit() 

    def get_potential_energy_mp2(self):
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 

        ene_ary = [] 
        pattern = re.compile(r"!.*MP2") 
        l = f.readline()
        while l:
            if pattern.search(l): ene_ary.append(float(l.split()[-1]))
            l = f.readline()

        if len(ene_ary) > 1: return ene_ary 
        elif len(ene_ary) == 1: return ene_ary[0]
        else: 
            print  "No the energy at mp2 level found"
            sys.exit() 


    def get_potential_energy_ccsdt(self):
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 

        ene_ary = [] 
        pattern = re.compile(r"!.*CCSD\(T\)") 
        l = f.readline()
        while l:
            if pattern.search(l): ene_ary.append(float(l.split()[-1]))
            l = f.readline()

        if len(ene_ary) > 1: return np.array(ene_ary)
        elif len(ene_ary) == 1: return ene_ary[0]
        else: 
            print  "No the energy at ccsd(t) level found"
            sys.exit() 

#    def get_potential_energy_mcscf(self):
#        return self.get_potential_energy__mcscf_multi(nrange=1) 

    def get_potential_energy_mcscf_multi(self,nrange):
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 

        ene = np.zeros(nrange) 
        l = f.readline()
        while l:
            if l.find('!MCSCF') > -1 and l.find('.1 Energy') > -1:
                ind = int(l[l.find('.1')-1]) - 1
                ene[ind] = float(l.split()[4])
            l = f.readline()
        
        try:
            return np.array(ene)
        except:
            print 'no convergence'
            sys.exit() 

    def get_potential_energy_caspt2_nonmix(self,nrange):
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 

        ene = np.zeros(nrange) 
        l = f.readline()
        while l:
            if l.find('!RSPT2') > -1 and l.find('.1 Energy') > -1:
                ind = int(l[l.find('.1')-1]) - 1
                ene[ind] = float(l.split()[4])
            
            l = f.readline()
        
        try:
            return np.array(ene)
        except:
            print 'no convergence'
            sys.exit() 

    
    def get_data_caspt2_multi(self,nrange):
        """ 
        substract 
           1 CASPT2 ENERGY of each states  
           2 MIXING COEFFICIENTS of each states  
           3 S(uperposition) matrix   
        from molpro output file 
        """
        
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 
        
        swc = False
        ene,co,smat = [], [], []
        
        l = f.readline()
        while l:
            # getting the S matrix value  
            if "S final" in l:
                for i in xrange(2): l = f.readline()
                ind = 1 
                while ind <= nrange:   
                    lary = l.split()
                    if lary[0] == str(ind):
                        smat.append(map(float,lary[1:])) 
                    ind += 1 
                    l = f.readline()
            
            # getting the caspt2 energy and coefficients  
            if "Scale diagonal effective Hamiltonian" in l: 
                l = f.readline()
                if "Heff(scaled) = Heff(unscaled) - LS*Sii" in l:swc = True
            if swc and "MS-CASPT2 energies" in l:
                for i in xrange(3): l = f.readline()
                ind = 1 
                while ind <= nrange:   
                    lary = l.split()  
                    if lary[0] == str(ind):
                        ene.append(float(lary[1]))
                        co1 = map(float,lary[2:nrange+2])
                        co.append(co1) 
                        ind += 1
                    l = f.readline()
                break 
            l = f.readline()
        
        try:
            return np.array(ene), np.array(co), np.array(smat) 
        except:
            assert 'no convergence'
            sys.exit()

#    def get_gradient(self):
#        try: 
#            f = open(self.wfile,'r')
#        except:
#            assert 'no file of ' + self.wfile
#            sys.exit() 
#
#        f1 = [] 
#        l = f.readline()
#        while l:
#            if l.find('SCF GRADIENT FOR STATE') > -1:
#                for _ in xrange(4):l = f.readline()
#                while l != '\n': 
#                    f1.append(map(float, l.split()[1:4]))
#                    l = f.readline()
#            l = f.readline()
#        
#        if f1 == []:
#            print 'cant get forces'
#            sys.exit() 
#        else:
#            return np.array(f1) 

#    def get_gradient(self):
#        return self.get_gradient_multi(now_state=0)

    def get_gradient_multi(self,now_state):

        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 
        
        f1 = [] 
        l = f.readline()
        while l:
            if l.find('GRADIENT FOR STATE') > -1:
                if int(l[l.find('.1')-1]) == now_state + 1:
                    for _ in xrange(4):l = f.readline()
                    while l != '\n': 
                        f1.append(map(float, l.split()[1:4]))
                        l = f.readline()
            l = f.readline()
        if f1 == []:
            print 'cant get gradient'
            sys.exit() 
        else:
            return np.array(f1) 

    def get_hessian_and_eqcoord(self):
        # Note: For Mass Weighted Hessian printed in molpro output 
        # atomic mass unit(amu) is used , not atomic unit(au) 
       
        if self.mol is None: self.mol = self.get_mol_info() 
        natom = len(self.mol)  
        f = open(self.wfile,'r')
        l = f.readline()
        n = 3*natom
        hess = np.zeros((n,n)) 
        coord = [] 
        chk = False 
        while l:
            l = f.readline()
            if l.find('Mass weighted Second Derivative Matrix') > -1: 
                chk = True 
                l = f.readline()
                str = 0 
                while l != "\n": 
                    i = 0; cnt = 1 
                    while str + i < n:
                        l = f.readline()
                        hess[str+i,str:str+cnt] = map(float, l.split()[1:cnt+1])
                        if cnt <= 4: cnt += 1 
                        i += 1 
                    str += cnt  
                    l = f.readline()
           
            if chk and l.find('Atomic Coordinates') > -1: 
                for i in xrange(3): l = f.readline()
                for i in xrange(natom):
                    l = f.readline()
                    coord.append(map(float, l.split()[3:]))
        return (hess + hess.T - np.diag(np.diag(hess)), coord) 

    def get_nacme(self,nrange):
        if self.mol is None: self.mol = self.get_mol_info() 
        natom = len(self.mol)  
        try: 
            f = open(self.wfile,'r')
        except:
            print 'no file of ' + self.wfile
            sys.exit() 
        nacme = np.zeros((nrange, nrange, \
                natom, 3))
        chk = False
        f1 = [] 
        l = f.readline()
        while l:
            if l.find('NACME FOR STATES') > -1:
                chk = True 
                i, j = map(int, re.search(r"(\d).1 - (\d).1", l).groups())
                for _ in xrange(4):l = f.readline()
                while l != '\n': 
                    f1.append(map(float, l.split()[1:4]))
                    l = f.readline()
                try:  
                    nacme[i-1,j-1] = np.array(f1)
                    nacme[j-1,i-1] = -nacme[i-1,j-1]
                except:
                    print 'incosistency of dimension'
                    sys.exit() 
                f1 = []
            l = f.readline()
        
        if not chk:
            print 'cant get nacme'
            sys.exit() 
        else:
            return nacme 

"""
def molpro_input_parser(wfile):
     
    inp = InputMOLPRO() 
    f = open(wfile)
    alist = [] 
    switch_option = False
    switch_method = False
    option, method, basis = "", "", ""
    l = True 

    while l:
        l = f.readline().lower()
        
        if l.find("memory") > -1:
            l1 = l.split(',')
            inp.set_memory(int(l1[1]),l1[2]) 
            switch_option = True 
            continue
        
        if l.find("nosym") > -1:
            inp.set_symmetry(l.strip()) 
            continue 
        
        if l.find("geomtype") > -1:
            l1 = l.split('=')
            inp.set_geomtype(l1[1]) 
            continue 
                
        if l.find("geometry") > -1:
            l = f.readline().split()
            natoms = int(l[0])   
            l = f.readline()
            inp.set_title(l.strip())  
            l = f.readline()
            sym, positions = [], []
            while l[0] != '}':
                atom = l.split() 
                sym.append(atom[0])
                positions.append(map(float, atom[1:4])) 
                l = f.readline()
            if natoms != len(sym): print "Caution!! You need to check the \
                number of atoms in geometry command" 
            positions = np.array(positions) * ang2bohr
            mol = Molecule(atomnames = sym, positions = positions)
            inp.set_molecule(mol)
            continue 
        
        if l.find("basis") > -1:
            if l.find("{") > -1:  
                l = f.readline()
                while l[0] != '}':
                    basis += l
                    l = f.readline()
                inp.set_basis_manual(basis)
            else: 
                basis = l.split('=')[1].strip()
                inp.set_basis_keyword(basis)
            switch_method = True
            continue 

        if switch_method:
            if l == '---\n':
                inp.set_method(method) 
                inp.set_method_save(method) 
                break
            method += l
    #print mol.get_positions() 
    return mol, inp 
"""

def molpro_input_parser(wfile):
     
    inp = InputMOLPRO() 
    f = open(wfile)
    alist = [] 
    switch_option = False
    switch_method = False
    option, method, basis = "", "", ""
    l = True 

    while l:
        l = f.readline().lower()
        
        if l.find("memory") > -1:
            l1 = l.split(',')
            inp.set_memory(int(l1[1]),l1[2]) 
            switch_option = True 
            continue
        
        if l.find("!<option>") > -1:
            l = f.readline().lower()
            txt = "" 
            while l.find("!<>") == -1:
                txt += l 
                l = f.readline().lower()
            inp.set_option(txt) 
            continue 
 
        if l.find("geomtype") > -1:
            l1 = l.split('=')
            inp.set_geomtype(l1[1].rstrip()) 
            continue 
                
        if l.find("geometry") > -1:
            l = f.readline().split()
            natoms = int(l[0])   
            l = f.readline()
            inp.set_title(l.strip())  
            l = f.readline()
            sym, positions = [], []
            while l[0] != '}':
                atom = l.split() 
                sym.append(atom[0])
                positions.append(map(float, atom[1:4])) 
                l = f.readline()
            if natoms != len(sym): print "Caution!! You need to check the \
                number of atoms in geometry command" 
            positions = np.array(positions) * ang2bohr
            mol = Molecule(atomnames = sym, positions = positions)
            inp.set_molecule(mol)
            continue 
        
        if l.find("basis") > -1:
            if l.find("{") > -1:  
                l = f.readline()
                while l[0] != '}':
                    basis += l
                    l = f.readline()
                inp.set_basis_manual(basis)
            else: 
                basis = l.split('=')[1].strip()
                inp.set_basis_keyword(basis)
            switch_method = True
            continue 

        if switch_method:
            if l == '---\n':
                inp.set_method(method) 
                inp.set_method_save(method) 
                break
            method += l
    #print mol.get_positions() 
    return mol, inp 

def test():
    a = InputMOLPRO('test.com')
    at1 = Atom(symbol='O', position=[0.0, 0.0, 0.0]) 
    at2 = Atom(symbol='H', position=[2.0, 0.0, 0.0]) 
    at3 = Atom(symbol='H', position=[0.0, 2.0, 0.0]) 
    a.set_option('symmetry,nosym\n') 
    a.set_method('hf\nforces\n') 
    a.set_coordinate([at1,at2,at3]) 
    a.make_input() 
    a.write()

def test1():
    b = OutputMOLPRO('mlpinp0_a.out')
    print b.get_gradient(1)
    print b.get_nacme(nrange=3)

def test2():
    mol, inp = molpro_input_parser('mlpinp1.com')
    print np.array_repr(mol.get_positions(), precision=3,suppress_small=True) 
    #inp.make_input()
    #inp.write()

if __name__ == '__main__': test2()
