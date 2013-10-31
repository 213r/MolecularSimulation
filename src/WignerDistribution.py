import numpy as np
from IO_MOLPRO import OutputMOLPRO, molpro_input_parser
from numpy import linalg
from Constants import au2cm,me2amu, amu2me
import sys, math

def wigner_distribution(fname_output,nsample, as_instance=False): 
    
    # read the hessian from output file 
    outp = OutputMOLPRO(fname_output)
    mol = outp.get_mol_info()
    hess_mw, q_eq = outp.get_hessian_and_eqcoord()
    q_eq = np.ravel(q_eq)

    # diagonalization
    eig, trs_mw = linalg.eigh(hess_mw)  
    eig = eig/amu2me
    hess_mw = hess_mw/amu2me
    eig = np.sqrt(eig[6:])
    #print eig*math.sqrt(me2amu)*au2cm  # print the frequency in cm-1

    # sampling mass-weighted normal-mode coordinates and momenta 
    q_rand = np.sqrt(0.5/eig)[:,np.newaxis]*np.random.randn(len(eig),nsample)
    p_rand = np.sqrt(0.5*eig)[:,np.newaxis]*np.random.randn(len(eig),nsample)
   
    # translate to mass-weighted Cartesian coordinates and momenta 
    trs_mw =  trs_mw[:,6:]
    q_mw = np.dot(trs_mw,q_rand) 
    p_mw = np.dot(trs_mw,p_rand) 

    # translate to mass-weighted Cartesian coordinates and momenta 
    #mass_3row = np.ravel(np.ones((len(mol),3))*mol.get_atomic_masses()[:,np.newaxis])
    mass_3row = np.ravel(np.ones((len(mol),3))*mol.get_masses()[:,np.newaxis])
    q = q_mw/np.sqrt(mass_3row)[:,np.newaxis] + q_eq[:,np.newaxis]
    p = p_mw*np.sqrt(mass_3row)[:,np.newaxis] 
    v = p/mass_3row[:,np.newaxis] 

    if as_instance:
        return q.T.reshape(nsample, len(mol),3),\
               q.T.reshape(nsample, len(mol),3)

    # print coordinates and velocities
    for i in xrange(nsample):
        fq_name = "coord" + str(i+1) 
        fv_name = "velocity" + str(i+1)  
        fq = open(fq_name,"w")
        fv = open(fv_name,"w")
        for j in xrange(len(mol)):
            fq.write("{0: 10.8f} {1: 10.8f} {2: 10.8f}\n".format(*q.T[i].reshape(len(mol),3)[j]))
            fv.write("{0: 10.8f} {1: 10.8f} {2: 10.8f}\n".format(*v.T[i].reshape(len(mol),3)[j]))
        fq.close()
        fv.close()

    # obtain the kintic and potential energy  
    ene_kin = 0.5*np.sum(p_mw*p_mw,axis=0)
    ene_pot = 0.5*np.sum(q_mw*np.dot(hess_mw,q_mw),axis=0)
    ene_tot = ene_kin + ene_pot
    #print 0.5*np.sum(eig)
    #print np.average(ene_kin)
    #print np.average(ene_pot)
    #print np.average(ene_pot)

def get_harmonic_frequency(wfile, mol):
    f = open(wfile,'r')
    l = f.readline()
    chk, chk2 = False, False 
    trs1,trs = [] ,[]
    freq = [] 
    while l:
        if l.find('Normal Modes') > -1: 
            chk = True
        
        if chk and l.find('Intensities [relative]') > -1: 
            for _ in xrange(len(mol)*3):
                l = f.readline()
                trs1.append(map(float,l.split()[1:]))
            if chk2: 
                trs = np.append(trs,trs1,axis=1)
            else:
                trs = trs1
                trs1 = []
                chk2 = True 
    
        if chk and l.find('Wavenumbers [cm-1]') > -1: 
            print l 
            freq.append(map(float,l.split()[2:]))
             
        if chk and l.find('Normal Modes of low/zero frequencies') > -1: 
            break 

        l = f.readline()
    
    if len(trs) == 0 or len(freq) == 0:
        print 'cant get gradient'
        sys.exit() 
    else:
        return trs.T,np.array(freq)

if __name__ == '__main__':
    mol, inp = molpro_input_parser("mlpinp0_cis.com") 
    wigner_distribution("mlpinp0_cis.out",mol,10,False)

