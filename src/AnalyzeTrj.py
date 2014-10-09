import numpy as np 
import sys 
import re
from Molecule import Molecule

def get_positions_at(file, time):
    f = open(file) 
    line = f.readline()
    positions = []  
    while line:
        if "time" in line:
            if float(line.split()[2]) == time:
                for _ in xrange(2): line = f.readline()
                while line:
                    if line == "\n": return np.array(positions)
                    positions.append(map(float,line.split()[1:]))
                    line = f.readline()
        line = f.readline()

def get_velocities_at(file, time):
    f = open(file) 
    line = f.readline()
    velocities = []  
    while line:
        if "time" in line:
            if float(line.split()[2]) == time:
                while line:
                    if "velocities" in line: break 
                    line = f.readline() 
                line = f.readline()
                while line:  
                    if line == "\n": return np.array(velocities)
                    velocities.append(map(float,line.split()[1:]))
                    line = f.readline()
        line = f.readline()

def get_data(fname, start, end=float("inf")):
    f = open(fname, "r") 
    line = f.readline() 
    while line: 
        regexp = re.search(r" Time = (\d+.\d+) ",line)
        if regexp:
            time = float(regexp.group(1))  
            if start <= time and time < end: 
                positions, velocities = [], [] 
                line = f.readline() 
                natom = int(line.split()[0]) 
                line = f.readline() 
                positions = np.array([map(float,f.readline().split()[1:]) for _ in xrange(natom)])       
                line = f.readline() 
                line = f.readline() 
                velocities = np.array([map(float,f.readline().split()[1:]) for _ in xrange(natom)])       
                yield time, np.array(positions), np.array(velocities) 
        line = f.readline() 

def trj2xyz(fname, start = 0.0, end=float("inf")):
    f = open(fname, "r") 
    line = f.readline() 
    txt = "" 
    check_1st = True 
    atoms = [] 
    while line: 
        regexp = re.search(r" Time = (\d+.\d+) ",line)
        if regexp:
            time = float(regexp.group(1))  
            if start <= time and time < end: 
                positions, velocities = [], [] 
                line = f.readline() 
                natom = int(line.split()[0]) 
                line = f.readline() 
                if check_1st:
                    for i in xrange(natom):  
                        aline = f.readline().split()
                        atoms.append(aline[0]) 
                        positions.append(map(float,aline[1:]))      
                    mol = Molecule(atoms) 
                    check_1st = False 
                else:
                    positions = np.array([map(float,f.readline().split()[1:]) for _ in xrange(natom)])       
                mol.set_positions(positions) 
                txt += mol.get_positions_formated(unit="ang",message = "Times={}".format(time)) 
                line = f.readline() 
                line = f.readline() 
        line = f.readline()
    return txt


