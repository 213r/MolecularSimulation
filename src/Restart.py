"""
    This module is for starting the simulation using restart document.
    Each method and its return values are adapted for the kind of simulation. 

    *Variables* 
    at_start => the time when the simulation is finished. Normally, tau is used as the unit.
    positions => Coordinates  
    velocities  => Velocities   
"""

def restart_mc(fname):
    with open(fname,'r') as f: 
        positions = [] 
        line = f.readline() 
        while line:
            if "Step" in line: nstep = float(line.split()[2])
            if "Coordinates" in line:
                natom = int(line.split()[0])  
                line = f.readline() 
                for _ in xrange(natom): positions.append(map(float, f.readline().split()[1:]))
            line = f.readline() 
        return positions 

def restart_md(fname):
    with open(fname,'r') as f: 
        positions, velocities = [], [] 
        line = f.readline() 
        while line:
            if "Time" in line: at_start = float(line.split()[2])
            if "Coordinates" in line:
                natom = int(line.split()[0])  
                line = f.readline() 
                for _ in xrange(natom): positions.append(map(float, f.readline().split()[1:]))
            if "Velocities" in line:
                natom = int(line.split()[0])  
                line = f.readline() 
                for _ in xrange(natom): velocities.append(map(float, f.readline().split()[1:]))
            line = f.readline()
    return at_start, positions, velocities 

def restart_md_qmmm(fname):
    with open(fname,'r') as f: 
        positions_qm, velocities_qm = [], [] 
        positions_mm, velocities_mm = [], [] 
        line = f.readline() 
        count_pot, count_vel = True , True
        while line:
            if "Time" in line: at_start = float(line.split()[2])
            if "#Coordinates" in line:
                natom = int(line.split()[0])  
                line = f.readline()
                if count_pot:
                    for _ in xrange(natom): 
                        positions_qm.append(map(float, f.readline().split()[1:]))
                    count_pot = False
                else:
                    for _ in xrange(natom): 
                        positions_mm.append(map(float, f.readline().split()[1:]))
            if "#Velocities" in line:
                natom = int(line.split()[0])  
                line = f.readline() 
                if count_vel:
                    for _ in xrange(natom): 
                        velocities_qm.append(map(float, f.readline().split()[1:]))
                    count_vel = False
                else:
                    for _ in xrange(natom): 
                        velocities_mm.append(map(float, f.readline().split()[1:]))
            line = f.readline()
    return at_start, positions_qm, velocities_qm, positions_mm, velocities_mm,

