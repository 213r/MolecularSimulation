import numpy as np 
import sys 
import re

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
                while line:
                    if "Coordinates" in line:
                        line = f.readline() 
                        while line != "\n":
                            positions.append(map(float,line.split()[1:]))       
                            line = f.readline() 
                    if "Velocities" in line:
                        line = f.readline() 
                        while line != "\n":
                            velocities.append(map(float,line.split()[1:]))       
                            line = f.readline() 
                        break 
                    line = f.readline() 
                yield time, positions, velocities  
        line = f.readline() 


