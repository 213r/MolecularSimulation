import numpy as np 
import re

def read_numerical_data(file):
    data = [] 
    with open(file) as f:
        line = True 
        line = f.readline()
        while line:
            if line.strip() == "\n" or line.strip() == "": continue 
            data.append(map(float, line.split())) 
            line = f.readline()
    return np.array(data,dtype="f")

def get_statistic_info(data):
    s = "average : {}\n".format(np.average(data))
    s += "variance : {}\n".format(np.var(data))
    s += "standard deviation : {}\n".format(np.std(data))
    return s 
