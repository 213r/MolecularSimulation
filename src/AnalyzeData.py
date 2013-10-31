import numpy as np 

def read_numerical_data(file):
    data = [] 
    with open(file) as f:
        line = True 
        while line:
            line = f.readline()
            if line == "\n" or line == "": continue 
            data.append(map(float, line.split())) 
    return np.array(data,dtype="f")

def get_statistic_info(data):
    s = "average : {}\n".format(np.average(data))
    s += "variance : {}\n".format(np.var(data))
    s += "standard deviation : {}\n".format(np.std(data))
    return s 
