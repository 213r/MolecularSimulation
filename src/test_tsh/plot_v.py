import  matplotlib.pyplot as plt
import re
import numpy as np

def convert_num_ary(fname):
    ary = []  
    with open(fname) as f:
        line = f.readline()  
        if re.match("\s*#",line): line = f.readline() 
        while line: 
            ary.append(map(float, line.split()))
            line = f.readline() 
    return np.array(ary).T 

data = convert_num_ary("md_tsh_qmmm_energy.dat") 
#data = convert_num_ary("length_min_arh2.dat") 
t, ene = data[0], data[1:] 
for i in xrange(len(ene)): 
#    if i == 5: continue 
#    if 2 <= i <= 6: plt.plot(t, ene[i])
    if  i == 1: plt.plot(t, ene[i])
    if  i == 7: plt.plot(t, ene[i])
plt.xlim(0,200)
plt.show()    
