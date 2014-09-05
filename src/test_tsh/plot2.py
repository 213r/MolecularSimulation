import  matplotlib.pyplot as plt
import re, sys
import numpy as np
from AnalyzeTrj import get_data 
from Molecule import Molecule
from math import pi

gen = get_data("md_tsh.trj", start=0.0)

t_ary = []
x_ary = []
v_ary = []

lim = 100000
for a in gen: 
    t, x, v = a
    if t > lim: break 
    t_ary.append(t) 
    x_ary.append(x[0][0]) 
    v_ary.append(v[0][0])
plt.plot(t_ary,x_ary)
#plt.plot(t_ary,v_ary)
plt.show()
