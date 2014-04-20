from Constants import kcal2hartree, ang2bohr 
import math
def amber_param_changer(a, b):
    return a * kcal2hartree, b * ang2bohr 


h_eps, h_sigma = amber_param_changer(0.6000, 0.0157)
c_eps, c_sigma = amber_param_changer(1.9080, 0.0860)
o_eps, o_sigma = amber_param_changer(1.6612, 0.2100)
ar_eps = 0.000379386396928
ar_sigma = 6.43640672

print math.sqrt(h_eps * ar_eps), (h_sigma + ar_sigma) / 2.0   
print math.sqrt(c_eps * ar_eps), (c_sigma + ar_sigma) / 2.0   
print math.sqrt(o_eps * ar_eps), (o_sigma + ar_sigma) / 2.0   


