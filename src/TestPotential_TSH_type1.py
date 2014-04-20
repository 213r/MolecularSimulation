import numpy as np
from math import exp
#import pylab 
import sys

def nac(x):
    tan2 = 2.0 * v12(x) / (v11(x) - v22(x))  
    dadr = 2.0 * ( dv12(x) * (v11(x) - v22(x)) -  v12(x) * \
            (dv11(x) - dv22(x))) / (v11(x) -v22(x))**2
    return -dadr / (2 * (1 + tan2 * tan2))

def adiabatic(x):
    a = 0.5*(v11(x) + v22(x))
    b = 0.5*np.sqrt((v11(x) - v22(x))**2 + 4*v12(x)**2)
    return a - b, a + b

def diff_adiabatic(x):
    a = 0.5*(dv11(x) + dv22(x))
    b = 0.5*((v11(x) - v22(x)) * (dv11(x) - dv22(x)) + \
            4 * v12(x) * dv12(x)) / \
            np.sqrt((v11(x) - v22(x))**2 + 4*v12(x)**2)
    return a - b, a + b


def v11(x):
    a, b = 0.01, 1.6
    if x <= 0:
        return -a*(1.0 - exp(b*x))
    else:
        return a*(1.0 - exp(-b*x))

def dv11(x):
    a, b = 0.01, 1.6
    if x <= 0:
        return a * b * exp(b*x)
    else:
        return a * b *  exp(-b*x)

def v12(x):
    c,d = 0.005, 1.0
    return c * exp(-d*x*x)

def dv12(x):
    c,d = 0.005, 1.0
    return -2 * x *c * d * exp(-d*x*x)

def v22(x): return -v11(x)

def dv22(x): return -dv11(x)

if __name__=='__main__':

    x = np.arange(-10,10,0.1)
    vec11 = np.vectorize(v11)(x)
    vec12 = np.vectorize(v12)(x)
    dvec11 = np.vectorize(dv11)(x)
    dvec12 = np.vectorize(dv12)(x)
    vec22 = - vec11
    dvec22 = - dvec11
    
    
    tan2 = 2.0 * vec12 / (vec11 - vec22)  
    dadr = 2.0 * ( dvec12 * (vec11 - vec22) -  vec12 * (dvec11 - dvec22)) / \
            (vec11 -vec22)**2
    nac12 = -dadr / (2 * (1 + tan2 * tan2))
    
    c = np.array([1,0]) 
    v = np.array([[1,2],[3,4]])
    p = np.array([1,1])
    d = np.array([[[1,1],[2,2]],[[3,3],[4,4]]])
    print diff_nac(c,v,p,d)

    sys.exit()
    e1, e2 = adiabatic(vec11, vec22, vec12)

    #pylab.plot(x,vec11(x))
    #pylab.plot(x,-vec11(x))
    pylab.plot(x,e1)
    pylab.plot(x,e2)
    pylab.plot(x,nac/50)
    pylab.xlim([-10,10])
    pylab.ylim([-0.04,0.04])
    pylab.show()
