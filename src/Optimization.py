import numpy as np 
import sys, math 

def read_data(name_file,separator=','):
    data = [] 
    with open(name_file) as f:
        data = []  
        l = f.readline()  
        while l:  
           data_l = map(float,l.split(separator))
           data.append(data_l)
           l = f.readline() 
        return np.array(data)

def normal_equation(G,y):
    """ 
        normal_equation(G,Y) computes the closed-form solution to linear 
        regression using the normal equations.
        G => (n,m) matrix   
        y => n vector 
        n => the number of sample 
        m => the number of coefficients  
    """ 
    return  np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),np.array([y]).T)

def steepest_descent(G, y, alpha=0.01, num_iters=20000):
    m = len(y)
    theta = np.ones((len(G[0,:]),1)) * 1
    #theta[1] = - theta[1]  
    history = []  
    for i in xrange(num_iters):
        #diff_theta = - alpha * np.dot((np.dot(G, theta) - np.array([y]).T).T, G).T / m
        diff_theta = - alpha * np.dot((np.dot(G, theta) - np.array([y]).T).T, 2 * theta.T * G).T / m
        theta += diff_theta          
        #print theta  
        history.append(math.sqrt(np.sum(diff_theta ** 2))) 
    return theta, history  

def normalize(G):
    mu = np.mean(G,axis=0) 
    sigma = np.std(G,axis=0) 
    #return mu, sigma, (G - mu) / sigma 
    return sigma, G / sigma 

if __name__ == "__main__":
    data = read_data("ex1data2.txt")
    G, y = data[:,:-1],data[:,-1] 
    G = np.array(G) 
    mu, sigma, G = normalize(G) 
    G = np.c_[np.array([np.ones(len(G[:,0]))]).T,G]
    theta = np.zeros((len(G[0,:]),1))
    print mu, sigma 
    theta,history =  steepest_descent(G,y) 
    theta[1:] = theta[1:] / sigma[:,np.newaxis]   
    print theta 
    theta[0] = theta[0] - np.sum(theta[1:] * mu) 
    print theta[:]
    G, y = data[:,:-1],data[:,-1] 
    G = np.array(G) 
    G = np.c_[np.array([np.ones(len(G[:,0]))]).T,G]
    theta = normal_equation(G,y)  
    print theta  
