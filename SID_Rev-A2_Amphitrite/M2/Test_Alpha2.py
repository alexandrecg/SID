# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 22:26:42 2021

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

def calc_Phi(A):
        Phi_0 = 0.1 #initial value
        
        def f_Phi(x):
            return (1-x)*np.sqrt(2)/x**1.5-A
        
        Phi = opt.fsolve(f_Phi, Phi_0)
        
        if(len(Phi) == 1):
            Phi = Phi[0]
            return Phi
            
        elif(len(Phi) == 0):
            print("ERROR: no Phi value found")
            
        else:
            print("ERROR: Phi multiple values")


rad_to_deg = 180/np.pi

#injector data
N = 6
r_in = 0.495
R_in = 4.61
C = 1.1002
A_i = 13.10
Ksi = 0.9576
Lambda = 0.0275


AK = []
K = []
Phi = []
Mu = []
A_range = np.arange(0.01, 15.00, 0.001)
    
for A_i in A_range:
    Phi_i = calc_Phi(A_i)
    Phi.append(Phi_i)
    Mu_i = Phi_i*np.sqrt(Phi_i/(2-Phi_i))
    Mu.append(Mu_i)
    AK.append((1-Phi_i)*np.sqrt(2) / (Phi_i*np.sqrt(Phi_i)))
    K.append(AK[-1]/A_i)


plt.plot(A_range,K, label = "K")
#plt.plot(A_range,AK, label = "AK")
#plt.plot(A_range,Phi, label = "Phi")
#plt.plot(A_range,Mu, label = "Mu")
plt.legend()
plt.grid()


'''

Phi_i = calc_Phi(A_i)
print("Phi_i = ", Phi_i)
AK = (1-Phi_i)*np.sqrt(2) / (Phi_i*np.sqrt(Phi_i))
print("AK = ", AK)
K = AK/A_i
print("K = ",K)
'''



'''                                       
Phi_w = calc_Phi(A_i*K)
print("Phi_w = ", Phi_w)


def calc_Alpha2_i(Mu_i, A_i, Phi_i):
    return 2*rad_to_deg * np.arcsin( 2*Mu_i*A_i / (1+np.sqrt(1-Phi_i)) )

def calc_Alpha2_w(Mu_w, A_i, Phi_i, K ,Ksi, C):
    return 2*rad_to_deg * np.arcsin( 2*Mu_w*A_i*K / ((1+np.sqrt(1-Phi_i)) * np.sqrt(1-Ksi*Mu_w**2*A_i**2/C**2)) )

print("Alpha2_i = ", calc_Alpha2_i(Mu_i, A_i, Phi_i))

print("Alpha2_w = ")#, calc_Alpha2_w(Mu_w, A_i, Phi_i, K ,Ksi, C))

print("\n sin = ",  2*Mu_w*A_i*K / ((1+np.sqrt(1-Phi_w)) * np.sqrt(1-Ksi*Mu_w**2*A_i**2/C**2)) )


Alpha2_i = []
Alpha2_w = []
var = np.arange(0.0,1.0, 0.001)

for Ksi in var:
    Alpha2_i.append(calc_Alpha2_i(Mu_i, A_i, Phi_i))
    Alpha2_w.append(calc_Alpha2_w(Mu_w, A_i, Phi_i, K ,Ksi, C))
    

#plt.plot(var, Alpha2_i, label = "Alpha2_i")
#plt.plot(var, Alpha2_w, label = "Alpha2_w")
#plt.legend()
#plt.grid()

passo = 0.001
AK = []
Phi = np.arange(0.05,0.95+passo,passo)

for phi in Phi:
    AK.append((1-phi)*np.sqrt(2) / (phi*np.sqrt(phi)))

plt.plot(Phi, AK)
plt.grid()
'''