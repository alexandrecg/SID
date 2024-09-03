# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 19:04:48 2021

@author: A. Goulart
"""
import numpy as np

import matplotlib.pyplot as plt


test = "test"
print("test = %s"%(test))


##### ideal injector method #####

angle_i = 80 #deg
angle_f = 140 #deg
angle_step = 0.01 #deg

angle = np.arange(angle_i,angle_f,angle_step)

phi_list = []
A_list = []
mu_list = []

for alpha2 in angle:
    
    phi = 2/(2+(np.tan(alpha2/2*np.pi/180))**2)
    phi_list.append(phi)
    
    A = (1-phi)*np.sqrt(2/phi)/phi
    A_list.append(A)
    
    mu = phi*np.sqrt(phi/(2-phi))
    mu_list.append(mu)
    

for i in range(len(angle)):
    angle[i] = angle[i]/2

fig, (ax1) = plt.subplots(2)    

ax1[0].plot(A_list, phi_list, label = "phi")
ax1[0].plot(A_list, mu_list, label = "mu")

ax1[0].set_ylabel("mu; phi")

ax1[1].plot(A_list, angle, label = "2alpha")
ax1[1].set_xlabel("A")
ax1[1].set_ylabel("alpha [deg]")

    
ax1[0].legend()
ax1[1].legend()
ax1[0].grid()
ax1[1].grid()


A_list_ideal = A_list
phi_list_ideal = phi_list

##### real injector method #####

phi_list = np.arange(0.70,0.80,0.01)

A_i = 0.50
A_f = 0.65
A_step = 0.05
A_list = np.arange(A_i, A_f + A_step, A_step)

def poly(A,phi):
    return A**2*phi**3 -2*phi**2 +4*phi -2

def calc_phi_eq(geom_char_eq):
    coefs = [geom_char_eq**2, -2, 4, -2]
    roots = np.roots(coefs)

    count = 0
    
    phi_eq = -1
    
    for root in roots:
        if(0 < root < 1 and root.imag == 0):
            count+=1
            phi_eq = root.real
            
    if(phi_eq == -1):
        print("Error in Calc_Phi_eq:\n No valid root found.")
        
    if(count > 1):
        print("Error in Calc_Phi_eq:\n Multiple valid values found (A = %.2f)"%(geom_char_eq))
        print(roots)
    
    return phi_eq




fig, (ax2) = plt.subplots(2)    

for A in A_list:
    y = []
    for phi in phi_list:
        y.append(poly(A,phi))
    
    ax2[0].plot(calc_phi_eq(A), 0, '.', markersize=10)    
    ax2[0].plot(phi_list, y, label = "A = %.2f"%(A))

ax2[0].set_xlabel("phi")
ax2[0].set_ylabel("poly(phi)")

ax2[0].legend()
ax2[0].grid()


A_i = 1
A_f = 35
A_step = 0.01
A_list = np.arange(A_i, A_f + A_step, A_step)

phi_list = []

for A in A_list:
    phi_list.append(calc_phi_eq(A))

ax2[1].plot(A_list, phi_list)
ax2[1].set_xlabel("A")
ax2[1].set_ylabel("phi")
ax2[1].grid()


##### method comparison #####

A_list = A_list_ideal
phi_list = []

error_abs = []
error_rel = []

for i in range(len(A_list)):
    
    phi_list.append(calc_phi_eq(A_list[i]))
    
    error_abs.append(phi_list_ideal[i]-phi_list[i])
    error_rel.append(abs((phi_list_ideal[i]-phi_list[i])/phi_list_ideal[i]))
    
fig, (ax3) = plt.subplots() 

ax3.plot(A_list, error_abs, label = "abs error")
ax3.plot(A_list, error_rel, label = "rel error")

ax3.set_xlabel("A")
ax3.set_ylabel("phi error")

ax3.legend()
ax3.grid()

