# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:57:11 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

n = 1000


A_min = 0.0

A_max = 25.0

A_step = (A_max-A_min)/n

A = np.arange(A_min,A_max+A_step,A_step)

Phi = np.arange(0.0,1.0+1/n,1/n)

Mi = []

Alpha2 = []

for a in A:
    Mi_i1 = 0
    Mi_i2 = 0
    for phi in Phi:
        Mi_i2 = 1/np.sqrt(a**2/(1-phi)+1/phi**2)
        
        if(Mi_i2>=Mi_i1):
            Mi_i1 = Mi_i2
            
            
        else:
            Mi.append(Mi_i1)
            Alpha2.append(2*np.arctan(np.sqrt(2*(1-phi)/phi))*180/np.pi)
            break


   
    
print("len(A):\t\t",len(A))
print("len(Phi):\t",len(Phi))
print("len(Mi):\t",len(Mi))
print("len(Alpha2):",len(Alpha2))

dMA = []

for i in range(n):
    if(i > 0):
        if(i == 1):
            dMA.append((Alpha2[i]-Alpha2[i-1])/(Mi[i]-Mi[i-1]))
            dMA.append((Alpha2[i]-Alpha2[i-1])/(Mi[i]-Mi[i-1]))
            dMA.append((Alpha2[i]-Alpha2[i-1])/(Mi[i]-Mi[i-1]))
            
        else:
            dMA.append((Alpha2[i]-Alpha2[i-1])/(Mi[i]-Mi[i-1]))
    
print(len(dMA))    
    
#plt.plot(A,Mi)
#plt.plot(A,Alpha2)
#plt.plot(Alpha2,Mi)
plt.plot(A,dMA)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#Propellant data

Rho_F = 800.0 #kg/m³

Rho_Ox = 1141.0 #kg/m³

m_F = 5.9 #g/s

m_Ox = 20.8 #g/s

Delta_P_F = 4.0 #bar

Delta_P_Ox = 5.0 #bar

Alpha2_F_0 = 130 #deg

Alpha2_Ox_0 = 140 #deg

t_w = 1.0 #mm

delta = 0.3 #mm

Dif_min = 1 #deg



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

step_max = 1000

Dif = 0

step = 0

Alpha2_Ox = Alpha2_Ox_0

check = 0

while Dif < Dif_min or step < step_max:
    
    for i in range(n):
        if(abs(Alpha2[i]-Alpha2_Ox)<0.5):
            Mi_Ox = Mi[i]
            check = 1
    if(check == 0):
        print("Search failed: Ox")
        break
    check = 0
    
    D_n_Ox = np.sqrt(4/np.pi*(m_Ox*1e-3/(Mi_Ox*np.sqrt(2*Rho_Ox*Delta_P_Ox*1e5))))*1e3 #mm
    
    D_n_F = D_n_Ox + 2*t_w + 2*delta
    
    Mi_F = m_F*1e-3/((np.pi/4*(D_n_F*1e-3)**2)*np.sqrt(2*Rho_F*Delta_P_F*1e5))
    
    for i in range(n):
        if(abs(Mi[i]-Mi_F)<0.05):
            Alpha2_F = Alpha2[i]
            check = 1
    if(check == 0):
        print("Search failed: F")
        break
    check = 0
    
    Dif = Alpha2_Ox-Alpha2_F
    
    print("Alpha2_Ox = %.1f deg; \tAlpha2_F = %.1f deg; \tDif = %.2f deg\n"%(Alpha2_Ox,Alpha2_F,Dif))
    print("D_n_Ox = %.2f mm; \tD_n_F = %.2f mm\n"%(D_n_Ox, D_n_F))
    print("Mi_Ox = %.4f; \tMi_F = %.4f\n"%(Mi_Ox,Mi_F))
    
    if(Dif < Dif_min):
        Alpha2_Ox = Alpha2_Ox - 1
    
    step += 1