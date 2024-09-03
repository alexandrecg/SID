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

#plt.plot(A,Mi)
#plt.plot(A,Alpha2)



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#Propellant data


Rho_F = 800.0 #kg/m³

Rho_Ox = 1141.0 #kg/m³


m_F = 5.9 #g/s

m_Ox = 20.8 #g/s





Inside = "F"

Outside = "Ox"

if(Inside == "Ox" and Outside == "F"):
    
    Rho_in = Rho_Ox
    Rho_out = Rho_F

    m_in = m_Ox
    m_out = m_F
    

    
elif(Inside == "F" and Outside == "Ox"):
    
    Rho_in = Rho_F
    Rho_out = Rho_Ox
    
    m_in = m_F
    m_out = m_Ox
    
else:
    print("ERROR")


Delta_P_out = 4.0 #bar

Delta_P_in = 5.0 #bar


Alpha2_out_0 = 130 #deg

Alpha2_in_0 = 140 #deg


t_w = 1.0 #mm

delta = 0.3 #mm

Dif_min = 1 #deg



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

step_max = 1000

Dif = 0

step = 0

Alpha2_in = Alpha2_in_0

check = 0


A2_in = []

A2_out = []


Mi_in_array = []

Mi_out_array = []

while Dif < Dif_min or step < step_max:
    
    for i in range(n):
        if(abs(Alpha2[i]-Alpha2_in)<0.5):
            Mi_in = Mi[i]
            check = 1
    if(check == 0):
        print("Search failed: in")
        break
    check = 0
    
    D_n_in = np.sqrt(4/np.pi*(m_in*1e-3/(Mi_in*np.sqrt(2*Rho_in*Delta_P_in*1e5))))*1e3 #mm
    
    D_n_out = D_n_in + 2*t_w + 2*delta
    
    Mi_out = m_out*1e-3/((np.pi/4*(D_n_out*1e-3)**2)*np.sqrt(2*Rho_out*Delta_P_out*1e5))
    
    for i in range(n):
        if(abs(Mi[i]-Mi_out)<0.05):
            Alpha2_out = Alpha2[i]
            check = 1
    if(check == 0):
        print("Search failed: out")
        break
    check = 0
    
    Dif = Alpha2_in-Alpha2_out
    
    #print("Alpha2_in = %.1f deg; \tAlpha2_out = %.1f deg; \tDif = %.2f deg\n"%(Alpha2_in,Alpha2_out,Dif))
    #print("D_n_in = %.2f mm; \tD_n_out = %.2f mm\n"%(D_n_in, D_n_out))
    #print("Mi_in = %.4f; \tMi_out = %.4f\n"%(Mi_in,Mi_out))
    
    A2_in.append(Alpha2_in)
    A2_out.append(Alpha2_out)
    
    Mi_in_array.append(Mi_in)
    Mi_out_array.append(Mi_out)
        
    if(Dif < Dif_min):
        Alpha2_in = Alpha2_in - 1
    
    step += 1
    
    


fig, ax1 = plt.subplots()

ax1.plot(range(len(A2_in)),A2_in,label="in")
ax1.plot(range(len(A2_in)),A2_out,label="out")

ax1.set_ylabel("Alpha2")

ax1.legend()
ax1.grid()




fig, ax2 = plt.subplots()

ax2.plot(range(len(Mi_in_array)),Mi_in_array,label="in")
ax2.plot(range(len(Mi_out_array)),Mi_out_array,label="out")

ax2.set_ylabel("Mi")

ax2.legend()
ax2.grid()



