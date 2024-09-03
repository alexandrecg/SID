# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 07:06:17 2020

@author: A. Goulart
"""

import numpy as np

#import scipy.optimize as opt

import matplotlib.pyplot as plt



#Propellant data


Rho_F = 800.0 #kg/m³

Rho_Ox = 1141.0 #kg/m³


Din_visc_F = 0.001095 #Pa.s

Din_visc_Ox = 2.21e-3 #Pa.s


m_F = 5.9 #g/s

m_Ox = 20.8 #g/s





ST1 = "F"

ST2 = "Ox"

if(ST1 == "Ox" and ST2 == "F"):
    
    Rho_1 = Rho_Ox
    Rho_2 = Rho_F

    m_1 = m_Ox
    m_2 = m_F
        
    Din_visc_1 = Din_visc_Ox
    Din_visc_2 = Din_visc_F
    

    
elif(ST1 == "F" and ST2 == "Ox"):
    
    Rho_1 = Rho_F
    Rho_2 = Rho_Ox
    
    m_1 = m_F
    m_2 = m_Ox
        
    Din_visc_1 = Din_visc_F
    Din_visc_2 = Din_visc_Ox
    
    
else:
    print("ERROR")


Delta_P_2 = 5.0 #bar

Delta_P_1 = 5.0 #bar


Alpha2_2 = 140 #deg

Alpha2_1 = 120 #deg


t_w = 1.0 #mm

delta = 0.3 #mm

###### ######

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


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### Parte 1 - Cálculo do ST1 #### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

contador_max = 10

contador = 0

Check = 0

for i in range(n):
    if(abs(Alpha2[i]-Alpha2_1)<0.5):
        Mi_1 = Mi[i]
        A_1 = A[i]
        Check = 1

if(Check == 0): 
    print("ERROR: Mi_1")
    Check = 0

R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_1*1e-3/(Mi_1*np.sqrt(2*Rho_1*Delta_P_1*1e5)))

print("\nA_1 = \t%.1f; R_n_1 = \t%.3f \tmm; 2Alpha = %.1f \tdeg"%(A_1,R_n_1,Alpha2_1))


inlet = []

outlet = []

outlet.append(R_n_1)


while contador < contador_max:
    
    n_1 = 3
    
    C = 6
    
    R_in_1 = C*R_n_1
    
    r_in_1 = np.sqrt(R_in_1*R_n_1/(n_1*A_1))
    
    inlet.append(r_in_1)
    
    k_Lin_1 = 5
    
    L_in_1 = k_Lin_1*r_in_1
    
    k_Ln_1 = 10
    
    L_n_1 = k_Ln_1*R_n_1
    
    k_Ls_1 = 2.0
    
    L_s_1 = k_Ls_1*R_in_1
    
    R_s_1 = R_in_1+r_in_1

    Re_in_1 = (2/np.pi)*m_1*1e-3/(n_1*r_in_1*1e-3*Din_visc_1)
    
    Lambda_1 = 0.3164/(Re_in_1**0.25)
    
    Alpha_in_1 = 90 - 180/np.pi * np.arctan(R_s_1/L_in_1)
    
    
    Ksi_in_1 = -1*Alpha_in_1/150 + 1.1
    
    Ksi_1 = Ksi_in_1 + Lambda_1*L_in_1/(2*r_in_1)
    
    
    A_eq_1 = R_in_1*R_n_1/(n_1*r_in_1**2+Lambda_1/2*R_in_1*(R_in_1-R_n_1))
    
    for i in range(n):
        if(abs(A[i]-A_eq_1)<0.1):
            Mi_eq_1 = Mi[i]
            Alpha2_eq_1 = Alpha2[i]
            Phi_eq_1 = Phi[i]
            Check = 1
    
    if(Check == 0): 
        print("ERROR: Mi_eq_1")
        Check = 0
    
    
    Mi_i_1 = Mi_eq_1/np.sqrt(1+Ksi_1*Mi_eq_1**2*A_1**2/C**2)
    
    Alpha1_eq_1_calc = 1*180/np.pi*np.arcsin(1*Mi_i_1*A_eq_1/((1+np.sqrt(1-Phi_eq_1))*np.sqrt(1-Ksi_1*Mi_i_1**1*A_1**1/C**1)))
    
    R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_1*1e-3/(Mi_i_1*np.sqrt(2*Rho_1*Delta_P_1*1e5)))
    
    A_1 = R_in_1*R_n_1/(n_1*r_in_1**2)
    
    outlet.append(R_n_1)
    
    print("\n\nIt. %d) A_1 = \t%.2f; R_n_1 = \t%.3f \tmm; 2Alpha = %.1f (%.1f) \tdeg"%(contador+1,A_1,R_n_1,Alpha2_eq_1,Alpha2_eq_1))

    contador += 1

fig, ax1 = plt.subplots()

ax1.plot(range(contador_max),inlet,label="inlet")
ax1.plot(range(len(outlet)),outlet,label="outlet")

ax1.set_xlabel("Iteration")
ax1.set_ylabel("[mm]")

ax1.legend()
ax1.grid()


print("\n\n\t ST1 Injector Geometry: \t")
print("Number of Inlet Channels = \t\t%.1f"%(n_1))

print("\nRadial Dimensions:")
print("Nozzle Radius = \t\t\t%.3f \tmm"%(R_n_1))
print("Inlet Radius =  \t\t\t%.3f \tmm"%(R_in_1))
print("Swirl chamber Radius = \t\t%.3f \tmm"%(R_s_1))
print("Inlet Channel Radius = \t\t%.3f \tmm"%(r_in_1))

print("\nLinear dimensions:")
print("Nozzle Length = \t\t\t%.3f \tmm"%(L_n_1))
print("Swirl Chamber Length = \t\t%.3f \tmm"%(L_s_1))
print("Inlet Channel Length = \t\t%.3f \tmm"%(L_in_1))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### Parte 2 - Cálculo do ST2 #### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


contador_max = 10

contador = 0

Check = 0

for i in range(n):
    if(abs(Alpha2[i]-Alpha2_2)<0.5):
        Mi_2 = Mi[i]
        A_2 = A[i]
        Check = 1

if(Check == 0): 
    print("ERROR: Mi_2")
    Check = 0

R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_2*1e-3/(Mi_2*np.sqrt(2*Rho_2*Delta_P_2*1e5)))

print("\n\n\nA_2 = \t%.1f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f \tdeg"%(A_2,R_n_2,Alpha2_2))

inlet = []

outlet = []

outlet.append(R_n_2)

while contador < contador_max:
    
    n_2 = 3
    
    C = 2
    
    R_in_2 = C*R_n_2
    
    r_in_2 = np.sqrt(R_in_2*R_n_2/(n_2*A_2))
    
    inlet.append(r_in_2)
    
    k_Lin_2 = 4
    
    L_in_2 = k_Lin_2*r_in_2
    
    k_Ln_2 = 0.5
    
    L_n_2 = k_Ln_2*R_n_2
    
    k_Ls_2 = 1.0
    
    L_s_2 = k_Ls_2*R_in_2
    
    R_s_2 = R_in_2+r_in_2

    Re_in_2 = (2/np.pi)*m_2*1e-3/(n_2*r_in_2*1e-3*Din_visc_2)
    
    Lambda_2 = 0.3164/(Re_in_2**0.25)
    
    Alpha_in_2 = 90 - 180/np.pi * np.arctan(R_s_2/L_in_2)
    
    
    Ksi_in_2 = -1*Alpha_in_2/150 + 1.1
    
    Ksi_2 = Ksi_in_2 + Lambda_2*L_in_2/(2*r_in_2)
    
    
    A_eq_2 = R_in_2*R_n_2/(n_2*r_in_2**2+Lambda_2/2*R_in_2*(R_in_2-R_n_2))
    
    for i in range(n):
        if(abs(A[i]-A_eq_2)<0.1):
            Mi_eq_2 = Mi[i]
            Alpha2_eq_2 = Alpha2[i]
            Phi_eq_2 = Phi[i]
            Check = 1
    
    if(Check == 0): 
        print("ERROR: Mi_eq_2")
        Check = 0
    
    
    Mi_i_2 = Mi_eq_2/np.sqrt(1+Ksi_2*Mi_eq_2**2*A_2**2/C**2)
    
    Alpha2_eq_2_calc = 2*180/np.pi*np.arcsin(2*Mi_i_2*A_eq_2/((1+np.sqrt(1-Phi_eq_2))*np.sqrt(1-Ksi_2*Mi_i_2**2*A_2**2/C**2)))
    
    R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_2*1e-3/(Mi_i_2*np.sqrt(2*Rho_2*Delta_P_2*1e5)))
    
    A_2 = R_in_2*R_n_2/(n_2*r_in_2**2)
    
    
    outlet.append(R_n_1)
    
    print("\n\nIt. %d) A_2 = \t%.2f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f (%.1f) \tdeg"%(contador+1,A_2,R_n_2,Alpha2_eq_2,Alpha2_eq_2_calc))

    contador += 1

fig, ax2 = plt.subplots()

ax2.plot(range(contador_max),inlet,label="inlet")
ax2.plot(range(len(outlet)),outlet,label="outlet")

ax2.set_xlabel("Iteration")
ax2.set_ylabel("[mm]")

ax2.legend()
ax2.grid()


print("\n\n\t ST2 Injector Geometry: \t")
print("Number of Inlet Channels = \t\t%.1f"%(n_2))

print("\nRadial Dimensions:")
print("Nozzle Radius = \t\t\t%.3f \tmm"%(R_n_2))
print("Inlet Radius =  \t\t\t%.3f \tmm"%(R_in_2))
print("Swirl chamber Radius = \t\t%.3f \tmm"%(R_s_2))
print("Inlet Channel Radius = \t\t%.3f \tmm"%(r_in_2))

print("\nLinear dimensions:")
print("Nozzle Length = \t\t\t%.3f \tmm"%(L_n_2))
print("Swirl Chamber Length = \t\t%.3f \tmm"%(L_s_2))
print("Inlet Channel Length = \t\t%.3f \tmm"%(L_in_2))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### Parte 3 - Integração dos Estágios## ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

print("\n\n\t Geometry Check: \t")

print("D_n_1 = %.3f mm"%(2*R_n_1))

print("D_w_1 = %.3f mm"%(2*R_n_1+2*t_w))

print("D_n_2 = %.3f mm"%(2*R_n_2))

if(R_n_2 >= R_n_1 + t_w + delta):
    print("GO! DeltaR = %.3f mm"%(R_n_2 - R_n_1 - t_w))
    
else:
    print("NOGO! DeltaR = %.3f mm"%(R_n_2 - R_n_1 - t_w))

##### #####

inj_H = [0,L_n_1]
inj_D_in = [R_n_1,R_n_1]
inj_D_in_mirror = [-1*R_n_1,-1*R_n_1]
inj_D_out = [R_n_1+t_w,R_n_1+t_w]
inj_D_out_mirror = [-1*R_n_1-t_w,-1*R_n_1-t_w]

inj_F_H = [0,0]
inj_F = [R_n_1,R_n_1+t_w]
inj_F_mirror = [-1*R_n_1,-1*R_n_1-t_w]


inj_2_H = [0,L_n_2]
inj_2_D_in = [R_n_2,R_n_2]
inj_2_D_in_mirror = [-1*R_n_2,-1*R_n_2]
inj_2_D_out = [R_n_2+t_w,R_n_2+t_w]
inj_2_D_out_mirror = [-1*R_n_2-t_w,-1*R_n_2-t_w]

inj_2_F_H = [0,0]
inj_2_F = [R_n_2,R_n_2+t_w]
inj_2_F_mirror = [-1*R_n_2,-1*R_n_2-t_w]

fig, ax3 = plt.subplots()

ax3.plot(inj_D_in,inj_H, linewidth = 2, color = '0.3')
ax3.plot(inj_D_in_mirror,inj_H, linewidth = 2, color = '0.3')
ax3.plot(inj_D_out,inj_H, linewidth = 2, color = '0.3')
ax3.plot(inj_D_out_mirror,inj_H, linewidth = 2, color = '0.3')
ax3.plot(inj_F,inj_F_H, linewidth = 2, color = '0.3')
ax3.plot(inj_F_mirror,inj_F_H, linewidth = 2, color = '0.3')

ax3.plot(inj_2_D_in,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_in_mirror,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_out,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_out_mirror,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_F,inj_2_F_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_F_mirror,inj_2_F_H, linewidth = 2, color = 'b')

ax3.grid()