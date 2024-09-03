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


m_F = 11.6 #5.9 #g/s

m_Ox = 23.1 #20.8 #g/s





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


Delta_P_2 = 4.0 #bar

Delta_P_1 = 6.0 #bar


Alpha2_2 = 142 #deg

Alpha2_1 = 145 #deg


t_w = 1.0 #mm

delta = 100 #% da espessura da película de fluido

recess = 0.25 #mm

angle_dif = 15 #deg

###### ######

n = 1000

'''
A_min = 0.0

A_max = 25.0

A_step = (A_max-A_min)/n

A = np.arange(A_min,A_max+A_step,A_step)
'''
Phi_min = 0.1

Phi_max = 0.90

Phi_step = (Phi_max-Phi_min)/n

Phi = np.arange(Phi_min,Phi_max+Phi_step,Phi_step)

Mi = []

A = []

Alpha2 = []

for phi in Phi:
    Mi_temp = phi*np.sqrt(phi/(2-phi))
    Alpha2_temp = 2*180/np.pi*np.arctan(np.sqrt(2*(1-phi)/phi))
    A_temp = (1-phi)/phi*np.sqrt(2/phi)
    
    Mi.append(Mi_temp)
    Alpha2.append(Alpha2_temp)
    A.append(A_temp)
    

'''
print("Mi(",len(Mi),")")
print(Mi)

print("\n\nPhi(",len(Phi),")")
print(Phi)

print("\n\nA(",len(A),")")
print(A)

fig, ax1 = plt.subplots()
ax1.plot(Phi,Mi)
ax1.set_xlabel("Phi")
ax1.set_ylabel("Mi")
ax1.grid()


fig, ax2 = plt.subplots()
ax2.plot(A,Phi)
ax2.set_xlabel("A")
ax2.set_ylabel("Phi")
ax2.grid()

fig, ax3 = plt.subplots()
ax3.plot(A,Mi)
ax3.set_xlabel("A")
ax3.set_ylabel("Mi")
ax3.grid()

fig, ax4 = plt.subplots()
ax4.plot(A,Alpha2)
ax4.set_xlabel("A")
ax4.set_ylabel("2Alpha")
ax4.grid()

'''
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### Parte 1 - Cálculo do ST1 #### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

contador_max = 10

contador = 0

Check = 0

for i in range(n):
    if(abs(Alpha2[i]-Alpha2_1)<0.1):
        Mi_1 = Mi[i]
        A_1 = A[i]
        Phi_1 = Phi[i]
        Check = 1

if(Check == 0): 
    print("ERROR: Mi_1")
    Check = 0

R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_1*1e-3/(Mi_1*np.sqrt(2*Rho_1*Delta_P_1*1e5)))

print("\nA_1 = \t%.1f; R_n_1 = \t%.3f \tmm; 2Alpha = %.1f \tdeg"%(A_1,R_n_1,Alpha2_1))


entrance = []


while contador < contador_max:
    
    n_1 = 3
    
    C = 4
    
    R_in_1 = C*R_n_1
    
    r_in_1 = np.sqrt(R_in_1*R_n_1/(n_1*A_1))
    
    entrance.append(r_in_1)
    
    k_Lin_1 = 5
    
    L_in_1 = k_Lin_1*r_in_1
    
    k_Ln_1 = 8
    
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
        if(abs(A[i]-A_eq_1)<0.05):
            Mi_eq_1 = Mi[i]
            Alpha2_eq_1 = Alpha2[i]
            Phi_eq_1 = Phi[i]
            Check = 1
    
    if(Check == 0): 
        print("ERROR: Mi_eq_1")
        Check = 0
    
    
    Mi_i_1 = Mi_eq_1/np.sqrt(1+Ksi_1*Mi_eq_1**2*A_1**2/C**2)
    
    Alpha2_eq_1_calc = 2*180/np.pi*np.arcsin(2*Mi_i_1*A_eq_1/((1+np.sqrt(1-Phi_eq_1))*np.sqrt(1-Ksi_1*Mi_i_1**2*A_1**2/C**2)))
    
    R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_1*1e-3/(Mi_i_1*np.sqrt(2*Rho_1*Delta_P_1*1e5)))
    
    A_1 = R_in_1*R_n_1/(n_1*r_in_1**2)
    
    r_mn_1 = R_n_1*np.sqrt(1-Phi_eq_1)
    
    t_fluid_1 = R_n_1-r_mn_1
    
    r_mk_1 = r_mn_1 * np.sqrt(2*(1-Phi_eq_1)/(2-Phi_eq_1))
    
    
    print("\n\nIt. %d) A_1 = \t%.2f; R_n_1 = \t%.3f \tmm; 2Alpha = %.1f (%.1f) \tdeg"%(contador+1,A_1,R_n_1,Alpha2_eq_1,Alpha2_eq_1_calc))

    contador += 1

fig, ax1 = plt.subplots()
ax1.plot(range(contador_max),entrance)
ax1.grid()

print("\nPhi = %.2f;Phi_eq = %.2f \n"%(Phi_1,Phi_eq_1))

print("r_mn_1 = %.2f mm"%(r_mn_1))
print("t_fluid_1 = %.2f mm"%(t_fluid_1))

print("\nr_mk_1 = %.2f mm"%(r_mk_1))


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
    if(abs(Alpha2[i]-Alpha2_2)<0.1):
        Mi_2 = Mi[i]
        A_2 = A[i]
        Phi_2 = Phi[i]
        Check = 1

if(Check == 0): 
    print("ERROR: Mi_2")
    Check = 0

R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_2*1e-3/(Mi_2*np.sqrt(2*Rho_2*Delta_P_2*1e5)))

print("\nA_2 = \t%.1f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f \tdeg"%(A_2,R_n_2,Alpha2_2))

entrance = []


while contador < contador_max:
    
    n_2 = 2
    
    C = 1
    
    R_in_2 = C*R_n_2
    
    r_in_2 = np.sqrt(R_in_2*R_n_2/(n_2*A_2))
    
    entrance.append(r_in_2)
    
    k_Lin_2 = 4
    
    L_in_2 = k_Lin_2*r_in_2
    
    k_Ln_2 = 1.5
    
    L_n_2 = k_Ln_2*R_n_2
    
    k_Ls_2 = 2.0
    
    L_s_2 = k_Ls_2*R_in_2
    
    R_s_2 = R_in_2+r_in_2

    Re_in_2 = (2/np.pi)*m_2*1e-3/(n_2*r_in_2*1e-3*Din_visc_2)
    
    Lambda_2 = 0.3164/(Re_in_2**0.25)
    
    Alpha_in_2 = 90 - 180/np.pi * np.arctan(R_s_2/L_in_2)
    
    
    Ksi_in_2 = -1*Alpha_in_2/150 + 1.1
    
    Ksi_2 = Ksi_in_2 + Lambda_2*L_in_2/(2*r_in_2)
    
    
    A_eq_2 = R_in_2*R_n_2/(n_2*r_in_2**2+Lambda_2/2*R_in_2*(R_in_2-R_n_2))
    
    for i in range(n):
        if(abs(A[i]-A_eq_2)<0.05):
            Mi_eq_2 = Mi[i]
            Alpha2_eq_2 = Alpha2[i]
            Phi_eq_2 = Phi[i]
            Check = 1
    
    if(Check == 0): 
        print("ERROR: Mi_eq_2")
        Check = 0
    
    #R_fluid = 
    
    Mi_i_2 = Mi_eq_2/np.sqrt(1+Ksi_2*Mi_eq_2**2*A_2**2/C**2)
    
    Alpha2_eq_2_calc = 2*180/np.pi*np.arcsin(2*Mi_i_2*A_eq_2/((1+np.sqrt(1-Phi_eq_2))*np.sqrt(1-Ksi_2*Mi_i_2**2*A_2**2/C**2)))
    
    R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(m_2*1e-3/(Mi_i_2*np.sqrt(2*Rho_2*Delta_P_2*1e5)))
    
    A_2 = R_in_2*R_n_2/(n_2*r_in_2**2)
    
    r_mn_2 = R_n_2*np.sqrt(1-Phi_eq_2)
    
    t_fluid_2 = R_n_2-r_mn_2
    
    r_mk_2 = r_mn_2 * np.sqrt(2*(1-Phi_eq_2)/(2-Phi_eq_2)) #R_n_2*np.sqrt(2*(1-Phi_eq_2)**2/(2-Phi_eq_2))
    
    
    print("\n\nIt. %d) A_2 = \t%.2f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f (%.1f) \tdeg"%(contador+1,A_2,R_n_2,Alpha2_eq_2,Alpha2_eq_2_calc))

    contador += 1

fig, ax2 = plt.subplots()
ax2.plot(range(contador_max),entrance)
ax2.grid()

print("\nPhi = %.2f;Phi_eq = %.2f \n"%(Phi_2,Phi_eq_2))

print("r_mn_2 = %.2f mm"%(r_mn_2))
print("t_fluid_2 = %.2f mm"%(t_fluid_2))

print("\nr_mk_2 = %.2f mm"%(r_mk_2))

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

if(R_n_2 - R_n_1 - t_w >= (1+delta/100)*t_fluid_2):
    print("GO! DeltaR = %.3f mm (mín: %.2f mm)"%(R_n_2 - R_n_1 - t_w,(1+delta/100)*t_fluid_2))
    
else:
    print("NOGO! DeltaR = %.3f mm (mín: %.2f mm)"%(R_n_2 - R_n_1 - t_w,(1+delta/100)*t_fluid_2))

##### #####
# Alteração do ângulo do ST2 #
##### #####

print("\n\n\t ST2 Modified Injector Geometry: \t")

Alpha2_2_goal = Alpha2_eq_1_calc - angle_dif

print("2Alpha_2_goal = \t%.1f deg"%(Alpha2_2_goal))

a_2 = 2*(1-Phi_eq_2)**2/(2-Phi_eq_2)

R_out_2 = R_n_2*np.sqrt(a_2*(1+(np.tan(Alpha2_2_goal/2*np.pi/180))**-2))

print("R_out_2 = \t\t\t%.3f mm"%(R_out_2))

L_out_2 = 0.5*R_out_2

print("L_out_2 = \t\t\t%.3f"%(L_out_2))

inj_3_H = [0,L_out_2]
inj_3_D_in = [R_out_2,R_out_2]
inj_3_D_in_mirror = [-1*R_out_2,-1*R_out_2]

y_i_3 = 0
y_f_3 = -5
y_step_3 = (y_f_3-y_i_3)/100
y_3 = np.arange(y_i_3,y_f_3+y_step_3,y_step_3)

x_3 = []
x_3_mirror = []

for y in y_3:
    x_temp = R_out_2*np.sqrt(1+(y)**2/(R_out_2*np.tan((90-0.5*Alpha2_2_goal)*np.pi/180))**2)
    x_3.append(x_temp)
    x_3_mirror.append(-1*x_temp)


##### #####
##### #####

inj_H = [0+recess,L_n_1+recess]
inj_D_in = [R_n_1,R_n_1]
inj_D_in_mirror = [-1*R_n_1,-1*R_n_1]
inj_D_out = [R_n_1+t_w,R_n_1+t_w]
inj_D_out_mirror = [-1*R_n_1-t_w,-1*R_n_1-t_w]

inj_F_H = [0+recess,0+recess]
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




fluid_D_1 = [r_mn_1,r_mn_1]
fluid_D_1_mirror = [-1*r_mn_1,-1*r_mn_1]

y_i_1 = 0+recess
y_f_1 = -5
y_step_1 = (y_f_1-y_i_1)/100
y_1 = np.arange(y_i_1,y_f_1+y_step_1,y_step_1)

x_1 = []
x_1_mirror = []

R_av_1 = (R_n_1 + r_mn_1)/2

for y in y_1:
    x_temp = R_av_1*np.sqrt(1+(y-recess)**2/(R_av_1*np.tan((90-0.5*Alpha2_eq_1_calc)*np.pi/180))**2)
    x_1.append(x_temp)
    x_1_mirror.append(-1*x_temp)


fluid_D_2 = [r_mn_2,r_mn_2]
fluid_D_2_mirror = [-1*r_mn_2,-1*r_mn_2]

y_i_2 = 0
y_f_2 = -5
y_step_2 = (y_f_2-y_i_2)/100
y_2 = np.arange(y_i_2,y_f_2+y_step_2,y_step_2)

x_2 = []
x_2_mirror = []

R_av_2 = (R_n_2 + r_mn_2)/2

print("\n")
print("1: ",Alpha2_eq_1_calc)
print("2: ",Alpha2_eq_2_calc)

for y in y_2:
    x_temp = R_av_2*np.sqrt(1+y**2/(R_av_2*np.tan((90-0.5*Alpha2_eq_2_calc)*np.pi/180))**2)
    x_2.append(x_temp)
    x_2_mirror.append(-1*x_temp)



fig, ax3 = plt.subplots()
ax3.grid()

ax3.plot(inj_D_in,inj_H, linewidth = 2, color = 'red')
ax3.plot(inj_D_in_mirror,inj_H, linewidth = 2, color = 'red')
ax3.plot(inj_D_out,inj_H, linewidth = 2, color = 'red')
ax3.plot(inj_D_out_mirror,inj_H, linewidth = 2, color = 'red')
ax3.plot(inj_F,inj_F_H, linewidth = 2, color = 'red')
ax3.plot(inj_F_mirror,inj_F_H, linewidth = 2, color = 'red')

ax3.plot(fluid_D_1,inj_H, linewidth = 1, color = 'red', linestyle = '--')
ax3.plot(fluid_D_1_mirror,inj_H, linewidth = 1, color = 'red', linestyle = '--')
ax3.plot(x_1,y_1, linewidth = 1, color = 'red', linestyle = '--')
ax3.plot(x_1_mirror,y_1, linewidth = 1, color = 'red', linestyle = '--')



ax3.plot(inj_2_D_in,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_in_mirror,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_out,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_D_out_mirror,inj_2_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_F,inj_2_F_H, linewidth = 2, color = 'b')
ax3.plot(inj_2_F_mirror,inj_2_F_H, linewidth = 2, color = 'b')

ax3.plot(fluid_D_2,inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
ax3.plot(fluid_D_2_mirror,inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
ax3.plot(x_2,y_2, linewidth = 1, color = 'blue', linestyle = '--')
ax3.plot(x_2_mirror,y_2, linewidth = 1, color = 'blue', linestyle = '--')


#ax3.plot(inj_3_D_in,inj_3_H, linewidth = 2, color = '0.3')
#ax3.plot(inj_3_D_in_mirror,inj_3_H, linewidth = 2, color = '0.3')

#ax3.plot(x_3,y_3, linewidth = 1, color = '0.3', linestyle = '--')
#ax3.plot(x_3_mirror,y_3, linewidth = 1, color = '0.3', linestyle = '--')

ax3.set_xlim(-8.0,8.0)
ax3.set_ylim(-5.0,10.0)

##### ##### ##### ##### #####
##### ##### ##### ##### #####

#### mod ####

trans = 0.5

inj_2_H = [L_out_2+trans,L_n_2]
inj_2_D_in = [R_n_2,R_n_2]
inj_2_D_in_mirror = [-1*R_n_2,-1*R_n_2]
inj_2_D_out = [R_n_2+t_w,R_n_2+t_w]
inj_2_D_out_mirror = [-1*R_n_2-t_w,-1*R_n_2-t_w]

inj_2_F_H = [0,0]
inj_2_F = [R_out_2,R_n_2+t_w]
inj_2_F_mirror = [-1*R_out_2,-1*R_n_2-t_w]


inj_3_H = [0,L_out_2-trans]
inj_3_D_in = [R_out_2,R_out_2]
inj_3_D_in_mirror = [-1*R_out_2,-1*R_out_2]

inj_23_H = [L_out_2-trans,L_out_2+trans]
inj_23_D_in = [R_out_2,R_n_2]
inj_23_D_in_mirror = [-1*R_out_2,-1*R_n_2]


##### ######


fig, ax4 = plt.subplots()
ax4.grid()

ax4.plot(inj_D_in,inj_H, linewidth = 2, color = 'red')
ax4.plot(inj_D_in_mirror,inj_H, linewidth = 2, color = 'red')
ax4.plot(inj_D_out,inj_H, linewidth = 2, color = 'red')
ax4.plot(inj_D_out_mirror,inj_H, linewidth = 2, color = 'red')
ax4.plot(inj_F,inj_F_H, linewidth = 2, color = 'red')
ax4.plot(inj_F_mirror,inj_F_H, linewidth = 2, color = 'red')

ax4.plot(fluid_D_1,inj_H, linewidth = 1, color = 'red', linestyle = '--')
ax4.plot(fluid_D_1_mirror,inj_H, linewidth = 1, color = 'red', linestyle = '--')
ax4.plot(x_1,y_1, linewidth = 1, color = 'red', linestyle = '--')
ax4.plot(x_1_mirror,y_1, linewidth = 1, color = 'red', linestyle = '--')



ax4.plot(inj_2_D_in,inj_2_H, linewidth = 2, color = '0.3')
ax4.plot(inj_2_D_in_mirror,inj_2_H, linewidth = 2, color = '0.3')
ax4.plot(inj_2_D_out,[0,L_n_2], linewidth = 2, color = '0.3')
ax4.plot(inj_2_D_out_mirror,[0,L_n_2], linewidth = 2, color = '0.3')
ax4.plot(inj_2_F,inj_2_F_H, linewidth = 2, color = '0.3')
ax4.plot(inj_2_F_mirror,inj_2_F_H, linewidth = 2, color = '0.3')

#ax4.plot(fluid_D_2,inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
#ax4.plot(fluid_D_2_mirror,inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
#ax4.plot(x_2,y_2, linewidth = 1, color = 'blue', linestyle = '--')
#ax4.plot(x_2_mirror,y_2, linewidth = 1, color = 'blue', linestyle = '--')


ax4.plot(inj_3_D_in,inj_3_H, linewidth = 2, color = '0.3')
ax4.plot(inj_3_D_in_mirror,inj_3_H, linewidth = 2, color = '0.3')

ax4.plot(x_3,y_3, linewidth = 1, color = '0.3', linestyle = '--')
ax4.plot(x_3_mirror,y_3, linewidth = 1, color = '0.3', linestyle = '--')


ax4.plot(inj_23_D_in,inj_23_H, linewidth = 2, color = '0.3')
ax4.plot(inj_23_D_in_mirror,inj_23_H, linewidth = 2, color = '0.3')

ax4.set_xlim(-8.0,8.0)
ax4.set_ylim(-5.0,10.0)